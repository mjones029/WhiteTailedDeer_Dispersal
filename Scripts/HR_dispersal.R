## HR_dispersal.R
# 
#========================================================	
# ---
### title: Home range-based dispersal analysis
# author: Marie Gilbertson
# date: "06/21/2021"
#---
###  Preamble	
# 
# What this code does:
# 1. Uses k-means clustering and home ranges to estimate dispersal from WTD collar data


##### Clear Environment #####
remove(list=ls())



#### set seed ####
set.seed(35416)


#### load libraries ####
library(adehabitatHR)
library(plyr)
library(dplyr)
library(rgeos) # centroids
library(ggplot2)
library(ggpubr)
library(bcp)
library(spatstat)
library(cluster) # for silhouette() function
library(purrr) # for map_dbl() function
library(lubridate) # for time intervals
library(raster)


#### load external functions ####

# silhouette method for determining number of clusters in k-means clustering
avg_sil <- function(k) {
  km.res <- kmeans(input, centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(input))
  mean(ss[, 3])
}



transform.clust <- function(key, old.ids){
  key.df <- key[,c("old.clust.id", "new.id")]
  key.df <- key.df[duplicated(key.df)==F,]
  new.ids <- old.ids
  for(j in 1:nrow(key.df)){
    new.ids[old.ids==key.df$old.clust.id[j]] <- key.df$new.id[j]
  }
  return(new.ids)
}

# extract final element(s) from character string
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}


# mode function 
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}



# function for verifying and updating dispersal classifications interactively
source("Scripts/verify_class.R")



#### load data ####
# load simulated metadata and movement data
meta <- get(load("simulated_metadata.Rdata"))
sim.dat <- get(load("simulated_movement_data.Rdata"))




#### convert movement data to adehabitatLT ltraj object ####
traj <- as.ltraj(sim.dat[,c("x", "y")], date = sim.dat$date, id = sim.dat$id)
traj

plot(traj)


#### DISPERSAL DETECTION LOOP START ####
full.out <- NULL # empty output for dispersal detection
seq.use.l <- vector("list", length(traj)) # empty output for sequence of cluster use
dur.use.l <- vector("list", length(traj)) # empty output for duration of cluster use

k.values <- 2:5 # number of clusters to optimize from in kmeans clustering (silhouette method)

for(i in 1:length(traj)){

  print(i)

  # to be on the safe side, null out udoi at the start of each loop
  udoi <- NULL

  #### k-means clustering of (x,y) locations ####
  
  deer1_utm <- traj[[i]][,c("x", "y")]
  input <- data.frame(x = coordinates(deer1_utm)[,1], y = coordinates(deer1_utm)[,2])

  
  ### how many clusters to evaluate? ###
  avg_sil_values <- map_dbl(k.values, avg_sil)
  
  ## optional: can plot silhouettes  
  # plot(k.values, avg_sil_values,
  #      type = "b", pch = 19, frame = FALSE,
  #      xlab = "Number of clusters K",
  #      ylab = "Average Silhouettes")

  
  # extract optimal number of clusters
  snc <- data.frame(k.values = k.values,
                    sil = avg_sil_values)
  nc <- snc$k.values[snc$sil==max(snc$sil)]
  
  
  # add date to deer1_utm
  deer1_utm$date <- traj[[i]][,c("date")]
  
  ### run kmeans clustering ###
  out <- kmeans(input, centers=nc, nstart = 25)

  input$cluster <- as.factor(out$cluster)
  centers <- data.frame(x=out$centers[,c("x")], y =out$centers[,c("y")], cluster = row.names(out$centers))
  
  
  input$time <- deer1_utm$date
  
  
  ### check clustering (optional) ###
  # plot points by central cluster
  # print(ggplot(data = input, aes(x=x, y=y, color=cluster)) + geom_point() +
  #   geom_point(data = centers, size=10, pch=21,color="black",
  #              aes(x=x, y=y, fill = cluster)))
  # 

  
  #### extract sequence of use of clusters ####
  use.seq <- input
  use.seq$cdiff <- c(NA, diff(use.seq$cluster)) # any non-zero is a cluster-use "switch"
  
  # get time intervals based on when switching occurs
  switches <- c(1, which(use.seq$cdiff != 0), nrow(use.seq)) 
  seq.out <- NULL
  
  for(j in 1:(length(switches)-1)){
    
    start <- use.seq$time[switches[j]]
    end <- use.seq$time[switches[j+1]-1]
    
    if(j==(length(switches)-1)){
      end <- tail(use.seq$time,1)
    }
    
    temp.int <- interval(start, end, tz  = "America/Chicago")
    temp.c <- use.seq$cluster[switches[j]]
    
    temp.seq.out <- data.frame(interval = temp.int,
                               old.clust.id = temp.c)
    
    seq.out <- rbind(seq.out, temp.seq.out)
  }
  
  
  # re-name clusters based on sequence of use
  seq.out$new.id <- NA
  seq.out$new.id[1] <- 1 # first cluster is always "1"
  
  for(j in 2:nrow(seq.out)){
      temp.old <- seq.out$old.clust.id[j]
      
      if(duplicated(seq.out$old.clust.id)[j]==F){
         temp.new <- max(na.omit(seq.out$new.id)) + 1
      }else{
        temp.new <- unique(na.omit(seq.out$new.id[seq.out$old.clust.id==temp.old],1))
      }
      
      seq.out$new.id[j] <- temp.new
  }
  
  # use the new id's in seq.out to update cluster id in "input" and "centers" dataframes
  input$new.id <- transform.clust(key = seq.out, old.ids = input$cluster)
  centers$new.id <- transform.clust(key = seq.out, old.ids = centers$cluster)
  
  input$cluster <- as.factor(input$new.id)

  centers$cluster <- as.factor(centers$new.id)
  centers <- centers[order(centers$cluster),]
  rownames(centers) <- seq(1, nrow(centers))
  
  
  # assemble exportable results for sequence of clusters, duration of use
  seq.vec <- paste(seq.out$new.id, collapse = "-")
  dur.vec <- as.numeric(as.duration(seq.out$interval), "days")
  dur.vec <- paste(dur.vec, collapse = "-")
  
  
  #### total duration of use for each cluster ####
  seq.out$dur.days <- as.numeric(as.duration(seq.out$interval), "days")
  
  total.dur <- data.frame(matrix(nrow = 5, ncol = 2))
  colnames(total.dur) <- c("cluster.id", "total.dur.days")
  total.dur$cluster.id <- seq(1, 5)
  
  for(j in 1:nc){
    temp.dur <- sum(seq.out$dur.days[seq.out$new.id==j])
    total.dur$total.dur.days[j] <- temp.dur
    
  }
  
  #### estimate home range vertices for identified clusters ####
  
  # make sure at least 30 locations for clusters used to estimate HRs
  c.v <- NULL
  hr.input <- input
  coverage <- ddply(hr.input, .(cluster), function(x) nrow(x))
  hr.input <- hr.input[hr.input$cluster %in% coverage$cluster[coverage$V1>=30],]
  
  c.ud <- NULL
  c.loc <- data.frame(x = hr.input$x, y = hr.input$y)
  c.sp <- SpatialPointsDataFrame(c.loc, data.frame(id = droplevels(as.factor(hr.input$cluster))))
  c.ud <- kernelUD(c.sp, h = "href", grid = 200)
  
  try(c.v <- getverticeshr(c.ud, percent = 95, unin = "m", unout = "km2"), silent = T)

  # sometimes the kernel isn't well estimated with defaults; if so, can try increasing the extent or dropping a cluster with minimal coverage (i.e. very few points)
  if(is.null(c.v)){
    repeat{
      # first, try increasing extent
      c.ud <- kernelUD(c.sp, h = "href", extent = 2)
      try(c.v <- getverticeshr(c.ud, percent = 95, unin = "m", unout = "km2"), silent = T)
      
      # if that doesn't work, remove a cluster
      if(is.null(c.v)){
        coverage2 <- coverage
        coverage2 <- coverage2[coverage2$V1!=min(coverage2$V1),]
        hr.input <- hr.input[hr.input$cluster %in% coverage2$cluster,]
        c.ud <- NULL
        c.loc <- data.frame(x = hr.input$x, y = hr.input$y)
        c.sp <- SpatialPointsDataFrame(c.loc, data.frame(id = droplevels(as.factor(hr.input$cluster))))
        c.ud <- kernelUD(c.sp, h = "href")
        try(c.v <- getverticeshr(c.ud, percent = 95, unin = "m", unout = "km2"), silent = T)
      }
      # if that still doesn't work, repeat and try again
      if(is.null(c.v)==F | nrow(coverage2)==0){
        break
      }
    }
  }
  
  ### plot home ranges with original points ###
  c.d <- suppressMessages(fortify(c.v))
  # save HR vertices for later use in estimating timing of dispersal events
  # save(c.d, file = file_name_here.Rdata))  
  
  #### plot clusters with centroids, use by time ####

  p1 <- ggplot() + geom_point(data = input, aes(x=x, y=y, color=cluster)) + geom_path(data = input, aes(x=x, y=y), alpha = 0.2)+
    geom_point(data = centers, size=10, pch=21,color="black", aes(x=x, y=y, fill = cluster)) +
    scale_color_manual(values = c("#377eb8", "#4daf4a", "#984ea3", "#e41a1c", "#ff7f00")) + 
    scale_fill_manual(values = c("#377eb8", "#4daf4a", "#984ea3", "#e41a1c", "#ff7f00")) +
    geom_polygon(data = c.d, aes(x = long, y = lat, group = group, linetype = id), colour = "black", fill = NA) +
    theme(legend.position = "none") + theme_bw()
  
  
  # plot cluster use by time
  p2 <- ggplot(data = input, aes(x=time, color = cluster, fill = cluster)) + geom_density(alpha=0.3) +
    scale_fill_manual(values = c("#377eb8", "#4daf4a", "#984ea3", "#e41a1c", "#ff7f00")) +
    scale_color_manual(values = c("#377eb8", "#4daf4a", "#984ea3", "#e41a1c", "#ff7f00")) +
    theme_bw()
  
  p <- ggarrange(p1, p2, nrow = 1)



  #### home range/cluster overlap ####
  o.df <- p.df <- matrix(nrow = 5, ncol = 5)
  rownames(p.df) <- colnames(p.df) <- rownames(o.df) <- colnames(o.df) <- seq(1, nrow(p.df))
  
  udoi <- suppressWarnings(kerneloverlap(c.sp, method = "UDOI", percent = 95, grid = 200))
  
  # if cluster(s) don't produce HR, test if those points are within other HR vertices
  no.hr.out <- NULL
  
  if(ncol(udoi)!=nc){

    # which clusters do not have a home range?
    no.hr <- centers$new.id[!centers$new.id %in% hr.input$new.id]
    
    # loop through clusters without home ranges, looping through existing home ranges
    for(k in 1:length(no.hr)){
      
      temp.nohr <- no.hr[k]
      temp.nohr.dat <- input[input$new.id==temp.nohr,]
      
      for(j in 1:length(unique(c.d$id))){
        temp.bound <- owin(poly = data.frame(x = rev(c.d$long[c.d$id==unique(c.d$id)[j]]), y = rev(c.d$lat[c.d$id==unique(c.d$id)[j]])))
        isin1 <- inside.owin(x = temp.nohr.dat$x, y = temp.nohr.dat$y, w = temp.bound)
        
        temp.out <- data.frame(hr.id = unique(c.d$id)[j],
                               nohr.id = temp.nohr,
                               isin = sum(isin1),
                               total.points = nrow(temp.nohr.dat))
        
        no.hr.out <- rbind(no.hr.out, temp.out)
      }
    }
  }
  
  # extract overlap results for later output
  for(j in 1:nrow(udoi)){
    for(k in 1:ncol(udoi)){
      ind1 <- rownames(udoi)[j]
      ind2 <- colnames(udoi)[k]
      o.df[ind1, ind2] <- udoi[j,k]
      
    }
  }
  
  # extract overlap results for cases where the points cannot produce a HR
  if(is.null(no.hr.out)==F){
    for(j in 1:nrow(p.df)){
      for(k in 1:ncol(p.df)){
        
        if(rownames(p.df)[j] %in% no.hr.out$hr.id & colnames(p.df)[k] %in% no.hr.out$nohr.id){
          temp <- no.hr.out[no.hr.out$hr.id==rownames(p.df)[j] & no.hr.out$nohr.id==colnames(p.df)[k],]
          p.df[j,k] <- p.df[k,j] <- temp$isin/temp$total.points
        }
        
        if(j==k & colnames(p.df)[k] %in% no.hr.out$nohr.id){
          p.df[j,k] <- 1 # want to assign overlap with self of 1 for screening for ranges with no overlap later
        }
        
      }
    }
  }
  
  #### estimate distance(s) between cluster centers ####
  d.df <- matrix(nrow = 5, ncol = 5)
  
  dist.centers <- as.matrix(dist(centers[,c("x", "y")], method = "euclidean")/1000)

  # extract distance results for output
  for(j in 1:nrow(dist.centers)){
    for(k in 1:ncol(dist.centers)){
      
      d.df[j,k] <- dist.centers[j,k]
      
    }
  }
  
  
  #### check if output matrices are all symmetrical ####
  if(any(c(isSymmetric(d.df), isSymmetric(o.df), isSymmetric(p.df))==F)){
    stop("assymetric matrices!")
  }
  
  #### ASSEMBLE OUTPUT ####
  temp.full.out <- data.frame(id = burst(traj[i]),
                              nc = paste(nc),
                              
                              # location of centers
                              c1.x = centers$x[1],
                              c1.y = centers$y[1],
                              c2.x = centers$x[2],
                              c2.y = centers$y[2],
                              c3.x = centers$x[3],
                              c3.y = centers$y[3],
                              c4.x = centers$x[4],
                              c4.y = centers$y[4],
                              c5.x = centers$x[5],
                              c5.y = centers$y[5],
                              
                              # distances between centers
                              dc.km.1_2 = d.df[1,2],
                              dc.km.1_3 = d.df[1,3],
                              dc.km.1_4 = d.df[1,4],
                              dc.km.1_5 = d.df[1,5],
                              dc.km.2_3 = d.df[2,3],
                              dc.km.2_4 = d.df[2,4],
                              dc.km.2_5 = d.df[2,5],
                              dc.km.3_4 = d.df[3,4],
                              dc.km.3_5 = d.df[3,5],
                              dc.km.4_5 = d.df[4,5],
                              
                              # overlap between home ranges
                              o.1_2 = o.df[1,2],
                              o.1_3 = o.df[1,3],
                              o.1_4 = o.df[1,4],
                              o.1_5 = o.df[1,5],
                              o.2_3 = o.df[2,3],
                              o.2_4 = o.df[2,4],
                              o.2_5 = o.df[2,5],
                              o.3_4 = o.df[3,4],
                              o.3_5 = o.df[3,5],
                              o.4_5 = o.df[4,5],
                              
                              # overlap of points if don't have home range
                              p.1_2 = p.df[1,2],
                              p.1_3 = p.df[1,3],
                              p.1_4 = p.df[1,4],
                              p.1_5 = p.df[1,5],
                              p.2_3 = p.df[2,3],
                              p.2_4 = p.df[2,4],
                              p.2_5 = p.df[2,5],
                              p.3_4 = p.df[3,4],
                              p.3_5 = p.df[3,5],
                              p.4_5 = p.df[4,5],
                              
                              
                              # duration of use by cluster (total use in days)
                              dur1.td = total.dur$total.dur.days[1],
                              dur2.td = total.dur$total.dur.days[2],
                              dur3.td = total.dur$total.dur.days[3],
                              dur4.td = total.dur$total.dur.days[4],
                              dur5.td = total.dur$total.dur.days[5]
                              
                              )

  # save sequence of use (vector)
  seq.use.l[[i]] <- seq.vec
  names(seq.use.l)[i] <- burst(traj[i])
  
  # save duration of use (vector)
  dur.use.l[[i]] <- dur.vec
  names(dur.use.l)[i] <- burst(traj[i])
  
  
  #### temp assign dispersal class ####
  class <- o.class <- end.class <- dur.class <- NA
  
  # if all ranges/points have some overlap: resident
  # if at least one range does not overlap: dispersal, excursion, or temporary range shift - more classification needed
  
  # if a range doesn't overlap any other ranges, it should have ov/pt = 1 (only overlap with itself)
  ov <- apply(o.df > 0, 1, sum, na.rm = T)
  pt <- apply(p.df > 0, 1, sum, na.rm = T)
  

  if(any(ov==1)| any(pt==1)){
    o.class <- "no.overlap"
  }else{
    o.class <- "overlap"
  }
  

  # if final range is different from initial range: dispersal
  # if final range is same as initial: excursion or temporary range shift - more classification needed
  init.range <- as.numeric(substr(seq.vec, 1, 1))
  final.range <- as.numeric(substrRight(seq.vec, 1))
  

  if(init.range!=final.range){
    # if first and final ranges are overlapping, consider them the "same range" for dispersal classification purposes
    o.i_f <- o.df[init.range, final.range]
    
    if(is.na(o.i_f)){
      o.i_f <- p.df[init.range, final.range]
    }
    
    if(o.i_f==0){
      end.class <- "diff.end"
    }else{
      end.class <- "same.end"
    }
  }else{
    end.class <- "same.end"
  }
  
  
  # if use of other range is short: excursion
  # if use of other range is long: temporary range shift
  
  # which range(s) don't overlap
  o.ind <- which(ov==1)
  p.ind <- which(pt==1)
  keep <- c(o.ind, p.ind)
  
  if(length(keep)>0){
    case.dur.i <- paste0("dur", keep, ".td") 
    case.dur <- temp.full.out[, case.dur.i]
  }else{
    durs <- na.omit(c(temp.full.out$dur1.td, temp.full.out$dur2.td, temp.full.out$dur3.td, temp.full.out$dur4.td, temp.full.out$dur5.td))
    case.dur <- min(durs)
  }

  
  if(any(case.dur < 7)){
   dur.class <- "short" 
  }else{
    dur.class <- "long"
  }
  
  
  
  # Now, will all the subclassifications, assign a dispersal class
  if(o.class=="overlap"){
    class <- "resident"
  }else if(o.class=="no.overlap"){
    
    if(end.class=="diff.end"){
      class <- "dispersal"
    }else if(end.class=="same.end"){
      
      if(dur.class=="short"){
        class <- "excursion"
      }else if(dur.class=="long"){
        class <- "temp.range.shift"
      }
      
    }
    
  }
  

  # add classification output to main output
  temp.class.out <- data.frame(o.class = o.class,
                               end.class = end.class,
                               dur.class = dur.class,
                               class = class
                               )

  temp.full.out <- cbind(temp.full.out, temp.class.out)
  

  
  
  # add this iteration's results to full output
  full.out <- rbind(full.out, temp.full.out)
  
  # add class and RSI estimate to plot, then save plot
  p <- annotate_figure(p, top = class)
  # save(p, file = file_name_here.Rdata)) # save the ggplot object itself for ease of making future changes
  # ggsave(file_name_here.jpeg, p, width = 16, height = 8, units = "in", dpi = 300) # save a jpeg of the plot for visualization
  
}


#### save output data ####
# save(full.out, file = file_name_here.Rdata)
# save(seq.use.l, file = file_name_here.Rdata)
# save(dur.use.l, file = file_name_here.Rdata)


#### VERIFY CLASSIFICATIONS ####

### NOTE: requires figure data to have been saved in dispersal detection algorithm
## this function will display a figure and ask the user to accept or update the dispersal determination
new.out <- verify_class(full.out,
                        dir = "HR_dispersal1")

# if skipping verification, let's just say that all classifications are correct:
full.out$new.class <- full.out$class


# save for replication/records
# save(new.out, file = "file_name_here.Rdata")




#### ESTIMATE TIMING OF DISPERSAL EVENTS ####

#### load functions ####
# function for identifying focal range/window when estimating timing of dispersal
fwind <- function(bounds = bound.list,
                  point = focus.point){
  windows <- lapply(bounds, function(z) inside.owin(x = point$x, y = point$y, w = z))
  fwind <- which(sapply(windows, '%in%', x = T))
  
  return(fwind)  
}


## look at just dispersers ##
disp.out <- full.out[full.out$new.class=="dispersal",] 
td <- traj[burst(traj) %in% disp.out$id]



full.time.out <- NULL

for(i in 1:length(td)){
  print(i)
  temp.traj <- td[i]
  nsd <- temp.traj[[1]]$R2n
  
  
  
  ## assign number of shifts to look for, based on number of centers ##
  shifts <- as.numeric(as.character(disp.out$nc[i]))-1
  
  ### Define dispersal timing as: ###
  ## Last location inside natal home range
  ## First location inside new home range
  
  
  # read in home range vertices
  c.d <- get(load(paste("File_path_", burst(td)[i], ".Rdata", sep = "")))
  
  bound.list <- vector(mode = "list", length = max(as.numeric(as.character(c.d$id))))
  for(j in 1:max(as.numeric(as.character(c.d$id)))){
    temp.id <- c.d[c.d$id==j,]
    
    # there are occasionally cases without HR vertices for a given range
    if(nrow(temp.id)==0){
      # in those cases, define a polygon based on a MCP enclosing the "best guess" of points in that cluster
      temp.meta <- full.out[full.out$id==burst(temp.traj),]
      
      # start by taking a subset of the trajectory within roughly the period in which a specific cluster was used
      # have to guesstimate by taking the total duration of use for each cluster; take +/- 5 days to get a broader window of the trajectory
      temp.dur <- temp.meta[,c("dur1.td", "dur2.td", "dur3.td", "dur4.td", "dur5.td")]
      temp.start.days <- ifelse(j==1, 0, sum(temp.dur[,c(1:(j-1))]))
      temp.end.days <- sum(temp.dur[,c(1:j)])
      temp.start <- as.POSIXct(ifelse(ceiling(temp.start.days)-5 < 0, paste(temp.meta$start), paste(temp.meta$start + days(ceiling(temp.start.days)-5))), tz = "America/Chicago")
      temp.end <- temp.meta$start + days(ceiling(temp.end.days)+5)
      temp.df <- ld(temp.traj)
      temp.df.sub <- temp.df[temp.df$date>= temp.start & temp.df$date <= temp.end,]
      
      # use kmeans to identify clusters of points within this subset of the trajectory
      temp.df.sub$cluster <- kmeans(temp.df.sub[,c("x", "y")], centers = 2, nstart = 25)$cluster
      
      
      # take the kmeans cluster for the point located closest to the original cluster center
      temp.df.sub$pt.dists <- pointDistance(temp.meta[,c(paste0("c", j, ".x"), paste0("c", j, ".y"))], temp.df.sub[,c("x", "y")], lonlat = F)
      temp.cluster <- temp.df.sub$cluster[which(temp.df.sub$pt.dists==min(temp.df.sub$pt.dists))]
      temp.df.sub.cluster <- temp.df.sub[temp.df.sub$cluster==temp.cluster,]
      
      # make a MCP for the points in this newly kmeans identified cluster of points
      coords <- SpatialPoints(temp.df.sub.cluster[,c("x", "y")], proj4string = CRS("insert_CRS_here"))
      temp.poly <- mcp(coords, unin = "m", unout = "m2")
      
      # update formatting to fit with rest of workflow
      temp.id <- fortify(temp.poly)
      temp.list <- vector(mode = "list", length = length(unique(temp.id$piece)))
      
    }else{
      temp.list <- vector(mode = "list", length = length(unique(temp.id$piece)))
    }

    for(l in 1:length(unique(temp.id$piece))){
      temp.bound <- NULL
      try(temp.bound <- owin(poly = data.frame(x = rev(temp.id$long[temp.id$piece==l]), y = rev(temp.id$lat[temp.id$piece==l]))), silent = T)
      if(is.null(temp.bound)){
        temp.bound <- owin(poly = data.frame(x = (temp.id$long[temp.id$piece==l]), y = (temp.id$lat[temp.id$piece==l])))
      }
      temp.list[[l]] <- temp.bound 
    }
    
    # take only the polygon with the largest area for a given range, unless secondary polygon(s) are at least 25% the size of the largest polygon
    areas <- unlist(lapply(temp.list, spatstat::area))
    props <- areas/max(areas)

    # take only the polygon (piece) with the highest proportion of locations for a given range (id), unless secondary polygon(s) have at least 33% of points as primary polygon
    temp.df <- ld(temp.traj)
    l.inside <- lapply(temp.list, inside.owin, x = temp.df$x, y = temp.df$y)
    n.inside <- unlist(lapply(l.inside, sum))
    p.inside <- n.inside/max(n.inside)
    
    multi.ids <- which(props >= 0.25 & p.inside/props >= (1/3))
    
    new.bound <- owin(poly = data.frame(x = rev(temp.id$long[temp.id$piece %in% multi.ids]), y = rev(temp.id$lat[temp.id$piece %in% multi.ids])))
    
    bound.list[[j]] <- new.bound
    
  }
  
  repeat{ # repeat if not happy with plotting
    ### loop through each shift and estimate timing ###
    temp.out <- NULL
    for(k in 1:shifts){
     repeat{ # repeat selection if selected points not within window
      ### interactively narrow time range for shift ###
      plot(nsd, type = "o")
      print(paste0("select timing of shift #", k, " of ", shifts))
      shift1 <- locator(2)
      
      shift1.dat <- temp.traj[[1]][floor(shift1$x[1]):ceiling(shift1$x[2]),]
      
      ### start of shift/dispersal ###
      ## which window/range is focal point within?
      focus.point <- shift1.dat[which(abs((shift1$y[1]-shift1.dat$R2n))==min(abs(shift1$y[1]-shift1.dat$R2n))),] # point at which NSD difference is minimized
      
      
      focus.window1 <- try(fwind(bounds = bound.list,
                            point = focus.point))
      
      # last point inside focal window/range?
      if(length(focus.window1)!=0){
        shift1.dat$isin1 <- inside.owin(x = shift1.dat$x, y = shift1.dat$y, w = bound.list[[focus.window1]])
        d1.s <- max(shift1.dat$date[shift1.dat$isin1==T])
      }
      
      
      ### end of shift/dispersal ###
      # which window/range is focal point within?

        focus.point <- shift1.dat[which(abs((shift1$y[2]-shift1.dat$R2n))==min(abs(shift1$y[2]-shift1.dat$R2n))),] # point at which NSD difference is minimized
        
        focus.window2 <- try(fwind(bounds = bound.list,
                              point = focus.point))
        
        # first point within focal window/range?
        if(length(focus.window2)!=0){
          
          if(length(focus.window2)>1){
            dates <- c(NULL)
            
            # find first date inside each of overlapping windows
            for(m in 1:length(focus.window2)){
              dates <- c(dates, paste(min(shift1.dat$date[inside.owin(x = shift1.dat$x, y = shift1.dat$y, w = bound.list[[focus.window2[m]]])])))
            }  
            
            # take the window with the earliest date, or if date is the same, just the lower id numbered range
            focus.window2 <- focus.window2[which(as.Date(dates)==min(as.Date(dates)))[1]]
            
            shift1.dat$isin2 <- inside.owin(x = shift1.dat$x, y = shift1.dat$y, w = bound.list[[focus.window2]])
            d1.e <- min(shift1.dat$date[shift1.dat$isin2==T & shift1.dat$isin1==F])
            
          }else if(length(focus.window2)==1){
            shift1.dat$isin2 <- inside.owin(x = shift1.dat$x, y = shift1.dat$y, w = bound.list[[focus.window2]])
            d1.e <- min(shift1.dat$date[shift1.dat$isin2==T & shift1.dat$isin1==F])
          }
        }
        
        if(length(focus.window1)!=0 & length(focus.window2)!=0 & (d1.e > d1.s))break
        print("Window fail. Select again.")
      }
      
      
      ### duration of shift/dispersal 1 ###
      dt1 <- difftime(d1.e, d1.s, tz = "America/Chicago", units = "hours")
      
      
  
      
      
      ### store results ###
      temp.temp.out <- data.frame(id = burst(temp.traj),
                             r1 = focus.window1, r2 = focus.window2,
                             d1 = d1.s, d2 = d1.e,
                             dt = dt1,
                             shift = k,
                             t.shifts = shifts)
      
      temp.out <- rbind(temp.out, temp.temp.out)
      
      # save shift1 for reproducibility
      # save(shift1, file = paste0("file_name.Rdata"))
      
      if(k < shifts){
        proceed <- readline("Proceed with another selection? Options are: yes or no ")
      }else{
        proceed <- "yes"
      }
      
      if(proceed=="no")break
    }
    ### screen shifts ###
    temp.dat <- temp.traj[[1]]
    plot.shift <- NULL
    
    for(z in 1:nrow(temp.out)){
      temp.shift <- temp.dat[temp.dat$date >= temp.out$d1[z] & temp.dat$date <= temp.out$d2[z],]
      
      plot.shift <- rbind(plot.shift, temp.shift)
    }

    plot.shift$shift <- NA
    
    for(z in 1:nrow(plot.shift)){
      for(y in 1:nrow(temp.out)){
        
        if(plot.shift$date[z] >= temp.out$d1[y] & plot.shift$date[z] <= temp.out$d2[y]){
          plot.shift$shift[z] <- temp.out$shift[y]
        }
      }
    }
    
    p <- ggplot() + geom_point(data = temp.dat, aes(x=x, y=y)) +
            geom_polygon(data = c.d, aes(x = long, y = lat, group = group, linetype = id), colour = "black", fill = NA) +
            # theme(legend.position = "none") +
            geom_point(data = plot.shift, aes(x = x, y = y, color = as.factor(shift))) +
            geom_path(data = plot.shift, aes(x = x, y = y, color = as.factor(shift)))
    
    print(p)
    pass <- readline("Happy with selection? Options are: yes or no ")
    
    if(pass=="yes")break
  }
  
  # save trajectory with timing
  # ggsave(filename = paste0("Figures/HR_dispersal/Disp_path/", burst(temp.traj), "_disperal timing path.jpeg"), p, width = 10, height = 8, units = "in", dpi = 300)
  full.time.out <- rbind(full.time.out, temp.out)

}



length(unique(full.time.out$id))


#### review and save output ####
# save(full.time.out, file = "file_name_here.Rdata")
