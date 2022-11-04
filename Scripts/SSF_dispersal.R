## SSF_dispersal.R
# 
#========================================================	
# ---
### title: Dispersal integrated step selection functions
# author: Marie Gilbertson
# date: "01/03/2022"
#---
###  Preamble	
# 
# What this code does:
# 1. Uses simulated movement data to demonstrate integrated step selection functions (iSSF) for dispersing white-tailed deer

##### Clear Environment #####
remove(list=ls())


#### set seed ####
set.seed(7184)


#### load libraries ####
library(plyr)
library(adehabitatLT)
library(amt)
library(raster)
library(ggplot2)
library(tidyverse)
library(rgeos) # for road intersections
library(pbapply) # progress bar for apply functions
library(lubridate)

#### LOAD FUNCTIONS ####
## NOTE: BCRW_sim() and w.circ.mean() functions are for simulating a biased correlated random walk
## Both functions originate from Long et al 2014 (Journal of Animal Ecology) and 
## were also used in Gilbertson et al 2021 (Methods in Ecology and Evolution; associated GitHub repo available at https://github.com/mjones029/Telemetry_Network_Simulations)

BCRW_sim <- function(
  n=100,          #number of movement steps (trajectory will be subset by "minute")
  h=0.25,         #step length parameter
  rho=0,          #bias correlation parameter (0-1, where 0 -> unbiased, uncorrelated random walk, and 1 -> biased, deterministic movement)
  b=1,            #bias strength parameter
  c=0,            #bias distance decay parameter
  y0=c(0,0),       #animal starting location
  ya=c(0,0)        #animal attraction location
){
  
  #---- Main Function ------
  y <- y0
  y.t <- y
  theta.y <- runif(1,-pi,pi)       #first direction is random
  
  
  for (i in 1:n){
    
    delta <- sqrt(sum((ya-y)^2))                              #distance to attraction point 
    psi <- atan2(ya[2]-y[2],ya[1]-y[1])                       #angle toward attraction point
    beta <- tanh(b*delta^c)                                   #bias effect     
    mu <- w.circ.mean(c(theta.y,psi),c(1-beta,beta))          #biased direction
    theta.y <- rwrpnorm(1,mu,rho)                             #"draw" actual turning angle based on "expected" angle, constrained by bias correlation parameter
    #step length from chi-squared distribution
    d.y <- h*rchi(1)
    y. <- y + c(d.y*cos(theta.y),d.y*sin(theta.y))            #calculate this "step"
    
    
    #Build the trajectory
    y.t <- rbind(y.t,y.)        
    #Save for next step!
    y <- y.
  }
  
  y.out <- data.frame(y.t,row.names=NULL)
  colnames(y.out) <- c("x","y")
  
  #add date/time to trajectory; considers date/time to be on a per-minute basis
  date <- seq(1, 60*(n+1), 60)
  y.out$date <- as_datetime(date)
  return(y.out)
}

#Weighted circular mean calculation
w.circ.mean <- function (x,w) 
{
  sinr <- sum(w*sin(x)) 
  cosr <- sum(w*cos(x)) 
  circmean <- atan2(sinr, cosr) 
  circmean 
}


## function for determining if steps cross a road
road.intersections <- function(step.data, road.data){
  lines.list <- vector(mode = "list", length = nrow(step.data))
  for(i in 1:length(lines.list)){
    y <- step.data[i,]
    p1 <- c(y$x1_, y$y1_)
    p2 <- c(y$x2_, y$y2_)
    temp.line <- Line(matrix(c(p1, p2), ncol = 2, byrow = T))
    lines.list[[i]] <- Lines(list(temp.line), ID = as.character(i))
  }
  
  sp.lines <- SpatialLines(lines.list, proj4string = crs(road.data))
  

  ## check for intersections with roads
  # returns cases where there is an intersection
  intersects.roads <- gIntersects(sp.lines, road.data, byid = T, returnDense = T)
  intersects.roads2 <- colSums(intersects.roads)
  step.data$int.road <- intersects.roads2

  return(step.data)
}



#### SIMULATE MOVEMENT DATA ####

## Create a vector of time stamps to represent the time stamps for simulated animal locations
# we'll simulate 60 days of movement data, with locations taken every 4 hours (but this is an arbitrary choice)
times <- seq(Sys.time(), Sys.time()+(3*24*60*60), by = (4*60*60)) 


# for simplicity, we'll simulate JUST the dispersal phase of movement
disp <- BCRW_sim(n = length(times)-1, rho = 0.8, y0 = c(0,0), ya = c(5,5))
ggplot(disp, aes(x = x, y = y)) + geom_path() + 
  coord_fixed() + theme_bw()
disp$date <- times
disp$id <- "A"


#### SIMULATE LANDSCAPE DATA ####

## let's simulate a highway as a SpatialLinesDataFrame object
# we'll intentionally make it run alongside our dispersing individual so it appears they avoid road crossings 
# (but the individual needs to cross the road sometimes, otherwise there's not enough heterogeneity to fit a model)
# we're also going to make this a SpatialLinesDataFrame object to mimic road data that is likely publicly available
l1 <- cbind(c(-2, disp$x[1:10]), c(-0.1, disp$y[1:10]-0.05))
l2 <- cbind(c(disp$x[10:nrow(disp)], 5), c(disp$y[10:nrow(disp)]-0.05, 6))
hw <- data.frame(rbind(l1, l2)) # data frame version for easy plotting
colnames(hw) <- c("x", "y")

Sl1 <- Line(l1)
Sl2 <- Line(l2)

S1 <- Lines(list(Sl1), ID = "a")
S2 <- Lines(list(Sl2), ID = "b")

Sl <- SpatialLines(list(S1, S2))
df <- data.frame(len = sapply(1:length(Sl), function(i) gLength(Sl[i, ])))
rownames(df) <- sapply(1:length(Sl), function(i) Sl@lines[[i]]@ID)

## SpatialLines to SpatialLinesDataFrame
roads <- SpatialLinesDataFrame(Sl, data = df)


ggplot() + geom_path(data = disp, aes(x = x, y = y), colour = "blue") + 
  geom_path(data = hw, aes(x = x, y = y)) +
  coord_fixed() + theme_bw()


## let's also simulate a habitat that is a mix of agricultural and forested land
land <- raster(nrows = 16, ncols = 16, xmn = -2, xmx = 6,
               ymn = -2, ymx = 6, resolution = 0.5)
land[] <- 2
land[c(0:2), c(1:3)] <- 1
land[c(2:4), c(4:5)] <- 1
land[c(7:8), c(2:6)] <- 1
land[c(6:8), c(10:14)] <- 1
land[c(12:14), c(2:10)] <- 1

land.spdf <- as(land, "SpatialPixelsDataFrame")
land_df <- as.data.frame(land.spdf)
colnames(land_df) <- c("value", "x", "y")
land_df$class <- ifelse(land_df$value==1, "ag", "forest")

## plot the trajectory, road, and landscape together
ggplot() + geom_raster(data = land_df, aes(x = x, y = y, fill = class)) +
  scale_fill_manual(values = c("#f7fcb9", "#31a354")) +
  theme_bw() + 
  geom_path(data = disp, aes(x = x, y = y), colour = "blue") + 
  geom_path(data = hw, aes(x = x, y = y)) +
  coord_fixed()



## since our landscape is just forest or ag, let's focus on just one of these land use types in our SSF
# we'll go with agricultural land
cu <- land
cu <- calc(cu, fun = function(x){x[x == 2] <- 0; return(x)})
cu
raster::plot(cu) # our raster is just 0/1 for forest/ag (i.e., 1 = agricultural land use)


## add names to our landscape variables (works better with SSF/AMT)
names(cu) <- "cultivation"
names(roads) <- "roads"




#### PREP MOVEMENT DATA FOR SSF ####
# make track object
dt1 <- disp %>% make_track(.x = x, .y = y, .t = date, id = id)

# if we have multiple individuals, this is a great way to have a table of tables
# we just have the one individual, but we'll keep this format because you'll typically be working with many individuals
dt2 <- dt1 %>% nest(data = c(-"id")) 

# we have regular data, but especially with irregular data, we'd want to resample
dt3 <- dt2 %>% 
  mutate(resamp = lapply(data, function(x){
    x %>% amt::track_resample(rate = hours(4), tolerance = minutes(20)) %>%
      amt::filter_min_n_burst(min_n = 3)}))

# transform from points to steps
ds1 <- dt3 %>%
  mutate(stp = lapply(resamp, function(x){
    x %>% amt::steps_by_burst() 
  }))

head(ds1$stp[[1]])



#### GENERATE RANDOM STEPS, EXTRACT COVARIATES ####

## generate random steps AND extract land use for each step at the same time
ds2 <- ds1 %>%
  mutate(restp = pblapply(stp, function(x){
    x %>% amt::random_steps(n_control = 10) %>%
      amt::extract_covariates(cu, where = "both") %>%
      mutate(cult_start = factor(cultivation_start)) %>%
      mutate(cult_end = factor(cultivation_end)) %>%
      mutate(log_sl_ = log(sl_)) %>%
      mutate(cos.ang = cos(ta_))
  }))

ds2
head(ds2$restp[[1]])
table(ds2$restp[[1]]$cult_end)

## let's plot everything together
ggplot() + geom_raster(data = land_df, aes(x = x, y = y, fill = class)) +
  scale_fill_manual(values = c("#f7fcb9", "#31a354")) +
  theme_bw() + 
  geom_path(data = hw, aes(x = x, y = y)) +
  geom_segment(data = ds2$restp[[1]], aes(x = x1_, y = y1_, xend = x2_, yend = y2_, color = case_)) +
  coord_fixed()



## which steps intersect with roads
ds3 <- ds2 %>%
  mutate(restp2 = pblapply(restp, road.intersections, road.data = roads))

print(head(ds3$restp2[[1]]), width = Inf)
table(ds3$restp2[[1]]$int.road)

## plot again, but color by whether a step crosses a road
ggplot() + geom_raster(data = land_df, aes(x = x, y = y, fill = class)) +
  scale_fill_manual(values = c("#f7fcb9", "#31a354")) +
  theme_bw() + 
  geom_path(data = hw, aes(x = x, y = y)) +
  geom_segment(data = ds3$restp2[[1]], aes(x = x1_, y = y1_, xend = x2_, yend = y2_, color = as.factor(int.road))) +
  coord_fixed()

# in some cases, a step could cross multiple roads, so make sure "int.road" is just binary data
ds4 <- ds3 %>%
  mutate(restp2 = lapply(restp2, function(x){
    x %>% mutate(int.road_end = factor(ifelse(int.road>0, 1, 0)))
  }))

print(head(ds4$restp2[[1]]), width = Inf)
table(ds4$restp2[[1]]$int.road_end)



#### BUILD AND RUN INDIVIDUAL-LEVEL MODELS ####

mod1 <- ds4$restp2[[1]] %>% 
  fit_issf(case_ ~ sl_ + log_sl_ + cos.ang + 
            cult_end + 
            int.road_end +
            strata(step_id_))
summary(mod1)
## we can see that we don't have apparent selection/avoidance of cultivation, but we do seem to avoid road crossings!
