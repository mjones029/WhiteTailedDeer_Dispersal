## InRange_cover.R
# 
#========================================================	
# ---
### title: In-range cover
# author: Marie Gilbertson
# date: "10/19/2021"
#---
###  Preamble	
# 
# What this code does:
# 1. Uses simulated movement data to demonstrate calculating proportion of cover types within home range polygons


##### Clear Environment #####
remove(list=ls())


#### set seed ####
set.seed(6846)


#### load libraries ####
library(adehabitatLT)
library(adehabitatHR)
library(raster)
library(MASS)
library(rgdal)
library(move)
library(ctmm)
library(ggplot2)
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


#### load external functions ####
source("Scripts/akde_poly.R")



#### SIMULATE MOVEMENT DATA ####
## Create a vector of time stamps to represent the time stamps for simulated animal locations
# we'll simulate 60 days of movement data, with locations taken every 4 hours (but this is an arbitrary choice)
times <- seq(Sys.time(), Sys.time()+(60*24*60*60), by = (4*60*60))


## We'll simulate a few resident individuals with a simple BCRW directed toward the range center
# we'll give these individuals coordinates in the Wisconsin Transverse Mercator system so we can assign a reasonable CRS when calculating aKDEs later

# simulate and view trajectory for individual 1
res1 <- BCRW_sim(n = length(times)-1, rho = 0.8, y0 = c(-10050000, 5325000), ya = c(-10050000, 5325000))
ggplot(res1, aes(x = x, y = y)) + geom_path() + 
  coord_fixed() + theme_bw()

# add time stamps and an identifier for this individual
res1$date <- times
res1$id <- "A"



# simulate and view trajectory for individual 1
res2 <- BCRW_sim(n = length(times)-1, rho = 0.8, y0 = c(-10050002, 5325002), ya = c(-10050002, 5325002))
ggplot(res2, aes(x = x, y = y)) + geom_path() + 
  coord_fixed() + theme_bw()

# add time stamps and an identifier for this individual
res2$date <- times
res2$id <- "B"



# simulate and view trajectory for individual 1
res3 <- BCRW_sim(n = length(times)-1, rho = 0.8, y0 = c(-10050000, 5325003), ya = c(-10050000, 5325003))
ggplot(res3, aes(x = x, y = y)) + geom_path() + 
  coord_fixed() + theme_bw()

# add time stamps and an identifier for this individual
res3$date <- times
res3$id <- "C"


## bind and view together
res <- rbind(res1, res2, res3)
ggplot(res, aes(x = x, y = y, colour = id)) + geom_path() + 
  coord_fixed() + theme_bw()


#### SIMULATE LANDSCAPE DATA ####
## let's simulate a habitat that is a mix of agricultural and forested land
land <- raster(nrows = 10, ncols = 10, xmn = -10050007, xmx = -10049997,
               ymn = 5324995, ymx = 5325005, resolution = 1)
land[] <- 2
land[c(0:2), c(1:3)] <- 1
land[c(4:5), c(5:8)] <- 1
land[c(9:10), c(1:7)] <- 1


land.spdf <- as(land, "SpatialPixelsDataFrame")
land_df <- as.data.frame(land.spdf)
colnames(land_df) <- c("value", "x", "y")
land_df$class <- ifelse(land_df$value==1, "ag", "forest")


## plot the trajectory, road, and landscape together
ggplot() + geom_raster(data = land_df, aes(x = x, y = y, fill = class)) +
  scale_fill_manual(values = c("#f7fcb9", "#31a354")) +
  theme_bw() + 
  geom_path(data = res, aes(x = x, y = y, colour = id)) + 
  coord_fixed()



#### CALCULATE PROPORTION OF LAND TYPE IN PRE-EXPECTED DISPERSAL RANGE ####

## list of individuals to evaluate
inds <- unique(res$id)


## empty object for output
rc.out <- data.frame(matrix(nrow = length(inds), ncol = 4))
colnames(rc.out) <- c("id", "area", "ag", "forest")
rc.out$id <- inds

for(i in 1:length(inds)){
  print(i)
  
  # data for one individual
  tempd.df <- res[res$id==inds[i],]
  
  
  # use CTMM to calculate aKDE
  aout <- akde_poly(tempd.df)
  new.area <- aout$area
  vert <- aout$pl.est
    
  ## get land cover within range polygon (vert)
  e <- suppressWarnings(extract(land, vert)) 
  class.counts <- lapply(e, table)
  class.prop <- data.frame(lapply(e, FUN = function(x) { prop.table(table(x)) }))
  
  # can update below if, for example, multiple land use classes would be considered "ag"
  ag <- ifelse(any(c(1) %in% class.prop$x), sum(class.prop$Freq[class.prop$x %in% c(1)]), 0)
  forest <- ifelse(any(c(2) %in% class.prop$x), sum(class.prop$Freq[class.prop$x %in% c(2)]), 0)


  # plot (optional)
  suppressMessages(vert.df <- fortify(vert))
  print(ggplot() + geom_raster(data = land_df, aes(x = x, y = y, fill = class)) +
    scale_fill_manual(values = c("#f7fcb9", "#31a354")) +
    geom_point(data = tempd.df, aes(x = x, y = y), alpha = 0.2) +
    geom_polygon(data = vert.df, aes(x = long, y = lat, group = group), color = "black", fill = "transparent") +
    theme_bw() + 
    coord_fixed()
  )
  
  #### save back key data ####
  rc.out$area[rc.out$id==inds[i]] <- new.area
  rc.out$ag[rc.out$id==inds[i]] <- ag
  rc.out$forest[rc.out$id==inds[i]] <- forest
  
}

rc.out

### we can store back information about land cover in home ranges and use them in statistical modeling of dispersal probabilities and distances

