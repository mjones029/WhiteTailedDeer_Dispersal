## ProxRate.R
# 
#========================================================	
# ---
### title: Proximity rate
# author: Marie Gilbertson
# date: "11/30/2021"
#---
###  Preamble	
# 
# What this code does:
# 1. Uses simulated movement data to demonstrate using WildlifeDI package to calculate "Proximity" scores for individuals

##### Clear Environment #####
remove(list=ls())


#### set seed ####
set.seed(97864)


#### load libraries ####
library(adehabitatHR)
library(wildlifeDI)
library(plyr)
library(dplyr)
library(ggplot2)


#### LOAD FUNCTIONS ####
## NOTE: BCRW_sim() and w.circ.mean() functions are for simulating a biased correlated random walk
## Both functions originate from Long et al 2014 (Journal of Animal Ecology) and 
## were also used in Gilbertson et al 2021 (Methods in Ecology and Evolution; associated GitHub repo available at https://github.com/mjones029/Telemetry_Network_Simulations)

### BCRW_sim is modified here for attraction to dynamic locations!
BCRW_sim <- function(
  n=100,          #number of movement steps (trajectory will be subset by "minute")
  h=0.25,         #step length parameter
  rho=0,          #bias correlation parameter (0-1, where 0 -> unbiased, uncorrelated random walk, and 1 -> biased, deterministic movement)
  b=1,            #bias strength parameter
  c=0,            #bias distance decay parameter
  y0=c(0,0),       #animal starting location
  ya=att.locs        #animal attraction location
){
  
  #---- Main Function ------
  y <- y0
  y.t <- y
  theta.y <- runif(1,-pi,pi)       #first direction is random
  
  
  for (i in 1:n){
    
    delta <- sqrt(sum((ya[i,]-y)^2))                              #distance to attraction point 
    psi <- atan2(ya[i,2]-y[2],ya[i,1]-y[1])                       #angle toward attraction point
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






#### SIMULATE DATA ####

## we'll simulate three individuals for demonstration purposes:
# 1. a resident individual that moves independently, 
# 2. a neighboring individual whose movement is temporarily biased toward individual #1
# 3. another independent resident individual

## Create a vector of time stamps to represent the time stamps for simulated animal locations
# we'll simulate 60 days of movement data, with locations taken every 4 hours (but this is an arbitrary choice)
times <- seq(Sys.time(), Sys.time()+(60*24*60*60), by = (4*60*60))


## Let's simulate individual #1: a resident individual with a simple BCRW directed toward the range center
# make an attraction matrix; this will just be the range center for this individual
att.locs <- matrix(5, nrow = length(times), ncol = 2)
# simulate and view trajectory
ind1 <- BCRW_sim(n = length(times)-1, rho = 0.8, y0 = c(5,5), ya = att.locs)
ggplot(ind1, aes(x = x, y = y)) + geom_path() + 
  coord_fixed() + theme_bw()

# add time stamps and an identifier for this individual
ind1$date <- times
ind1$id <- "ind1"


## Let's simulate individual #2: a resident that is temporarily attracted to individual #1
# make an attraction matrix, but now the attraction point at the start and end of the trajectory will be a range center at (10,10),
# but in between, will be the current location for individual 1
att.locs <- as.matrix(ind1[,c("x", "y")])
att.locs[c(1:100, 200:length(times)),] <- 6
# simulate and view trajectory
ind2 <- BCRW_sim(n = length(times)-1, rho = 0.8, y0 = c(6,6), ya = att.locs)
ggplot(ind2, aes(x = x, y = y)) + geom_path() + 
  coord_fixed() + theme_bw()

# add time stamps and an identifier for this individual
ind2$date <- times
ind2$id <- "ind2"



## Let's simulate individual #3: another resident individual with a simple BCRW directed toward the range center
# make an attraction matrix; this will just be the range center for this individual
att.locs <- matrix(c(5,6), nrow = length(times), ncol = 2, byrow = T)
# simulate and view trajectory
ind3 <- BCRW_sim(n = length(times)-1, rho = 0.8, y0 = c(5,6), ya = att.locs)
ggplot(ind3, aes(x = x, y = y)) + geom_path() + 
  coord_fixed() + theme_bw()

# add time stamps and an identifier for this individual
ind3$date <- times
ind3$id <- "ind3"


## bind and view together
inds <- rbind(ind1, ind2, ind3)
ggplot(inds, aes(x = x, y = y, color = id)) +
  geom_path() + 
  theme_bw() + 
  coord_fixed()





#### CALCULATE PROXIMITY METRICS ####
## Use WildlifeDI package to calculate proximity (Prox() function)
# need data to be in ltraj format
t1 <- dl(ind1)
t2 <- dl(ind2)
t3 <- dl(ind3)

# we'll calculate "Prox" with a 60 minute time threshold and a 0.5 spatial units distance threshold
tp.12 <- Prox(t1, t2, tc = 60*60, dc = 0.5, local = F, GetSimultaneous = T)
tp.12


tp.13 <- Prox(t1, t3, tc = 60*60, dc = 0.5, local = F, GetSimultaneous = T)
tp.13

## if we have many neighbors, we can calculate Prox for each one to get the number of proximate (Prox>0) individuals    
    