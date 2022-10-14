## Simulate_data.R
# 
#========================================================	
# ---
### title: Simulate movement data
# author: Marie Gilbertson
# date: "10/14/2022"
#---
###  Preamble	
# 
# What this code does:
# 1. Simulates animal movement data for demonstration of dispersal detection algorithm.



##### CLEAR ENVIRONMENT #####
remove(list=ls())


#### SET SEED ####
set.seed(3218)


#### LOAD LIBRARIES ####
library(adehabitatHR)
library(lubridate)
library(ggplot2)



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



#### SIMULATE TRAJECTORIES ####

## Create a vector of time stamps to represent the time stamps for simulated animal locations
# we'll simulate 60 days of movement data, with locations taken every 4 hours (but this is an arbitrary choice)
times <- seq(Sys.time(), Sys.time()+(60*24*60*60), by = (4*60*60))



#### Simulate Dispersal ####
## We'll simulate dispersal in two "phases"
## Phase 1: biased correlated random walk (BCRW) directed toward "natal" range center
## Phase 2: BCRW directed toward new range center

# set duration (in number of locations/fixes) for phase 1
phase1.n <- 150
# simulate and view phase 1 movement
disp1 <- BCRW_sim(n = phase1.n, rho = 0.8, y0 = c(5,5), ya = c(5,5))
ggplot(disp1, aes(x = x, y = y)) + geom_path() + 
  coord_fixed() + theme_bw()

# take the final location from phase 1 to use as the starting location for phase 2
end1 <- tail(disp1[,c("x", "y")],1)
# simulate and view phase 2 movement
disp2 <- BCRW_sim(n = length(times)-(phase1.n+1), rho = 0.8, y0 = c(end1$x, end1$y), ya = c(10, 10))
ggplot(disp2, aes(x = x, y = y)) + geom_path() + 
  coord_fixed() + theme_bw()

# remove last step of phase 1 trajectory so no duplicated locations when both phases are combined
disp1 <- disp1[-c(nrow(disp1)),] 
# combine phases and view
disp <- rbind(disp1, disp2)
ggplot(disp, aes(x = x, y = y)) + geom_path() + 
  coord_fixed() + theme_bw()

# add time stamps and an identifier for this individual
disp$date <- times
disp$id <- "A"


#### Simulate Residency ####
## We'll simulate a resident individual with a simple BCRW directed toward the range center

# simulate and view trajectory
res <- BCRW_sim(n = length(times)-1, rho = 0.8, y0 = c(5,5), ya = c(5,5))
ggplot(res, aes(x = x, y = y)) + geom_path() + 
  coord_fixed() + theme_bw()

# add time stamps and an identifier for this individual
res$date <- times
res$id <- "B"



#### Simulate Excursion ####
## We'll simulate an excursion in three phases:
## Phase 1: BCRW directed toward natal range
## Phase 2: Short duration BCRW directed toward an excursion target
## Phase 3: BCRW directed toward natal range

# Using the same "phase 1 duration" as dispersal movement, simulate and view phase 1 movement
ex1 <- BCRW_sim(n = phase1.n, rho = 0.8, y0 = c(5,5), ya = c(5,5))
ggplot(ex1, aes(x = x, y = y)) + geom_path() + 
  coord_fixed() + theme_bw()

# assign a short duration for the excursion behavior (in this case, 8 time steps)
ex.n <- 8
# take the final phase 1 location to use as the start of phase 2 movement
end1 <- tail(ex1[,c("x", "y")],1)
# simulate and view phase 2 movement
ex2 <- BCRW_sim(n = ex.n, rho = 0.8, y0 = c(end1$x, end1$y), ya = c(10, 10))
ggplot(ex2, aes(x = x, y = y)) + geom_path() + 
  coord_fixed() + theme_bw()

# take the final phase 2 location to use as the start of phase 3 movement
end2 <- tail(ex2[,c("x", "y")],1)
# simulate and view phase 3 movement
ex3 <- BCRW_sim(n = length(times) - (phase1.n + ex.n + 1), rho = 0.8, y0 = c(end2$x, end2$y), ya = c(5, 5))
ggplot(ex3, aes(x = x, y = y)) + geom_path() + 
  coord_fixed() + theme_bw()

# remove last step of phase 1 and 2 trajectories so no duplicated locations when all phases are combined
ex1 <- ex1[-c(nrow(ex1)),]
ex2 <- ex2[-c(nrow(ex2)),]
# combine and view
ex <- rbind(ex1, ex2, ex3)
ggplot(ex, aes(x = x, y = y)) + geom_path() + 
  coord_fixed() + theme_bw()

# add time stamps and an identifier for this individual
ex$date <- times
ex$id <- "C"



#### Simulate temporary range shift ####
## Like simulating excursions, we'll simulate a temporary range shift in three phases:
## Phase 1: BCRW directed toward natal range
## Phase 2: Longer duration BCRW directed toward an excursion target
## Phase 3: BCRW directed toward natal range

# Using the same "phase 1 duration" as dispersal movement, simulate and view phase 1 movement
rs1 <- BCRW_sim(n = phase1.n, rho = 0.8, y0 = c(5,5), ya = c(5,5))
ggplot(rs1, aes(x = x, y = y)) + geom_path() + 
  coord_fixed() + theme_bw()

# Assign a longer duration for the temporary range shift movement (in this case, 84 times steps)
rs.n <- 84
# take the final phase 1 location to use as the start of phase 2 movement
end1 <- tail(rs1[,c("x", "y")],1)
# simulate and view phase 2 movement
rs2 <- BCRW_sim(n = rs.n, rho = 0.8, y0 = c(end1$x, end1$y), ya = c(10, 10))
ggplot(rs2, aes(x = x, y = y)) + geom_path() + 
  coord_fixed() + theme_bw()

# take the final phase 2 location to use as the start of phase 3 movement
end2 <- tail(rs2[,c("x", "y")],1)
# simulate and view phase 3 movement
rs3 <- BCRW_sim(n = length(times) - (phase1.n + rs.n + 1), rho = 0.8, y0 = c(end2$x, end2$y), ya = c(5, 5))
ggplot(rs3, aes(x = x, y = y)) + geom_path() + 
  coord_fixed() + theme_bw()

# remove last step of phase 1 and 2 trajectories so no duplicated locations when all phases are combined
rs1 <- rs1[-c(nrow(rs1)),]
rs2 <- rs2[-c(nrow(rs2)),]
# combine and view
rs <- rbind(rs1, rs2, rs3)
ggplot(rs, aes(x = x, y = y)) + geom_path() + 
  coord_fixed() + theme_bw()

# add time stamps and an identifier for this individual
rs$date <- times
rs$id <- "D"


#### GENERATE SEPARATE METADATA FILE ####
meta <- data.frame(id = c("A", "B", "C", "D"),
                   collar1.on = paste(rep(times[1], 4)),
                   collar1.off = paste(rep(times[length(times)], 4)),
                   sex = c("Male", "Female", "Male", "Female"),
                   age = rep("8mo", 4)
                   )

meta


#### COMBINE AND SAVE ####
# combine all simulated individuals into one object, view, and save; this object can be used to explore the "HR_dispersal.R" script for dispersal detection 
sim.dat <- rbind(disp, res, ex, rs)

ggplot(sim.dat, aes(x = x, y = y, color = id)) + geom_path() +
  facet_grid(~id) +
  coord_fixed() + theme_bw()

save(sim.dat, file = "simulated_movement_data.Rdata")
save(meta, file = "simulated_metadata.Rdata")
