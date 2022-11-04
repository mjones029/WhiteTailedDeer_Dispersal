## Traversability.R
# 
#========================================================	
# ---
### title: Dispersal traversability
# author: Marie Gilbertson
# date: "10/25/21"
#---
###  Preamble	
# 
# What this code does:
# 1. Estimates "traversability" of environment in and around pre-dispersal range using simulated potential dispersal paths
# Note that simulated dispersal paths are based on hidden Markov models

##### Clear Environment #####
remove(list=ls())


#### set seed ####
set.seed(2651)

#### load libraries ####
library(momentuHMM)
library(adehabitatLT) # for dispersing deer example
library(ggplot2)
library(plyr)
library(dplyr)
library(raster)
library(rgdal)
library(rgeos)
library(sp)
library(raster)
library(amt)


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


#### SIMULATE DATA ####
## we'll simulate some new dispersal data so we have two distinct "movement phases" for hidden Markov modeling (HMM):
# (1) a residency state with short step lengths and low concentration of turning angles, and 
# (2) a dispersal state with long step lengths and more concentrated turning angles

# set duration (in number of locations/fixes) for phase 1
times <- seq(Sys.time(), Sys.time()+(60*24*60*60), by = (4*60*60))
phase1.n <- 150

# simulate and view phase 1 movement
disp1 <- BCRW_sim(n = phase1.n, h = 0.25, rho = 0.8, y0 = c(5,5), ya = c(5,5))
ggplot(disp1, aes(x = x, y = y)) + geom_path() + 
  coord_fixed() + theme_bw()

# take the final location from phase 1 to use as the starting location for phase 2
end1 <- tail(disp1[,c("x", "y")],1)

# assign a short duration for the dispersal behavior (in this case, 10 time steps)
ex.n <- 10
# simulate and view phase 2 movement
disp2 <- BCRW_sim(n = ex.n, h = 2, rho = 0.8, y0 = c(end1$x, end1$y), ya = c(15, 15))
ggplot(disp2, aes(x = x, y = y)) + geom_path() + 
  coord_fixed() + theme_bw()

# take the final phase 2 location to use as the start of phase 3 movement
end2 <- tail(disp2[,c("x", "y")],1)
# simulate and view phase 3 movement
disp3 <- BCRW_sim(n = length(times) - (phase1.n + ex.n + 1), h = 0.25, rho = 0.8, y0 = c(end2$x, end2$y), ya = c(15, 15))
ggplot(disp3, aes(x = x, y = y)) + geom_path() + 
  coord_fixed() + theme_bw()


disp1 <- disp1[-c(nrow(disp1)),]
disp2 <- disp2[-c(nrow(disp2)),]
# combine and view
disp <- rbind(disp1, disp2, disp3)
ggplot(disp, aes(x = x, y = y)) + geom_path() + 
  coord_fixed() + theme_bw()


# add time stamps and an identifier for this individual
disp$date <- times
disp$id <- "A"






#### PREP DATA for HMM ####

## Now let's prep this data for HMM
traj.df <- disp

## for illustrative purposes, let's drop a handful of locations randomly from the trajectory to mimic the fact that GPS collar data typically has missing locations
traj.df <- traj.df[-sample(1:nrow(traj.df), size = 25, replace = F),]
ggplot(traj.df, aes(x = x, y = y, colour = id)) + geom_path() + coord_fixed() + theme_bw()


### regularization: interpolate location gaps with CRW ###
## check how often our locations are typically recorded
table(as.numeric(difftime(traj.df$date[2:nrow(traj.df)], traj.df$date[1:(nrow(traj.df)-1)], units = "hours")))

# locations are typically recorded every 4 hours, so we'll use crawlWrap() to fit CTCRW and predict locations every 4 hours
colnames(traj.df) <- c("x", "y", "time", "ID") # required column names for crawlWrap
# crawl expects time stamps to be rounded to the top of the hour
traj.df$time <- as.POSIXct(round(traj.df$time,"hours"), tz = "America/Chicago")
d.crw <- crawlWrap(obsData=traj.df, timeStep="4 hours",
                   theta=c(8, 2), fixPar=c(NA,NA))

# make a fresh momentuHMM object with the interpolated, fully regularized data
cd <- prepData(data = d.crw$crwPredict[,c("ID", "time", "mu.x", "mu.y")], type = "UTM", coordNames = c("mu.x", "mu.y"))
head(cd)
ggplot(cd, aes(x = x, y = y)) + geom_path() + coord_fixed() + theme_bw()


#### FIT HMMs ####

#### fit super simple two-state model with no covariates (correlated random walk?) ####

### build and fit the model ###
# label states
stateNames <- c("encamped","dispersing")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# initial parameters
Par0_m1 <- list(step=c(0.3,0.1,2,0.1),angle=c(3, 0, 0.5, 0.8)) # mean, sd of gamma for each state (mean, mean, sd, sd); mean and concentration of wrapped Cauchy for each state
# fit model
m1 <- fitHMM(data = cd, nbStates = 2, dist = dist, Par0 = Par0_m1,
             estAngleMean = list(angle=TRUE), stateNames = stateNames)
m1

### view/check model outputs ###
## For a given model, the function viterbi() computes the most likely state sequence
states <- viterbi(m1)
# derive percentage of time spent in each state 
print(table(states)/nrow(cd))


print(ggplot(cd, aes(x = x, y = y, color = factor(states), group = 1)) +
  geom_path() +
  scale_color_manual(values = c("#fec44f", "#4eb3d3"))) +
  theme_bw()


#### fit super simple three-state model with no covariates ####
### build and fit the model ###
# update label states
stateNames <- c("encamped", "intermediate", "dispersing")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# initial parameters
Par0_m2 <- list(step=c(0.1, 0.5, 2, 0.1, 1, 1), # mean, sd of gamma for each state (mean, mean, mean, sd, sd, sd);
                angle=c(3, 3, 0, 0.5, 0.5, 0.8))  # mean and concentration of wrapped Cauchy for each state (mean, mean, mean, conc, conc, conc)
# fit model
m2 <- fitHMM(data = cd, nbStates = 3, dist = dist, Par0 = Par0_m2,
             estAngleMean = list(angle=TRUE), stateNames = stateNames)
m2


### view/check model output ###
## For a given model, the function viterbi() computes the most likely state sequence
states <- viterbi(m2)
# derive percentage of time spent in each state
print(table(states)/nrow(cd))

print(ggplot(cd, aes(x = x, y = y, color = factor(states), group = 1)) +
        geom_path() +
        scale_color_manual(values = c("#fec44f", "#4eb3d3", "#238b45"))) +
        theme_bw()


print(AIC(m1, m2)) # the two state model is clearly the better fit in this simple example


#### save the dispersal movement parameters from our top model ####
step <- m1$mle$step
step
step.mean <- step[rownames(step)=="mean", colnames(step)=="dispersing"]
step.sd <- step[rownames(step)=="sd", colnames(step)=="dispersing"]
angle <- m1$mle$angle
angle
angle.mean <- angle[rownames(angle)=="mean", colnames(angle)=="dispersing"]
angle.conc <- angle[rownames(angle)=="concentration", colnames(angle)=="dispersing"]





#### TRAVERSABILITY ####
### now we can simulate dispersal trajectories and calculate some metrics of "traversability"
## in this case, we'll calculate:
# 1. the proportion of each simulated path that falls in agricultural land
# 2. the proportion of simulated paths that cross a major road

## we'll simulate dispersal paths from our test individual's natal or pre-dispersal range (which we simulated earlier)
pre.disp <- disp1
ggplot(pre.disp, aes(x = x, y = y)) + geom_path() + coord_fixed() + theme_bw()


#### simulate landscape data ####

## let's simulate a highway as a SpatialLines object
x <- c(-15,4,7,10,15,20,35)
y <- c( 0, 0,4,10,12,20,30)
hw <- SpatialLines(list(Lines(Line(cbind(x,y)), ID="hw")))

ggplot() + geom_path(data = pre.disp, aes(x = x, y = y), colour = "blue") + 
  geom_path(data = data.frame(x = x, y = y), aes(x = x, y = y)) +
  coord_fixed() + theme_bw()


## let's also simulate a habitat that is a mix of agricultural and forested land
land <- raster(nrows = 50, ncols = 50, xmn = -15, xmx = 35,
               ymn = -15, ymx = 35, resolution = 1)
land[] <- 2
land[c(30:45), c(25:45)] <- 1
land[c(15:25), c(10:20)] <- 1
land[c(40:50), c(0:30)] <- 1
land[c(7:12), c(0:40)] <- 1

land.spdf <- as(land, "SpatialPixelsDataFrame")
land_df <- as.data.frame(land.spdf)
colnames(land_df) <- c("value", "x", "y")
land_df$class <- ifelse(land_df$value==1, "ag", "forest")

## plot the trajectory, road, and landscape together
ggplot() + geom_raster(data = land_df, aes(x = x, y = y, fill = class)) +
  scale_fill_manual(values = c("#f7fcb9", "#31a354")) +
  theme_bw() + 
  geom_path(data = pre.disp, aes(x = x, y = y), colour = "blue") + 
  geom_path(data = data.frame(x = x, y = y), aes(x = x, y = y)) +
  coord_fixed()


#### simulate dispersal paths ####
### now we can simulate our dispersal paths
## take a sample of starting locations from empirical trajectory
# we'll just do 15 here, but I recommend doing many more with a real dataset
starts <- pre.disp[sample(1:nrow(pre.disp), size = 15, replace = T), c("x", "y")]

# make a list with the randomly selected starting points  
starts.list <- split(starts, seq(nrow(starts)))
starts.list <- lapply(starts.list, function(y) c(y$x, y$y))
  
## simulate dispersal trajectories based on HMM movement parameters
# we just have the one set of dispersal movement parameters for this example, 
# but with a full dataset, we'd take the average movement parameters across individuals
sim.par <- list(step = c(step.mean, step.sd),
                  angle = c(angle.mean, angle.conc))
dist = list(step = "gamma", angle = "vm")

# now we'll make the actual simulated movements
# we'll simulate 10 timesteps per path, since that was the length of our original simulated dispersal event
set.seed(531)
temp.sims <- simData(nbAnimals = 15, nbStates = 1, dist = dist,
                       Par = sim.par, obsPerAnimal = 10, initialPosition = starts.list)
  
# view simulated movements
ggplot() + geom_raster(data = land_df, aes(x = x, y = y, fill = class)) +
  scale_fill_manual(values = c("#f7fcb9", "#31a354")) +
  geom_path(data = pre.disp, aes(x = x, y = y), colour = "blue") + 
  geom_path(data = data.frame(x = x, y = y), aes(x = x, y = y)) +
  geom_path(data = temp.sims, aes(x = x, y = y, group = ID, colour = ID), alpha = 0.8) +
  theme(legend.position = "none") +
  theme_bw() +  coord_fixed()

  
# create dataframe that also links each simulate point to land use (1 = ag, 2 = forest)
temp.simraster <- prepData(temp.sims[,c("ID", "x", "y")], type = "UTM", coordNames = c("x", "y"), spatialCovs = list(land = land))
  
#### traversability metrics ####
### now we calculate our traversability metrics with our simulated dispersal paths
## calculate proportion of each simulated path that was in ag
# we can make this a function for flexibility (e.g., if more than one land use class in a raster is classified as "ag," etc.)
prop.crop.point <- function(x){
  crop.points <- 1
  if(any(x %in% crop.points)){
    out <- length(which(x %in% crop.points))/length(x)
  }else{
    out <- 0
  }
  
  return(out)
}
prop.crop.df <- ddply(temp.simraster, .(ID), function(y) prop.crop.point(y$land))
prop.crop.df  


## calculate proportion of simulated dispersal paths that intersects our simulated highway
# make simulated trajectories into one SpatialLines object
lines.list <- vector(mode = "list", length = length(unique(temp.simraster$ID)))
for(j in 1:length(unique(temp.simraster$ID))){
  temp.path <- as.data.frame(temp.simraster[temp.simraster$ID==j,])
  temp.line <- Line(temp.path[,c("x", "y")])
  lines.list[[j]] <- Lines(list(temp.line), ID = as.character(j))
  names(lines.list)[j] <- as.character(j)
  
}
sp.lines <- SpatialLines(lines.list)
plot(sp.lines)
  
## check which paths intersect our simulated highway
intersects.roads <- gIntersects(sp.lines, hw, byid = T, returnDense = T)
intersects.roads <- colSums(intersects.roads)
num.intersects.roads <- length(intersects.roads[intersects.roads>0])
num.intersects.roads/15
  
# we can add this information to our simulated dataset and plot it
temp.simraster$int.road <- ifelse(temp.simraster$ID%in%which(intersects.roads>0), "intersects", "no intersection")

ggplot() + geom_raster(data = land_df, aes(x = x, y = y, fill = class)) +
  scale_fill_manual(values = c("#f7fcb9", "#31a354")) +
  geom_path(data = pre.disp, aes(x = x, y = y), colour = "blue") + 
  geom_path(data = data.frame(x = x, y = y), aes(x = x, y = y)) +
  geom_path(data = temp.simraster, aes(x = x, y = y, group = ID, colour = int.road), alpha = 0.8) +
  scale_color_manual(values = c("red", "blue")) + 
  theme(legend.position = "none") +
  theme_bw() +  coord_fixed()


### we can store back information about these traversability metrics and use them in statistical modeling of dispersal probabilities and distances