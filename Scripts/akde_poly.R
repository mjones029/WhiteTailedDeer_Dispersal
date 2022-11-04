## akde_poly.R
# 
#========================================================	
# ---
### title: AKDE polygon
# author: Marie Gilbertson
# date: "12/01/2021"
#---
###  Preamble	
# 
# What this code does:
# 1. Function that uses ctmm to estimate AKDE and outputs area and polygon

akde_poly <- function(pre.data){

  pre.m <- move(x=pre.data$x, y=pre.data$y, 
                time=pre.data$date, 
                proj=CRS("+init=epsg:3071"), 
                data=pre.data, animal=pre.data$id)
  # plot(pre.m)
  
  
  pre.t <- suppressMessages(as.telemetry(pre.m))
  
 
  # -------------------------------------------------------------------------
  ## Step 3. Selecting the best-fit movement model through model selection
  
  # Calculate an automated model guesstimate:
  GUESS1 <- ctmm.guess(pre.t, interactive = FALSE)
  
  # Automated model selection, starting from GUESS:
  FIT1_ML <- ctmm.select(pre.t, GUESS1, method = 'ML', cores = 4)
  
  
  # -------------------------------------------------------------------------
  ## Step 4. -- Feeding a movement model into the home range estimator
  
  # estimate aKDE (can modify this to best suit your data)
  UD1_ML <- akde(pre.t, FIT1_ML)
  
  area <- summary(UD1_ML, units = F)$CI[2] # home range area estimation
  
  
  pl <- SpatialPolygonsDataFrame.UD(UD1_ML, level.UD = 0.95, level = 0.95)
  pl.est <- pl[2,]
  # plot(pl.est)
  
  out <- list(area, pl.est)
  names(out) <- c("area", "pl.est")
  return(out)
}
