## verify_class.R
# 
#========================================================	
# ---
### title: Verify classification
# author: Marie Gilbertson
# date: ""
#---
###  Preamble	
# 
# What this code does:
# 1. Function that loops through all HR_dispersal figures to check and update dispersal classifications


verify_class <- function(fd = full.out, # dispersal detection output data frame
                         dir = "HR_dispersal_test" # directory where dispersal detection figures are saved
                         ){

  fd$new.class <- NA
  
  for(i in 1:nrow(fd)){
    print(i)
    figdat_file <- paste("Figures/", dir, "/Figure_data/figuredata_id_", fd$id[i], ".Rdata", sep = "")
    figdat <- get(load(figdat_file))
    print(figdat)
    
    
    
    new.class <- fd$class[i]
    # check classification and update if needed
    repeat{
      raw <- readline("Class choice good? Choices are: ok, not  ")
    
      if(raw=="not"){
          new <- readline("New class? Choices are: re, di, ex, ts, uk  ")
    
          if(new=="re"){
            new.class <- "resident"
          }else if(new=="di"){
            new.class <- "dispersal"
          }else if(new=="ex"){
            new.class <- "excursion"
          }else if(new=="ts"){
            new.class <- "temp.range.shift"
          }else if(new=="uk"){
            new.class <- "unknown"
          }else{
            new.class <- "ERROR"
          }
    
        }
      if(raw=="ok")break
    }
    
    fd$new.class[i] <- paste(new.class)
    figdat <- annotate_figure(figdat, top = paste0(fd$class[i], " -> ", fd$new.class[i]))
    
    save(figdat, file = paste("Figures/", dir, "/Reclass/Reclass_figure_data/reclassfiguredata_id_", fd$id[i], ".Rdata", sep = ""))
    ggsave(paste("Figures/", dir, "/Reclass/reclass_id_", fd$id[i], ".jpg", sep = ""), figdat, width = 16, height = 8, units = "in", dpi = 300)
  }
  return(fd)
}