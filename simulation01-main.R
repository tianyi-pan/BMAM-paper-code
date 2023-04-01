#!/usr/bin/env Rscript

## Section 3.4 & Section 5.3
## Simulation 1: main simulation. 

library(parallel)
library(doParallel)
library(doRNG)

num_for_cores <- parallel::detectCores()-1
# num_for_cores <- 85

cl <- makeCluster(num_for_cores) 
registerDoParallel(cl)



for(K in c(300,500)){
for (Nk in c(10,20)) {
for (modeltype in c("conditional","marginal")){
  R <- 500
  clusterExport(cl, c("K","Nk"))
  foreach(i = 1:R, .errorhandling = "pass") %dorng%{
    system(paste("./do_simu_main.R",as.numeric(K), as.numeric(Nk), as.numeric(i), as.character(modeltype), sep = " "))
  }
}  
}
}
stopCluster(cl)