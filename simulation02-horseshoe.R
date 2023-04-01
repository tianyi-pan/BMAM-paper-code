#!/usr/bin/env Rscript

## Appendix: Extension marginal variable selection with the horseshoe prior
## Simulation 2: horseshoe simulation

library(parallel)
library(doParallel)
library(doRNG)

num_for_cores <- parallel::detectCores()-1

cl <- makeCluster(num_for_cores)
registerDoParallel(cl)

for(K in c(100)){
  for (betachoose in c("sparsity1","sparsity4")) {
    for (indept in c(1)){
      R <- 500
      foreach(i = 1:R, .errorhandling = "pass") %dorng%{
        system(paste("./do_simu_horseshoe.R",as.numeric(K), as.character(betachoose), 
                     as.numeric(indept), as.numeric(i), sep = " "))
      }
    }  
  }
}
stopCluster(cl)

