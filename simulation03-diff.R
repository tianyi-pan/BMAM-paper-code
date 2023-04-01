#!/usr/bin/env Rscript

## Appendix: Simulation Heterogeneity Among Cluster
## Simulation 3: difference between marginal and conditional models

library(parallel)
library(doParallel)
library(doRNG)

num_for_cores <- parallel::detectCores() - 1

cl <- makeCluster(num_for_cores)
registerDoParallel(cl)

for(K in c(100)){
for (Nk in c(10)) {
for (modeltype in c("conditional","marginal")){
  R <- 500
  clusterExport(cl, c("K","Nk"))
  foreach(i = 1:R, .errorhandling = "pass") %dorng%{
    system(paste("./do_simu_diff.R",as.numeric(K), as.numeric(Nk), as.numeric(i), as.character(modeltype), sep = " "))
  }
}  
}
}
stopCluster(cl)
