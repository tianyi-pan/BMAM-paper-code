#!/usr/bin/env Rscript
### load packages #####################
# devtools::install_github('tianyi-pan/mam')
library(mam)
library(mgcv)
library(gamm4)
library(tidyverse)
library(brms)
library(brmsmargins)
library(data.table)
library(bmam)
library(splines)
# library(tmvtnorm) # package for truncated norm
library(simstudy) # pacakge for correlated uniform distribution


## point constrain. at x = 0. 

### function for data generation ###################

## smooth function setting 
f1temp <- function(x){(1.5*sin(pi*x)-2*x)}
f1 <- function(x) f1temp(x)-f1temp(0)


f2temp <- function(x){(1.5*cos(2*pi*(x+1/4)))}
f2 <- function(x) f2temp(x)-f2temp(0)

pc_x1 <- 0
pc_x2 <- 0

## function
get_data_marginal <- function(SSmat,x1,x2,x3,K,Nk, beta3){
  fx1 <- f1(x1)
  fx2 <- f2(x2)
  beta <- c(0,1,1,beta3)
  
  
  
  mean.formula <- ~fx1+fx2+x3
  lv.formula <- ~1+x3
  
  nclust <- rep(Nk,K) # number of units in each cluster
  id <- rep(seq(K), nclust) # id of each unit
  data  <- data.frame(id, x1, fx1, x2, fx2, x3) # generate data (without y)
  data  <- data[order(data$id),] # order data by id and time
  
  
  data_gen <- bmam::GenBinaryY(mean.formula = mean.formula, lv.formula = lv.formula, 
                               beta=beta, Sigma=SSmat, id=id, data=data, q=120, 
                               Yname = "y")
  names(data_gen)[c(8,9)] <- c("intercepts", "slopes")
  
  return(data_gen)
}


get_data_conditional <- function(SSmat,x1,x2,x3,K,Nk, beta3){
  fx1 <- f1(x1)
  fx2 <- f2(x2)
  beta <- c(0,1,1,beta3)
  
  
  
  mean.formula <- ~fx1+fx2+x3
  lv.formula <- ~1+x3
  
  nclust <- rep(Nk,K) # number of units in each cluster
  id <- rep(seq(K), nclust) # id of each unit
  
  X <- cbind(1,fx1,fx2,x3)
  Z <- cbind(1,x3)
  
  ranef <- mvtnorm::rmvnorm(n = K, mean = rep(0,dim(SSmat)[1]), sigma = SSmat)
    
  
  ranef <- do.call("rbind",
                   apply(ranef, 1, 
                         function(row.) matrix(rep(row., times = Nk), nrow = Nk, byrow = TRUE), simplify = F))
  
  

  expit <- function(aa) {
    exp_aa <- exp(aa)
    exp_aa/(1+exp_aa)
  }
  
  MuC <- expit(X %*% beta + apply(Z * ranef, 1, sum))
  Y   = as.integer(runif(K*Nk) < MuC)
  
  data_gen <- data.frame("id" = id, 
                         "x1" = x1, "fx1" = fx1, 
                         "x2" = x2, "fx2" = fx2,
                         "x3" = x3,  "y" = Y,
                         "intercepts" = ranef[,1], "slopes" = ranef[,2])
  return(data_gen)
}


### functions for saving results ##################
save_results <- function(filename, mc.summary){
  # save results of simulation
  for( i in seq_len(nrow(mc.summary)) ){
    sline <- mc.summary[i,]
    if(file.exists(filename)){
      write.table(sline, file = filename, sep = ",", append = TRUE, quote = FALSE,
                  col.names = FALSE, row.names = FALSE)
    }else{
      write.csv(sline, filename, row.names = FALSE)
    }
  }
}



### simulation function ###########################

do_simulation_multiple <- function(K,Nk,beta3,sigma0,sigma3,rho,iter,verbose=FALSE, modeltype) {
  
  set.seed(iter)
  
  ## choose where the data are generated
  if(modeltype == "conditional"){
    get_data <- get_data_conditional # from conditional model
  } else if(modeltype == "marginal"){
    get_data <- get_data_marginal # from marginal model
  } else{
    stop("provide correct model type.")
  }
  
  ## form the true covariance matrix
  trueSigma <- matrix(c(sigma0^2,
                        rho*sigma0*sigma3,rho*sigma0*sigma3,
                        sigma3^2),nrow=2,ncol=2)
  ## generate data
  
  x3 <- round(runif(K*Nk,-1,1),3)
  
  #### cluster covariates
  # x3 <- rep(sample(x3, K), each = Nk)
  x1 <- round(runif(K*Nk,-1,1),3)
  x1 <- rep(sample(x1, K), each = Nk)
  
  # for uniform distribution params1 is minimum, params2 is maximum
  # https://stats.stackexchange.com/questions/66610/generate-pairs-of-random-numbers-uniformly-distributed-and-correlated
  
  df2 = genCorGen(K, nvars = Nk, params1 = rep(-1,Nk), 
                  params2 = rep(1,Nk), dist = "uniform", rho = 0.5, corstr = "cs", wide = TRUE)
  df2 = as.data.frame(df2[,-1])
  x2 <- c(t(df2))
  
  dat <- get_data(trueSigma,x1,x2,x3,K,Nk, beta3)
  
  
  

  f1_pc_c <- f1(pc_x1)
  f2_pc_c <- f2(pc_x2)
  
  ## prepare data
  dat2pred <- data.frame(
    x1 = c(seq(-1,1,length = 100),rep(0,100)),
    x2 = c(rep(0,100),seq(-1,1,length = 100)),
    x3 = 0,
    varname = c(rep("x1",100), rep("x2", 100))
  )
  
  dat_x0 <- data.frame(
    x1 = 0,
    x2 = 0,
    x3 = 0,
    varname = c(rep("x1",100), rep("x2", 100))
  )
  fx1 <- f1(dat2pred$x1) - f1_pc_c
  fx2 <- f2(dat2pred$x2) - f2_pc_c
  
  ### fit models
  # 1. BMAM
  model_brms <- brms::brm(brms::bf(y ~  x3 + s(x1, pc = pc_x1) + s(x2, pc = pc_x2) + (1+x3|id)),
                          data = dat, family = "bernoulli", cores = 1, seed = 4321,
                          save_pars = save_pars(all = TRUE),
                          warmup = 1000, iter = 2000, chains = 4, refresh=0, backend = "cmdstanr")
  
  bmamfit <- bmam::bmam(object = model_brms, centered = F, posterior = FALSE,
                        k=100, CIType="ETI", CI = 0.95, preddat = dat2pred, est.fun = TRUE)
  

  
  # 2. Bayesian GAMM (conditional model)
  bgammfit <-  
    do.call(rbind, apply(
    t(fitted(
    object = model_brms, newdata = dat2pred,
    re_formula = NA, scale = "linear",
    summary=FALSE) - 
    fitted(
      object = model_brms, newdata = dat_x0,
      re_formula = NA, scale = "linear",
      summary=FALSE)),
    1, brmsmargins::bsummary, CIType="ETI", CI = 0.95))
  
    

  
  
  ### get fitted values
  fitteddat <- tidyr::tibble(
    ## Simulation parameters
    K=K,Nk=Nk,beta3=beta3, sigma0=sigma0,sigma3=sigma3,rho=rho,iter=iter,
    ## true values
    fx1 = fx1, fx2 = fx2,
    var = rep(c("x1","x2"),each = 100),
    ## 1. BMAM
    bmamfitted = bmamfit$Predicted_Summary$M,
    bmamlower = bmamfit$Predicted_Summary$LL,
    bmamupper = bmamfit$Predicted_Summary$UL,
    bmamcovrx1 = fx1 >= bmamlower & fx1 <= bmamupper,
    bmamcovrx2 = fx2 >= bmamlower & fx2 <= bmamupper,
    bmambiasx1 = bmamfitted - fx1,
    bmambiasx2 = bmamfitted - fx2,
    
    ## 2. Bayesian GAMM (conditional model)
    bgammfitted = bgammfit$M,
    bgammlower = bgammfit$LL,
    bgammupper = bgammfit$UL,
    bgammcovrx1 = fx1 >= bgammlower & fx1 <= bgammupper,
    bgammcovrx2 = fx2 >= bgammlower & fx2 <= bgammupper,
    bgammbiasx1 = bgammfitted - fx1,
    bgammbiasx2 = bgammfitted - fx2
  )
  
  
  
  
  filepath <- paste0("results/","C3_diff/",modeltype)
  if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)
  fileclass <- file.path(filepath,paste0("K",K,"_Nk",Nk,"_beta3",beta3,"_sigma0",sigma0,"_sigma3",sigma3,"_rho",rho))
  
  save_results(paste0(fileclass,"_smooth.csv"),fitteddat)
  

}

args <- commandArgs(trailingOnly = TRUE)
tryCatch(
  expr = {
    do_simulation_multiple(K = as.numeric(args[1]), Nk = as.numeric(args[2]), 
                           beta3 = 1,sigma0 = 1,sigma3 = 0.5,rho=0.5,iter = as.numeric(args[3]),
                           modeltype = as.character(args[4]))
    do_simulation_multiple(K = as.numeric(args[1]), Nk = as.numeric(args[2]), 
                           beta3 = 1,sigma0 = 2,sigma3 = 1,rho=0.5,iter = as.numeric(args[3]),
                           modeltype = as.character(args[4]))
  },
  error = function(e){
    message("error")
  }
)  

