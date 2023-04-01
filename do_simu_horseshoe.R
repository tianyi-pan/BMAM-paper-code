#!/usr/bin/env Rscript
### load packages #####################
library(brms)
library(bmam)
library(dplyr)
library(simstudy) # pacakge for correlated uniform distribution

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

## point constrain! at x = 0. 
## smooth function setting 
f1temp <- function(x){(1.5*sin(pi*x)-2*x)}
f1 <- function(x) f1temp(x)-f1temp(0)


f2temp <- function(x){(1.5*cos(2*pi*(x+1/4)))}
f2 <- function(x) f2temp(x)-f2temp(0)



pc_x1 <- 0
pc_x2 <- 0


do_simulation_multiple <- function(K,betachoose, indept, iter, verbose=FALSE) {
  set.seed(iter)
  Nk <- 10
  sigma0 = 2;sigma3 = 1;rho=0.5
  if(betachoose == "sparsity1") {
    beta <- c(0,1,1,
              2,0,0,0,0,0,0,0)
    betal <- c(2,0,0,0,0,0,0,0)
  }else{
    beta <- c(0,1,1,
              1,1,1,1,0,0,0,0)
    betal <- c(1,1,1,1,0,0,0,0)
  }
  
  
  ## form the true covariance matrix
  trueSigma <- matrix(c(sigma0^2,
                        rho*sigma0*sigma3,rho*sigma0*sigma3,
                        sigma3^2),nrow=2,ncol=2)
  if(indept){
    ### 1. x1,x2 are independent with other predictors

    x1 <- round(runif(K*Nk,-1,1),3)
    x1 <- rep(sample(x1, K), each = Nk)
    
    
    df2 = genCorGen(K, nvars = Nk, params1 = rep(-1,Nk), 
                    params2 = rep(1,Nk), dist = "uniform", rho = 0.5, corstr = "cs", wide = TRUE)
    df2 = as.data.frame(df2[,-1])
    x2 <- c(t(df2))
    
    corr = .5
    p <- 8
    Covmatrix <- outer(
      1:p, 1:p,
      function(x, y) {
        corr^abs(x - y)
      }
    )
    X <- MASS::mvrnorm(K*Nk, rep(0, p), Covmatrix)
    colnames(X) <- c("x3","x4","x5","x6","x7","x8","x9","x10")
    
    # corr <- 0
    # X <- sapply(1:8, function(r.) round(runif(K*Nk,-1,1),3))
    # colnames(X) <- c("x3","x4","x5","x6","x7","x8","x9","x10")
    
  }else{
    corr <- 0.5
    ### 2. x1,x2 are correlated with other predictors
    X = genCorGen(K*Nk, nvars = 10, 
                   params1 = rep(-1,10), params2 = rep(1,10), 
                   dist = "uniform", rho = corr, corstr = "cs", wide = TRUE)
    X <- as.data.frame(X[,-1])
    x1 <- X[,1]
    x1 <- rep(sample(x1, K), each = Nk)
    x2 <- X[,2]
    X <- X[,3:10]
    colnames(X) <- c("x3","x4","x5","x6","x7","x8","x9","x10")
  }
  
  
  fx1 <- f1(x1)
  fx2 <- f2(x2)
  

  ## generate data based on marginal model  
  mean.formula <- ~fx1+fx2+x3+x4+x5+x6+x7+x8+x9+x10
  # mean.formula <- ~fx1+fx2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17
  # +x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+x21+x22
  lv.formula <- ~1+x3
  
  nclust <- rep(Nk,K) # number of units in each cluster
  id <- rep(seq(K), nclust) # id of each unit
  data  <- data.frame(id, x1, fx1,x2,fx2)
  data <- cbind(data, X)
  data  <- data[order(data$id),] # order data by id and time
  
  
  data_gen <- bmam::GenBinaryY(mean.formula = mean.formula, lv.formula = lv.formula, 
                               beta=beta, Sigma=trueSigma, id=id, data=data, q=120, 
                               Yname = "y")
  
  ## prepare data
  dat2pred <- data.frame(
    x1 = c(seq(-1,1,length = 100),rep(0,100)),
    x2 = c(rep(0,100),seq(-1,1,length = 100)),
    x3 = 0,
    x4 = 0, 
    x5 = 0,
    x6 = 0, 
    x7 = 0,
    x8 = 0, 
    x9 = 0,
    x10 = 0, 
    varname = c(rep("x1",100), rep("x2", 100))
  )
  
  # pc_x1 <- 0
  # pc_x2 <- 0
  f1_pc_c <- f1(pc_x1)
  f2_pc_c <- f2(pc_x2)
  
  fx1_true <- f1(dat2pred$x1) - f1_pc_c
  fx2_true <- f2(dat2pred$x2) - f2_pc_c
  
  ### fit model
  
  horseshoe_term <- brms.horseshoe(shrinkage.term = y~x3+x4+x5+x6+x7+x8+x9+x10, 
                                   nonshrinkage.term = ~s(x1, pc = pc_x1) + s(x2, pc = pc_x1)+ (1+x3|id))
  
  # 
  # horseshoe_term <- brms.horseshoe(shrinkage.term = y~x3+x4+x5+x6+x7+x8+x9+x10, 
  #                                  nonshrinkage.term = ~s(x1) + s(x2)+ (1+x3+x4|id))
  
  
  
  model_brms_smooth <- brm(horseshoe_term$formula,
                           data = data_gen, family = bernoulli(link = "logit"), cores = 1, seed = 4321,
                           warmup = 1000, iter = 2000, chains = 4, refresh=0, backend = "cmdstanr",
                           save_pars = save_pars(all = TRUE),
                           prior=c(horseshoe_term$prior))
  # bmam.fit <- bmam(model_brms_smooth)
  bmam.fit <- bmam(object = model_brms_smooth, centered = F, posterior = FALSE,
                   k=100, CIType="ETI", CI = 0.95, preddat = dat2pred, est.fun = TRUE)
  
  
  
  model_brms_flat <- brm(bf(y ~ x3+x4+x5+x6+x7+x8+x9+x10+ s(x1, pc = pc_x1)+s(x2, pc = pc_x1)+ (1+x3|id)),
                         data = data_gen, family = bernoulli(link = "logit"), cores = 1, seed = 4321,
                         save_pars = save_pars(all = TRUE),
                         warmup = 1000, iter = 2000, chains = 4, refresh=0, backend = "cmdstanr")
  # bmam.fit_flat <- bmam(model_brms_flat)
  bmam.fit_flat <- bmam(object = model_brms_flat, centered = F, posterior = FALSE,
                        k=100, CIType="ETI", CI = 0.95, preddat = dat2pred, est.fun = TRUE)
  
  
  
  ### save results
  
  # dat2pred_flat <- bmam.fit_flat$Preddat
  fitteddat <- tidyr::tibble(
    ## Simulation parameters
    iter = iter, 
    corr = corr,
    var = rep(c("x1","x2"),each = 100),
    ## 1. BMAM
    bmamfitted = bmam.fit$Predicted_Summary$M,
    bmamlower = bmam.fit$Predicted_Summary$LL,
    bmamupper = bmam.fit$Predicted_Summary$UL,
    bmamcovrx1 = fx1_true >= bmamlower &  fx1_true <= bmamupper,
    bmamcovrx2 =  fx2_true >= bmamlower &  fx2_true <= bmamupper,
    bmambiasx1 = bmamfitted -  fx1_true,
    bmambiasx2 = bmamfitted -  fx2_true,
    
    ## 1. BMAM flat
    bmamfitted_flat = bmam.fit_flat$Predicted_Summary$M,
    bmamlower_flat = bmam.fit_flat$Predicted_Summary$LL,
    bmamupper_flat = bmam.fit_flat$Predicted_Summary$UL,
    bmamcovrx1_flat = fx1_true >= bmamlower_flat & fx1_true <= bmamupper_flat,
    bmamcovrx2_flat =  fx2_true >= bmamlower_flat &  fx2_true <= bmamupper_flat,
    bmambiasx1_flat = bmamfitted_flat -  fx1_true,
    bmambiasx2_flat = bmamfitted_flat -  fx2_true,
  )
  
  label_horseshoe <- c("a_x3","a_x4","a_x5","a_x6","a_x7","a_x8","a_x9","a_x10")
  label_flat <- c(c("x3","x4","x5","x6","x7","x8","x9","x10"))
  cond_horseshoe <- lapply(label_horseshoe, function(name){
    out <- brmsmargins::bsummary(fixef(model_brms_smooth,summary = FALSE)[,name], CIType="ETI", CI = 0.95)
    out[, Label := name]
    out
  }) 
  cond_horseshoe <- do.call(rbind, cond_horseshoe)  
  cond_flat <- lapply(label_flat, function(name){
    out <- brmsmargins::bsummary(fixef(model_brms_flat,summary = FALSE)[,name], CIType="ETI", CI = 0.95)
    out[, Label := name]
    out
  }) 
  cond_flat <- do.call(rbind, cond_flat)  
  
  linear.estimate <- tidyr::tibble(
    iter = iter, 
    indept = indept,
    corr = corr,
    label = bmam.fit_flat$Summary$Label[2:9],
    true = c(betal),
    # horseshoe = bmam.fit$Summary$M[1:8],
    horseshoe = filter(bmam.fit$Summary, Label %in% label)$M,
    # horseshoebias = bmam.fit$Summary$M[1:8] - betal,
    horseshoebias = horseshoe - betal,
    # flat = bmam.fit_flat$Summary$M[2:9],
    flat = filter(bmam.fit_flat$Summary, Label %in% label)$M,
    flatbias = flat - betal,
    
    ## conditional
    horseshoe_cond = cond_horseshoe$M,
    horseshoebias_cond = horseshoe_cond - betal,
    flat_cond = cond_flat$M,
    flatbias_cond = flat_cond - betal,
  )
  
  filepath <- paste0("results/","C2_horseshoe")
  if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)

  
  fileclass <- paste0(filepath,"/U_","K",as.character(K),"_corr",as.character(corr),"_p10_",betachoose, ifelse(indept, "_indept","_correlated"))
  save_results(paste0(fileclass,"_smooth.csv"),fitteddat)
  save_results(paste0(fileclass,"_linear.csv"),linear.estimate)
  
}


args <- commandArgs(trailingOnly = TRUE)
tryCatch(
  expr = {
    do_simulation_multiple(K = as.numeric(args[1]), 
                           betachoose = as.character(args[2]), 
                           indept = as.numeric(args[3]),
                           iter = as.numeric(args[4]))
  },
  error = function(e){
    message("error")
  }
)  