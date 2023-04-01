#!/usr/bin/env Rscript
### load packages #####################
# **Note: install the mam package:
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

## point constrain! at x = 0. 

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
  
    

  # 3. Bayesian GAM
  bgamfit <- brms::brm(brms::bf(y ~  x3 + s(x1, pc = pc_x1) + s(x2, pc = pc_x2)),
                       data = dat, family = "bernoulli", cores = 1, seed = 4321,
                       warmup = 1000, iter = 2000, chains = 4, refresh=0, backend = "cmdstanr")
  bgamfitted_sample <- fitted(object = bgamfit, newdata = dat2pred, scale = "linear", summary = FALSE)
  bgamfitted_sample_x0 <- fitted(object = bgamfit, newdata = dat_x0, scale = "linear", summary = FALSE)
  bgamfitted_summary <- data.table::as.data.table(do.call(rbind, apply(t(bgamfitted_sample-bgamfitted_sample_x0), 
                                                                       1, brmsmargins::bsummary, 
                                                                       CIType="ETI", CI = 0.95)))
  
  
  # 4. MAM 
  # get estimates
  tryCatch(
    expr = {
      themam <- NULL
      try( ## error handling due to optimizer
        themam <- mam::mam(smooth = list(s(x1, pc = pc_x1), s(x2, pc = pc_x2)),
                           re = y ~ (1+x3|id),
                           fe = ~ x3,
                           dat = dat,
                           margdat = dat,
                           preddat = dat2pred,
                           est.fun = TRUE, # change mam code
                           control = mam_control(
                             method = 'trust',
                             varmethod = 1,
                             verbose = FALSE,
                             retcond = TRUE))
      )
      if(is.null(themam)){
        themam <- mam::mam(smooth = list(s(x1, pc = pc_x1), s(x2, pc = pc_x2)),
                           re = y ~ (1+x3|id),
                           fe = ~ x3,
                           dat = dat,
                           margdat = dat,
                           preddat = dat2pred ,
                           est.fun = TRUE, # change mam code
                           control = mam_control(
                             method = 'BFGS',
                             varmethod = 1,
                             verbose = FALSE,
                             retcond = TRUE))
      }
    },
    error = function(e){
      themam <- NULL
    }
  )  

  
  
  # 5. GAMM (conditional)
  gamm <- gamm4::gamm4(y ~ s(x1, pc = pc_x1) + s(x2, pc = pc_x2) + x3,random=~(1+x3|id), data = dat,family = binomial())
  # gammfit <- predict(gamm$gam, se.fit=TRUE, type="terms")
  gammfit <- predict(gamm$gam, newdata=dat2pred, se.fit=TRUE, type="terms")
  # gammfit_x0 <- predict(gamm$gam,newdata=dat_x0,se.fit=TRUE,type="terms")
  gammfitted <- apply(gammfit$fit,1,sum)
  gammfitted_se <- apply(gammfit$se.fit,1,sum)
  gammpred <- ranef(gamm$mer)$id
  
  
  # 6. GAM
  gam1 <- mgcv::gam(y ~ s(x1, pc = pc_x1) + s(x2, pc = pc_x2) + x3,data=dat,family=binomial(),method="REML")
  gamfit <- predict(gam1,newdata=dat2pred,se.fit=TRUE,type="terms")
  gamfitted <- apply(gamfit$fit,1,sum)
  gamfitted_se <- apply(gamfit$se.fit,1,sum)
  
  
  
  
  ### get fitted values
  fitteddat <- tidyr::tibble(
    ## Simulation parameters
    K=K,Nk=Nk,beta3=beta3, sigma0=sigma0,sigma3=sigma3,rho=rho,iter=iter,
    ## true values
    fx1 = fx1, fx2 = fx2,
    ## mam
    ifmam = as.numeric(!is.null(themam)),
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
    bgammbiasx2 = bgammfitted - fx2,
    
    ## 3. Bayesian GAM 
    bgamfitted = bgamfitted_summary$M,
    bgamlower =  bgamfitted_summary$LL,
    bgamupper =  bgamfitted_summary$UL,
    bgamcovrx1 = fx1 >= bgamlower & fx1 <= bgamupper,
    bgamcovrx2 = fx2 >= bgamlower & fx2 <= bgamupper,
    bgambiasx1 =  bgamfitted_summary$M - fx1,
    bgambiasx2 =  bgamfitted_summary$M - fx2,
    
    ## 4-1. mam marginal
    mammargfitted = if(is.null(themam)) rep(0,200) else as.numeric(themam$mam$fitted),
    mammargse = if(is.null(themam)) rep(0,200) else as.numeric(themam$mam$fitted_se),
    mammarglower = mammargfitted - 1.96*mammargse,
    mammargupper = mammargfitted + 1.96*mammargse,
    mammargcovrx1 = fx1 >= mammarglower & fx1 <= mammargupper,
    mammargcovrx2 = fx2 >= mammarglower & fx2 <= mammargupper,
    mammargbiasx1 = mammargfitted - fx1,
    mammargbiasx2 = mammargfitted - fx2,
    
    ## 4-2. mam conditional
    mamcondfitted = if(is.null(themam)) rep(0,200) else as.numeric(themam$conditional$fitted),
    mamcondse = if(is.null(themam)) rep(0,200) else as.numeric(themam$conditional$fitted_se),
    mamcondlower = mamcondfitted - 1.96*mamcondse,
    mamcondupper = mamcondfitted + 1.96*mamcondse,
    mamcondcovrx1 = fx1 >= mamcondlower & fx1 <= mamcondupper,
    mamcondcovrx2 = fx2 >= mamcondlower & fx2 <= mamcondupper,
    mamcondbiasx1 = mamcondfitted - fx1,
    mamcondbiasx2 = mamcondfitted - fx2,
    
    ## 5. gamm
    gammcondfitted = gammfitted,
    gammcondse = gammfitted_se,
    gammcondlower = gammcondfitted - 1.96*gammcondse,
    gammcondupper = gammcondfitted + 1.96*gammcondse,
    gammcondcovrx1 = fx1 >= gammcondlower & fx1 <= gammcondupper,
    gammcondcovrx2 = fx2 >= gammcondlower & fx2 <= gammcondupper,
    gammcondbiasx1 = gammcondfitted - fx1,
    gammcondbiasx2 = gammcondfitted - fx2,
    
    ## 6. gam
    gammargfitted = as.numeric(gamfitted),
    gammargse = as.numeric(gamfitted_se),
    gammarglower = gammargfitted - 1.96*gammargse,
    gammargupper = gammargfitted + 1.96*gammargse,
    gammargcovrx1 = fx1 >= gammarglower & fx1 <= gammargupper,
    gammargcovrx2 = fx2 >= gammarglower & fx2 <= gammargupper,
    gammargbiasx1 = gammargfitted - fx1,
    gammargbiasx2 = gammargfitted - fx2,
    
  )
  
  
  
  
  ### random effects
  b0 <- dat$intercepts[seq(1,nrow(dat),by=Nk)] ## cluster level random intercept
  b1 <- dat$slopes[seq(1,nrow(dat),by=Nk)] ## cluster level random slope
  bki <- dat$intercepts+dat$slopes*dat$x3 ## unit level random effect: b0+x3*b1
  
  ## mam
  if(is.null(themam)){
    mamMSEPb0 <- 0
    mamMSEPb1 <- 0
    mamMSEPbki <- 0
  }else{
    mamMSEPb0 <- mean((themam$conditional$predU[,1]-b0)^2)
    mamMSEPb1 <- mean((themam$conditional$predU[,2]-b1)^2)
    mamMSEPbki <- mean( ( apply(themam$conditional$predU[dat$id,]*cbind(1,dat$x3),1,sum)-bki )^2)
  }
  
  ## gamm
  gammMSEPb0 <- mean((gammpred[,1]-b0)^2)
  gammMSEPb1 <- mean((gammpred[,2]-b1)^2)
  gammMSEPbki <- mean( c(( apply(gammpred[dat$id,]*cbind(1,dat$x3),1,sum)-bki )^2))
  
  ## bmam & bgamm
  bgammpred <- brms::ranef(model_brms)$id
  bgammMSEPb0 <- mean((bgammpred[,1,"Intercept"]-b0)^2)
  bgammMSEPb1 <- mean((bgammpred[,1,"x3"]-b1)^2)
  bgammMSEPbki <- mean( c(( apply(bgammpred[dat$id,1,c("Intercept","x3")]*cbind(1,dat$x3),1,sum)-bki )^2))
  
  
  
  
  preddat <- tidyr::tibble(
    K=K,Nk=Nk,beta3=beta3, sigma0=sigma0,sigma3=sigma3,rho=rho,iter=iter,
    
    ifmam = as.numeric(!is.null(themam)),
    
    mamMSEPb0=mamMSEPb0,  
    mamMSEPb1=mamMSEPb1,  
    mamMSEPbki=mamMSEPbki,
    
    gammMSEPb0=gammMSEPb0,
    gammMSEPb1=gammMSEPb1,
    gammMSEPbki=gammMSEPbki,
    
    bgammMSEPb0=gammMSEPb0,
    bgammMSEPb1=gammMSEPb1,
    bgammMSEPbki=gammMSEPbki
  )
  
  
  
  
  ### estimate variance components
  
  ## mam
  if(is.null(themam)){
    Sig <- 0
    mamtheta <- c(0,0,0,0)
  }else{
    Sig <- themam$conditional$theta
    mamtheta <- c(sqrt(Sig[1,1]),Sig[1,2]/sqrt(Sig[1,1]*Sig[2,2]),sqrt(Sig[2,2]))
  }
  
  ## gamm
  gammV <- VarCorr(gamm$mer)$id
  gammtheta <- c(sqrt(gammV[1,1]),gammV[1,2]/sqrt(gammV[1,1]*gammV[2,2]),sqrt(gammV[2,2]))
  ## bgamm
  bgammV <- VarCorr(model_brms)$id
  bgammtheta <- c(bgammV$sd[1,1],bgammV$cor[1,1,"x3"],bgammV$sd[2,1])
  
  names(mamtheta) <- names(gammtheta) <- names(bgammtheta) <- c("sd-Intercept","cor","sd-x3")
  
  
  
  varcomp <- tidyr::tibble(
    K=K,Nk=Nk,beta3=beta3,sigma0=sigma0,sigma3=sigma3,rho=rho,iter=iter,
    
    ifmam = as.numeric(!is.null(themam)),
    
    mam_sd_int=mamtheta['sd-Intercept'],
    mam_sd_slope=mamtheta['sd-x3'],
    mam_cor=mamtheta['cor'],
    gamm_sd_int=gammtheta['sd-Intercept'],
    gamm_sd_slope=gammtheta['sd-x3'],
    gamm_cor=gammtheta['cor'],
    bgamm_sd_int=bgammtheta['sd-Intercept'],
    bgamm_sd_slope=bgammtheta['sd-x3'],
    bgamm_cor=bgammtheta['cor']
  )
  
  
  
  
  ### x3
  bgammx3 <- brmsmargins::bsummary(fixef(model_brms,summary = FALSE)[,"x3"], CIType="ETI", CI = 0.95)
  bgamx3 <- brmsmargins::bsummary(fixef(bgamfit,summary = FALSE)[,"x3"], CIType="ETI", CI = 0.95)
  
  
  x3estimate <- tidyr::tibble(
    K=K,Nk=Nk,beta3=beta3,sigma0=sigma0,sigma3=sigma3,rho=rho,iter=iter,
    
    ifmam = as.numeric(!is.null(themam)),
    
    # 1. BMAM
    bmambiasx3 = bmamfit$Summary[2,'M'] - beta3,
    bmamlowerx3 = bmamfit$Summary[2,'LL'],
    bmamupperx3 = bmamfit$Summary[2,'UL'],
    bmamcovrx3 = beta3 >= bmamlowerx3 & beta3 <= bmamupperx3,
    
    # 2. Bayesian GAMM (conditional model)
    bgammbiasx3 = bgammx3$M - beta3,
    bgammlowerx3 = bgammx3$LL,
    bgammupperx3 = bgammx3$UL,
    bgammcovrx3 = beta3 >= bgammlowerx3 & beta3 <= bgammupperx3,
    
    # 3. Bayesian GAM
    bgambiasx3 = bgamx3$M - beta3,
    bgamlowerx3 = bgamx3$LL,
    bgamupperx3 = bgamx3$UL,
    bgamcovrx3 = beta3 >= bgamlowerx3 & beta3 <= bgamupperx3,
    
    # 4. mam
    mambiasx3 = ifelse(is.null(themam), NA, themam$mam$coefsmooth[1] - beta3),
    mamlowerx3 = ifelse(is.null(themam), NA, themam$mam$coefsmooth[1] - 1.96*themam$mam$coef_se[2]),
    mamupperx3 = ifelse(is.null(themam), NA, themam$mam$coefsmooth[1] + 1.96*themam$mam$coef_se[2]),
    mamcovrx3 = beta3 >= mamlowerx3 & beta3 <= mamupperx3,
    
    # 5. gamm
    gammbiasx3 = fixef(gamm$mer)[2] - beta3,
    gammlowerx3 = fixef(gamm$mer)[2] - 1.96*summary(gamm$mer)$coefficients[2,2],
    gammupperx3 = fixef(gamm$mer)[2] + 1.96*summary(gamm$mer)$coefficients[2,2],
    gammcovrx3 = beta3 >= gammlowerx3 & beta3 <= gammupperx3,
    
    # 6. gam
    gambiasx3 = gam1$coefficients[2] - beta3,
    gamlowerx3 = gam1$coefficients[2] - 1.96*summary(gam1)$se[2],
    gamupperx3 = gam1$coefficients[2] + 1.96*summary(gam1)$se[2],
    gamcovrx3 = beta3 >= gamlowerx3 & beta3 <= gamupperx3
  )
  
  
  rm(list = c("model_brms", "bmamfit", "themam", "bgamfit", "bgammfit",
              "bgamfitted_sample","bgamfitted_summary", "dat","dat2pred", "gam1",
              "gamfit","gammfit","gammpred"))
  
  
  filepath <- paste0("results/","C1/",modeltype)
  if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)
  fileclass <- file.path(filepath,paste0("K",K,"_Nk",Nk,"_beta3",beta3,"_sigma0",sigma0,"_sigma3",sigma3,"_rho",rho))
  
  save_results(paste0(fileclass,"_smooth.csv"),fitteddat)
  save_results(paste0(fileclass,"_randomeffect.csv"),preddat)
  save_results(paste0(fileclass,"_variancecomponents.csv"),varcomp)
  save_results(paste0(fileclass,"_x3.csv"),x3estimate)
  
  rm(list = c("fitteddat", "preddat", "varcomp", "x3estimate"))

}

args <- commandArgs(trailingOnly = TRUE)
tryCatch(
  expr = {
    do_simulation_multiple(K = as.numeric(args[1]), Nk = as.numeric(args[2]), 
                           beta3 = 1,sigma0 = 2,sigma3 = 1,rho=0.5,iter = as.numeric(args[3]),
                           modeltype = as.character(args[4]))
  },
  error = function(e){
    message("error")
  }
)  