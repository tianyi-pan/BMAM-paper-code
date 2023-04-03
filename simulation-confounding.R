## Confounding Simulation
### Countering the omitted variable argument of Muff et al (2016)
### Treating z_k as a confounder (or precision variable)

## load packages
library(lme4)
library(geepack)
library(ggplot2)
library(tidyverse)

## inverse of logit function
expit <- function(x){
  exp(x)/(1+exp(x))
}

set.seed(1234)
K <- 400 ## number of clusters
Nk <- 10 ## cluster size
RR <- 400 ## repetitions
confounding <- TRUE ## set to FALSE to run "precision variable" simulation

ests <- c()
for(rr in 1:RR){
  
  id <- rep(seq(1,K),each=Nk) ## cluster ids
  
  ## generate z as confounder or "precision variable"
  if(confounding==TRUE){
    # ## z_k as confounder; x_k is unit level
    # z <- rep(rnorm(K,0,1),each=Nk) 
    # x <- rnorm(K*Nk,0.5*z,1) ## this version can be corrected via poor mans version
    # x <- rbinom(K*Nk,1,expit(-1+(2*z)^2))  ## this one cannot be corrected even with poor mans version
    
    ## z_k as confounder; x_k is cluster level
    ### yields even worse bias for GLMM
    z_k <- rnorm(K,0,1)
    z <- rep(z_k,each=Nk)
    x <- rep(rnorm(K,0.5*z_k,1),each=Nk)
  }else{
    ## z_k as "precision variable"
    z <- rep(rnorm(K,0,1),each=Nk)
    x <- rnorm(K*Nk,0,1)
  }
  x2 <- rnorm(K*Nk,0,1)
  
  # ## z_k as mediator
  # xk <- rnorm(K,0,1)
  # x <- rep(xk,each=Nk)
  # z <- rep(rnorm(K,xk,1),each=Nk) 
  
  

  ## generate outcomes from GLM
  y <- rbinom(K*Nk,1,expit(-1+x+z+x2))
  
  ## fit models
  adjusted <- glm(y~x+x2+z,family=binomial) ## adjusting for z_k
  omit <- geeglm(y~x+x2,id=id,family=binomial,corstr="independence")   ## omitting z_k
  glmm <- glmer(y~x+x2+(1|id),family=binomial) ## GLMM with random intercept instead of z_k
  
  ests <- rbind(ests,c(adjusted$coef[2],omit$coef[2],summary(glmm)$coef[2,1]))
  
  print(rr)
}
colnames(ests) <- c("Adjusted","Marginal","GLMM")

## summarize estimates of logOR for X
round(apply(ests,2,mean),2)
round(apply(ests,2,sd),3)

## boxplot
boxplot(ests);abline(h=1,lty=2)

## make ggplot for paper
ests_tab <- data.frame(ests) %>% 
  pivot_longer(cols=everything(),
               names_to = "Model", 
               values_to = "Estimates")
ests_tab$Model <- factor(ests_tab$Model,levels=c("Adjusted","Marginal","GLMM"))


ests_box <- ggplot(ests_tab,aes(y=Estimates,x=Model))+
  geom_boxplot(fill="grey") +
  ylim(c(1-0.6,1+0.6))+
  geom_hline(yintercept=1,linetype = 'dashed',alpha=0.65)+
  labs(x="", y="Estimate")+
  theme_classic() 
ests_box


filepath <- paste0("figures/confounding/")
if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)

if(confounding==TRUE){
  ests_box <- ests_box+ggtitle(expression(paste(z[k]," as confounder")))
  ggsave("figures/confounding/box_confounding.pdf",ests_box,width=4,height=4)
}else{
  ests_box <- ests_box+ggtitle(expression(paste(z[k]," as precision variable")))
  ggsave("figures/confounding/box_precision.pdf",ests_box,width=4,height=4)
}

