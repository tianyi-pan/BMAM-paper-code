# Marginal and Conditional Inference for Bayesian Hierarchical Additive Models in Ecology
Source Code to replicate the results in the paper *Marginal and Conditional Inference for Bayesian Hierarchical Additive Models in Ecology*

## Packages
[BMAM package](https://github.com/tianyi-pan/BMAM): Bayesian Marginal Additive Model (BMAM) based on [brms](https://github.com/paul-buerkner/brms) and [brmsmargins](https://github.com/JWiley/brmsmargins). 

+ Install the development version from Github:
```R
# install.packages('devtools')
devtools::install_github('tianyi-pan/BMAM')
```
+ A brief [vignette](https://tianyi-pan.github.io/BMAM)

## Simulation
Each part contains three R code files, including 1) a function for the simulation 2) a script for running, and 3) a script for summarizing the results.

+ Section 3.4 & Section 5.3
    - do_simu_main.R
    - simulation01-main.R
    - summary01-main.R
+ Appendix D:
    - do_simu_horseshoe.R
    - simulation02-horseshoe.R
    - summary02-horseshoe.R
+ Appendix H:
    - do_simu_diff.R
    - simulation03-diff.R
    - summary03-diff.R

The following script is for running the simulation and summarizing the results. 
+ Appendix A:
    - ---.R



## Application
+ Section 2: beavers_summary.R
+ Section 5.4: application.R

All the results will be saved in the current directory(`pwd` in Shell; `getwd()` in R). The scripts will create directories ./results and ./figures for the complete results and figures.