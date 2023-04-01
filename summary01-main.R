## Section 3.4 & Section 5.3
## Simulation 1: main simulation

### summarize results ##########
library(dplyr)
library(tidyr)
library(ggplot2)
library(latex2exp)


#### load data ########
## read and reorganize the data
modeltypelist <- c("conditional","marginal")
folderlist <- c("results")
for (modeltype in modeltypelist) {
  subfolderlist <- c(paste0("C1/",modeltype))
  resulttypes <- c("smooth.csv", "randomeffect.csv", "variancecomponents.csv","x3.csv")
  
  getresults <- function(folder, subfolder, resulttype){
    file_list <- list.files(paste(folder, subfolder, sep = "/"), full.names = T)
    target_list <- file_list[grepl(resulttype, file_list)]
    results_list <- lapply(target_list, read.csv)
    results <- do.call("rbind", results_list)
    return(results)
  }
  
  lst <- expand.grid(folderlist, subfolderlist,stringsAsFactors = F)
  lsttodo <- vector(mode='list',length=nrow(lst))
  for (i in 1:length(lsttodo)) lsttodo[[i]] <- lst[i, ]
  
  
  smooth <- lapply(lsttodo, function(tmp){
    getresults(tmp[1], tmp[2], "smooth.csv")  
  }) %>% do.call(what = "rbind")
  
  randomeffect <- lapply(lsttodo, function(tmp){
    getresults(tmp[1], tmp[2], "randomeffect.csv")  
  }) %>% do.call(what = "rbind")
  
  variancecomponents <- lapply(lsttodo, function(tmp){
    getresults(tmp[1], tmp[2], "variancecomponents.csv")  
  }) %>% do.call(what = "rbind")
  
  
  x3 <- lapply(lsttodo, function(tmp){
    getresults(tmp[1], tmp[2], "x3.csv")  
  }) %>% do.call(what = "rbind")
  
  save(smooth,randomeffect,variancecomponents, x3, file = paste0("results/C1_",modeltype,".RData"))
  
}


#### Tables ##########

for (modeltype in c("conditional", "marginal")) {
  cat(modeltype)
  filename <- paste0("results/C1_", modeltype, ".RData") # cluster level
  load(filename)
  Nk_value <- c(10,20)
  K_list <- c(300,500)
  
  smooth_id <- smooth
  
  ### x1 ##############
  value_ifmam <- 1
  estimationtable_x1 <- smooth_id %>%
    filter(var=="x1", K %in% K_list) %>%     # choose x1 
    filter(ifmam == value_ifmam) %>% # choose mam 
    group_by(K,Nk,sigma0,sigma3,rho) %>%
    summarize(
      bias_x1_bmam = mean(abs(bmambiasx1)),
      covr_x1_bmam = mean(abs(bmamcovrx1)),
      
      bias_x1_bgamm = mean(abs(bgammbiasx1)),
      covr_x1_bgamm = mean(bgammcovrx1),
      
      bias_x1_bgam = mean(abs(bgambiasx1)),
      covr_x1_bgam = mean(bgamcovrx1),
      
      bias_x1_mam = mean(abs(mammargbiasx1), na.rm = TRUE),
      covr_x1_mam = mean(mammargcovrx1, na.rm = TRUE),
      
      bias_x1_mam_cond = mean(abs(mamcondbiasx1)),
      covr_x1_mam_cond = mean(mamcondcovrx1),
      
      bias_x1_gamm = mean(abs(gammcondbiasx1)),
      covr_x1_gamm = mean(gammcondcovrx1),
      
      bias_x1_gam = mean(abs(gammargbiasx1)),
      covr_x1_gam = mean(gammargcovrx1),
      
      count = n()/100
    )
  failure <- select(estimationtable_x1, K, Nk, count) %>% mutate(failure = (1-count/500)*100)
  print(knitr::kable(failure, digits = 2))
  df_plot_bias_x1 <- filter(estimationtable_x1, Nk %in% Nk_value) %>% 
    pivot_longer(., cols = c(bias_x1_bmam, bias_x1_bgamm,bias_x1_bgam,bias_x1_mam,bias_x1_mam_cond,bias_x1_gamm,bias_x1_gam), # bias_x1_gees
                 names_to = "Methods", values_to = "Bias") %>%
    select(K, Nk, Methods, Bias)
  df_plot_covr_x1 <- filter(estimationtable_x1, Nk %in% Nk_value) %>% 
    pivot_longer(., cols = c(covr_x1_bmam, covr_x1_bgamm,covr_x1_bgam,covr_x1_mam,covr_x1_mam_cond,covr_x1_gamm,covr_x1_gam),
                 names_to = "Methods", values_to = "CoverageRate")%>%
    select(K, Nk, Methods, CoverageRate)
  
  cat("x1: bias")
  print(knitr::kable(df_plot_bias_x1, digits = 2))
  cat("x1: coverage rate")
  print(knitr::kable(df_plot_covr_x1, digits = 2))
  
  ### x2 ##############
  value_ifmam <- 1
  estimationtable_x2 <- smooth_id %>%
    filter(var=="x2", K %in% K_list) %>%     # choose x2 
    filter(ifmam == value_ifmam) %>% # choose mam
    group_by(K,Nk,sigma0,sigma3,rho) %>%
    summarize(
      bias_x2_bmam = mean(abs(bmambiasx2)),
      covr_x2_bmam = mean(abs(bmamcovrx2)),
      
      bias_x2_bgamm = mean(abs(bgammbiasx2)),
      covr_x2_bgamm = mean(bgammcovrx2),
      
      bias_x2_bgam = mean(abs(bgambiasx2)),
      covr_x2_bgam = mean(bgamcovrx2),
      
      bias_x2_mam = mean(abs(mammargbiasx2), na.rm = TRUE),
      covr_x2_mam = mean(mammargcovrx2, na.rm = TRUE),
      
      bias_x2_mam_cond = mean(abs(mamcondbiasx2)),
      covr_x2_mam_cond = mean(mamcondcovrx2),
      
      bias_x2_gamm = mean(abs(gammcondbiasx2)),
      covr_x2_gamm = mean(gammcondcovrx2),
      
      bias_x2_gam = mean(abs(gammargbiasx2)),
      covr_x2_gam = mean(gammargcovrx2),
      
      count = n()/100
    )
  df_plot_bias_x2 <- filter(estimationtable_x2, Nk %in% Nk_value) %>% 
    pivot_longer(., cols = c(bias_x2_bmam, bias_x2_bgamm,bias_x2_bgam,bias_x2_mam,bias_x2_mam_cond,bias_x2_gamm,bias_x2_gam), # bias_x2_gees
                 names_to = "Methods", values_to = "Bias") %>%
    select(K, Nk, Methods, Bias)
  df_plot_covr_x2 <- filter(estimationtable_x2, Nk %in% Nk_value) %>% 
    pivot_longer(., cols = c(covr_x2_bmam, covr_x2_bgamm,covr_x2_bgam,covr_x2_mam,covr_x2_mam_cond,covr_x2_gamm,covr_x2_gam),
                 names_to = "Methods", values_to = "CoverageRate")%>%
    select(K, Nk, Methods, CoverageRate)
  cat("x2: bias")
  print(knitr::kable(df_plot_bias_x2, digits = 2))
  cat("x2: coverage rate")
  print(knitr::kable(df_plot_covr_x2, digits = 2))

  ### x3 ##############
  value_ifmam <- 1
  estimationtable_x3 <- x3 %>%
    filter(K %in% K_list) %>%     # choose x2 
    filter(ifmam == value_ifmam) %>% # choose mam
    group_by(K,Nk,sigma0,sigma3,rho) %>%
    summarize(
      bias_x3_bmam = mean(abs(bmambiasx3)),
      covr_x3_bmam = mean(abs(bmamcovrx3)),
      
      bias_x3_bgamm = mean(abs(bgammbiasx3)),
      covr_x3_bgamm = mean(bgammcovrx3),
      
      bias_x3_bgam = mean(abs(bgambiasx3)),
      covr_x3_bgam = mean(bgamcovrx3),
      
      bias_x3_mam = mean(abs(mambiasx3), na.rm = TRUE),
      covr_x3_mam = mean(mamcovrx3, na.rm = TRUE),
      
      bias_x3_gamm = mean(abs(gammbiasx3)),
      covr_x3_gamm = mean(gammcovrx3),
      
      bias_x3_gam = mean(abs(gambiasx3)),
      covr_x3_gam = mean(gamcovrx3)
    )
  df_plot_bias_x3 <- filter(estimationtable_x3, Nk %in% Nk_value) %>% 
    pivot_longer(., cols = c(bias_x3_bmam, bias_x3_bgamm,bias_x3_bgam,bias_x3_mam,bias_x3_gamm,bias_x3_gam), # bias_x2_gees
                 names_to = "Methods", values_to = "Bias") %>%
    select(K, Nk, Methods, Bias)
  df_plot_covr_x3 <- filter(estimationtable_x3, Nk %in% Nk_value) %>% 
    pivot_longer(., cols = c(covr_x3_bmam, covr_x3_bgamm,covr_x3_bgam,covr_x3_mam,covr_x3_gamm,covr_x3_gam),
                 names_to = "Methods", values_to = "CoverageRate")%>%
    select(K, Nk, Methods, CoverageRate)
  cat("x3: bias")
  print(knitr::kable(df_plot_bias_x3, digits = 2))
  cat("x3: coverage rate")
  print(knitr::kable(df_plot_covr_x3, digits = 2))
  
  
  ### variance component ###########
  vc_bias <- variancecomponents %>%
    filter(K %in% K_list, Nk %in% Nk_value) %>%     
    filter(ifmam == value_ifmam) %>% # choose mam 
    group_by(K,Nk,sigma0,sigma3,rho) %>%
    summarize(
      mam_sd_int_bias = mean(mam_sd_int - sigma0),
      mam_sd_slope_bias = mean(mam_sd_slope - sigma3),
      mam_cor_bias = mean(mam_cor - rho),
      
      gamm_sd_int_bias = mean(gamm_sd_int - sigma0),
      gamm_sd_slope_bias = mean(gamm_sd_slope - sigma3),
      gamm_cor_bias = mean(gamm_cor - rho),
      
      bgamm_sd_int_bias = mean(bgamm_sd_int - sigma0),
      bgamm_sd_slope_bias = mean(bgamm_sd_slope - sigma3),
      bgamm_cor_bias = mean(bgamm_cor - rho)
    )
  cat("variance component")
  knitr::kable(t(vc_bias), digits = 2)
}





### Section 3.1 Figures ############
library(dplyr)
library(tidyr)
library(ggplot2)
library(latex2exp)
library(gridExtra)


f1temp <- function(x){(1.5*sin(pi*x)-2*x)}
f1 <- function(x) f1temp(x)-f1temp(0)
f2temp <- function(x){(1.5*cos(2*pi*(x+1/4)))}
f2 <- function(x) f2temp(x)-f2temp(0)
x <- seq(-1,1,length = 100)

value_ifmam <- 1
K_choose <- 500
Nk_choose <- 10

theme_set(theme_bw())
theme_replace(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
GGPLOTTEXTSIZE <- 25
GGPLOTSTSIZE <- 24
GGPLOTLEGENDSIZE <- 23

PLOTWIDTH <- 12
PLOTHEIGHT <- 9

cols <- c("gamm"="deepskyblue",
          "gam"="coral",
          "f" = "black")

ltys <- c("gamm"="solid",
          "gam"="solid",
          "f" = "dashed")

labs <- c("GAMM",
          "GAM",
          "True function")

alphas <- c("gamm"=0.8,
            "gam"=0.8,
            "f" = 1)

sizes <- c("gamm"=2.5,
           "gam"=2.5,
           "f" = 2.5)

### marginal models #########
filename <- "results/C1_marginal.RData"

load(filename)
smooth_id <- smooth
estimationtable_x1 <- smooth_id %>%
  filter(var=="x1", K == K_choose, Nk == Nk_choose) %>%    
  filter(ifmam == value_ifmam)
# estimationtable_x1$f_test <- round(estimationtable_x1$bgamfitted - estimationtable_x1$bgambiasx1, 3)
estimationtable_x1 <- estimationtable_x1 %>% group_by(fx1) %>%
  summarise(gamm = mean(gammcondfitted),
            gam = mean(gammargfitted))
estimationtable_x1 <- estimationtable_x1[match(round(f1(x), 3), round(estimationtable_x1$fx1, 3)),]
estimationtable_x1$x <- x
estimationtable_x1$f <- f1(x)


df_plot_x1_small <- pivot_longer(estimationtable_x1, cols = c(f, gamm, gam),
                                 names_to = "models", values_to = "estimate")
df_plot_x1_small$models <- factor(df_plot_x1_small$models, levels = c("gamm","gam","f"), labels = c("gamm", "gam", "f"))

p1_marg <- ggplot(data = df_plot_x1_small) + 
  geom_line(aes(x = x, y = estimate, color = models, linetype = models, alpha = models, size = models))+
  scale_color_manual(name = NULL, values = cols, labels = labs)+
  scale_alpha_manual(name = NULL, values = alphas, labels = labs)+
  scale_linetype_manual(name = NULL, values = ltys, labels = labs) +
  scale_size_manual(name = NULL, values = sizes, labels = labs) +
  # ggtitle(TeX('$f^{M}_{1}$ $(X)$')) + 
  ylim(-3,3) +
  theme(text = element_text(size = GGPLOTTEXTSIZE),
        legend.text=element_text(size=GGPLOTLEGENDSIZE),
        plot.subtitle=element_text(size=GGPLOTSTSIZE))+
  xlab("X") + ylab(NULL)

p1_marg



smooth_id <- smooth
estimationtable_x2 <- smooth_id %>%
  filter(var=="x2", K == K_choose, Nk == Nk_choose) %>%     # choose x2
  filter(ifmam == value_ifmam)
estimationtable_x2 <- estimationtable_x2 %>% group_by(fx2) %>%
  summarise(gamm = mean(gammcondfitted),
            gam = mean(gammargfitted))
estimationtable_x2 <- estimationtable_x2[match(round(f2(x), 3), round(estimationtable_x2$fx2, 3)),]
estimationtable_x2$x <- x
estimationtable_x2$f <- f2(x)

df_plot_x2_small <- pivot_longer(estimationtable_x2, cols = c(gamm, gam, f),
                                 names_to = "models", values_to = "estimate")
df_plot_x2_small$models <- factor(df_plot_x2_small$models, levels = c("gamm","gam","f"), labels = c("gamm", "gam", "f"))


p2_marg <- ggplot(data = df_plot_x2_small) + 
  geom_line(aes(x = x, y = estimate, color = models, linetype = models, alpha = models, size = models))+
  scale_color_manual(name = NULL, values = cols, labels = labs)+
  scale_alpha_manual(name = NULL, values = alphas, labels = labs)+
  scale_size_manual(name = NULL, values = sizes, labels = labs) +
  scale_linetype_manual(name = NULL, values = ltys, labels = labs) +
  # ggtitle(TeX('$f^{M}_{2}$ $(X)$')) + 
  ylim(-3,3) +
  theme(text = element_text(size = GGPLOTTEXTSIZE),
        legend.text=element_text(size=GGPLOTLEGENDSIZE),
        plot.subtitle=element_text(size=GGPLOTSTSIZE))+
  xlab("X") + ylab(NULL)

p2_marg




### conditional models #########
filename <- "results/C1_conditional.RData"

load(filename)
smooth_id <- smooth
estimationtable_x1 <- smooth_id %>%
  filter(var=="x1", K == K_choose, Nk == Nk_choose) %>%    
  filter(ifmam == value_ifmam)
estimationtable_x1 <- estimationtable_x1 %>% group_by(fx1) %>%
  summarise(gamm = mean(gammcondfitted),
            gam = mean(gammargfitted))
estimationtable_x1 <- estimationtable_x1[match(round(f1(x), 3), round(estimationtable_x1$fx1, 3)),]
estimationtable_x1$x <- x
estimationtable_x1$f <- f1(x)

df_plot_x1_small <- pivot_longer(estimationtable_x1, cols = c(gamm, gam, f),
                                 names_to = "models", values_to = "estimate")
df_plot_x1_small$models <- factor(df_plot_x1_small$models, levels = c("gamm","gam","f"), labels = c("gamm", "gam", "f"))


p1_cond <- ggplot(data = df_plot_x1_small) + 
  geom_line(aes(x = x, y = estimate, color = models, linetype = models, alpha = models, size = models))+
  scale_color_manual(name = NULL, values = cols, labels = labs)+
  scale_alpha_manual(name = NULL, values = alphas, labels = labs)+
  scale_linetype_manual(name = NULL, values = ltys, labels = labs) +
  scale_size_manual(name = NULL, values = sizes, labels = labs) +
  # ggtitle(TeX('$f^{C}_{1}$ $(X)$')) + 
  ylim(-3,3) +
  theme(text = element_text(size = GGPLOTTEXTSIZE),
        legend.text=element_text(size=GGPLOTLEGENDSIZE),
        plot.subtitle=element_text(size=GGPLOTSTSIZE))+
  xlab("X") + ylab(NULL)

p1_cond



smooth_id <- smooth
estimationtable_x2 <- smooth_id %>%
  filter(var=="x2", K == K_choose, Nk == Nk_choose) %>%     # choose x2
  filter(ifmam == value_ifmam)
estimationtable_x2 <- estimationtable_x2 %>% group_by(fx2) %>%
  summarise(gamm = mean(gammcondfitted),
            gam = mean(gammargfitted))
estimationtable_x2 <- estimationtable_x2[match(round(f2(x), 3), round(estimationtable_x2$fx2, 3)),]
estimationtable_x2$x <- x
estimationtable_x2$f <- f2(x)

df_plot_x2_small <- pivot_longer(estimationtable_x2, cols = c(gamm, gam, f),
                                 names_to = "models", values_to = "estimate")
df_plot_x2_small$models <- factor(df_plot_x2_small$models, levels = c("gamm","gam","f"), labels = c("gamm", "gam", "f"))


p2_cond <- ggplot(data = df_plot_x2_small) + 
  geom_line(aes(x = x, y = estimate, color = models, linetype = models, alpha = models, size = models))+
  scale_color_manual(name = NULL, values = cols, labels = labs)+
  scale_alpha_manual(name = NULL, values = alphas, labels = labs)+
  scale_linetype_manual(name = NULL, values = ltys, labels = labs) +
  scale_size_manual(name = NULL, values = sizes, labels = labs) +
  # ggtitle(TeX('$f^{C}_{2}$ $(X)$')) + 
  ylim(-3,3) +
  theme(text = element_text(size = GGPLOTTEXTSIZE),
        legend.text=element_text(size=GGPLOTLEGENDSIZE),
        plot.subtitle=element_text(size=GGPLOTSTSIZE))+
  xlab("X") + ylab(NULL)

p2_cond

filepath <- paste0("figures/S3.1/")
if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)

ggsave(p1_marg, filename = "figures/S3.1/f1_marg.pdf", width=PLOTWIDTH, height=PLOTHEIGHT)
ggsave(p2_marg, filename = "figures/S3.1/f2_marg.pdf", width=PLOTWIDTH, height=PLOTHEIGHT)
ggsave(p1_cond, filename = "figures/S3.1/f1_cond.pdf", width=PLOTWIDTH, height=PLOTHEIGHT)
ggsave(p2_cond, filename = "figures/S3.1/f2_cond.pdf", width=PLOTWIDTH, height=PLOTHEIGHT)



