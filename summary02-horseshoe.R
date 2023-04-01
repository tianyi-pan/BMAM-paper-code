## Appendix: Extension marginal variable selection with the horseshoe prior
## Simulation 2: horseshoe simulation

### summarize results ##########
library(dplyr)
library(ggplot2)
library(latex2exp)

theme_set(theme_bw())
GGPLOTTEXTSIZE <- 20
PLOTWIDTH <- 8
PLOTHEIGHT <- 6

### case 1 #########

indept_linear <- read.csv("results/C2_horseshoe/U_K100_corr0.5_p10_sparsity4_indept_linear.csv")
indept_smooth <- read.csv("results/C2_horseshoe/U_K100_corr0.5_p10_sparsity4_indept_smooth.csv")

filepath <- paste0("figures/","C2_horseshoe/")
if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)

## smooth

estimationtable_indept_smooth_x1 <- indept_smooth %>%
  filter(var=="x1") %>%     # choose x1 
  summarize(
    bias_x1_bmam = mean(abs(bmambiasx1)),
    covr_x1_bmam = mean(bmamcovrx1),
    al_x1_bmam = mean(bmamupper - bmamlower),
    
    bias_x1_bmam_flat = mean(abs(bmambiasx1_flat), na.rm = T),
    covr_x1_bmam_flat = mean(bmamcovrx1_flat, na.rm = T),
    al_x1_bmam_flat = mean(bmamupper_flat - bmamlower_flat, na.rm = T),
    
    count = n()/100
  )


re_x1 <- estimationtable_indept_smooth_x1
row.names(re_x1) <- c("independent")
knitr::kable(t(re_x1), digit = 2)

estimationtable_indept_smooth_x2 <- indept_smooth %>%
  filter(var=="x2") %>%     # choose x1 
  summarize(
    bias_x2_bmam = mean(abs(bmambiasx2)),
    covr_x2_bmam = mean(bmamcovrx2),
    al_x2_bmam = mean(bmamupper - bmamlower),
    
    bias_x2_bmam_flat = mean(abs(bmambiasx2_flat)),
    covr_x2_bmam_flat = mean(bmamcovrx2_flat),
    al_x2_bmam_flat = mean(bmamupper_flat - bmamlower_flat),
    
    count = n()/100
  )


re_x2 <- estimationtable_indept_smooth_x2
row.names(re_x2) <- c("independent")
knitr::kable(t(re_x2), digit = 2)


## linear
estimationtable_indept_linear <- indept_linear %>%
  group_by(label) %>%
  summarize(
    TRUE_value = mean(true),
    bias_horseshoe = mean(horseshoebias),
    sd_horseshoe = sd(horseshoebias),
    
    bias_flat = mean(flatbias),
    sd_flat = sd(flatbias)
  )

estimationtable_indept_linear$relative <- (estimationtable_indept_linear$sd_flat/estimationtable_indept_linear$sd_horseshoe)^2
knitr::kable(estimationtable_indept_linear, digit = 2)

plot_df_indept <- data.frame(label = rep(indept_linear$label,2), 
                             true = rep(indept_linear$true, 2),
                             estimate = c(indept_linear$horseshoe, indept_linear$flat),
                             prior = rep(c("horseshoe", "flat"), each = nrow(indept_linear)))
plot_df_indept$label <- factor(plot_df_indept$label, level = c("x3","x4", "x5", "x6", "x7", "x8", "x9", "x10"))



p1 <- ggplot(plot_df_indept, aes(x=label, y=estimate, fill = prior))
p1 <- p1 + geom_boxplot(width = 0.6, alpha = 0.9) + 
  geom_errorbar(aes(ymin = true,ymax = true), width = 0.8, size = 0.8)

p1 <- p1 + theme(text = element_text(size = GGPLOTTEXTSIZE)) + 
  xlab(TeX('$\\beta$')) + ylab("Estimate") + 
  scale_x_discrete(limits = c("x3","x4","x5","x6","x7","x8","x9","x10"),
                   labels=c("x3" = TeX('$X_3$'), "x4" = TeX('$X_4$'), 
                            "x5" = TeX('$X_5$'), "x6" = TeX('$X_6$'), 
                            "x7" = TeX('$X_7$'), "x8" = TeX('$X_8$'),
                            "x9" = TeX('$X_9$'), "x10" = TeX('$X_{10}$')))
p1 <- p1 + ylim(-1, 3.1)
p1
ggsave(p1, filename = paste0(filepath,"indept_linear_sparsity4.pdf"), width = PLOTWIDTH, height = PLOTHEIGHT)



### case 2 ####
indept_linear <- read.csv("results/C2_horseshoe/U_K100_corr0.5_p10_sparsity1_indept_linear.csv")
indept_smooth <- read.csv("results/C2_horseshoe/U_K100_corr0.5_p10_sparsity1_indept_smooth.csv")
## smooth
estimationtable_indept_smooth_x1 <- indept_smooth %>%
  filter(var=="x1") %>%     # choose x1 
  summarize(
    bias_x1_bmam = mean(abs(bmambiasx1)),
    covr_x1_bmam = mean(bmamcovrx1),
    al_x1_bmam = mean(bmamupper - bmamlower),
    
    bias_x1_bmam_flat = mean(abs(bmambiasx1_flat)),
    covr_x1_bmam_flat = mean(bmamcovrx1_flat),
    al_x1_bmam_flat = mean(bmamupper_flat - bmamlower_flat, na.rm = T),
    
    count = n()/100
  )


re_x1 <- estimationtable_indept_smooth_x1
row.names(re_x1) <- c("independent")

re_x1 <- estimationtable_indept_smooth_x1
knitr::kable(t(re_x1), digit = 2)



estimationtable_indept_smooth_x2 <- indept_smooth %>%
  filter(var=="x2") %>%     # choose x1 
  summarize(
    bias_x2_bmam = mean(abs(bmambiasx2)),
    covr_x2_bmam = mean(bmamcovrx2),
    al_x2_bmam = mean(bmamupper - bmamlower),
    
    bias_x2_bmam_flat = mean(abs(bmambiasx2_flat)),
    covr_x2_bmam_flat = mean(bmamcovrx2_flat),
    al_x2_bmam_flat = mean(bmamupper_flat - bmamlower_flat),
    
    count = n()/100
  )

re_x2 <- estimationtable_indept_smooth_x2
row.names(re_x2) <- c("independent")
re_x2 <- estimationtable_indept_smooth_x2
knitr::kable(t(re_x2), digit = 2)


## linear
estimationtable_indept_linear <- indept_linear %>%
  group_by(label) %>%
  summarize(
    TRUE_value = mean(true),
    bias_horseshoe = mean(horseshoebias),
    sd_horseshoe = sd(horseshoebias),
    
    bias_flat = mean(flatbias),
    sd_flat = sd(flatbias)
  )

estimationtable_indept_linear$rel <- (estimationtable_indept_linear$sd_flat/estimationtable_indept_linear$sd_horseshoe)^2
knitr::kable(estimationtable_indept_linear, digit = 2)

plot_df_indept <- data.frame(label = rep(indept_linear$label,2), 
                      true = rep(indept_linear$true, 2),
                      estimate = c(indept_linear$horseshoe, indept_linear$flat),
                      prior = rep(c("horseshoe", "flat"), each = nrow(indept_linear)))
plot_df_indept$label <- factor(plot_df_indept$label, level = c("x3","x4", "x5", "x6", "x7", "x8", "x9", "x10"))


p1 <- ggplot(plot_df_indept, aes(x=label, y=estimate, fill = prior))
p1 <- p1 + geom_boxplot(width = 0.6, alpha = 0.9) + 
  geom_errorbar(aes(ymin = true,ymax = true), width = 0.8, size = 0.8)

p1 <- p1 + theme(text = element_text(size = GGPLOTTEXTSIZE)) + 
  xlab(TeX('$\\beta$')) + ylab("Estimate") + 
  scale_x_discrete(labels=c("x3" = TeX('$X_3$'), "x4" = TeX('$X_4$'), 
                            "x5" = TeX('$X_5$'), "x6" = TeX('$X_6$'), 
                            "x7" = TeX('$X_7$'), "x8" = TeX('$X_8$'),
                            "x9" = TeX('$X_9$'), "x10" = TeX('$X_{10}$')))
p1 <- p1 + ylim(-1, 3.1)
p1
ggsave(p1, filename = paste0(filepath, "indept_linear_sparsity1.pdf"), width = PLOTWIDTH, height = PLOTHEIGHT)
