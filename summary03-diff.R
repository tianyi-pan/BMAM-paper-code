## Appendix: Simulation Heterogeneity Among Cluster
## Simulation 3: difference between marginal and conditional models


### summarize results ##########
library(dplyr)
library(tidyr)
library(ggplot2)
library(latex2exp)
library(gridExtra)

#### load data ###########
## read and reorganize the data
modeltypelist <- c("marginal","conditional")
folderlist <- c("results")
for (modeltype in modeltypelist) {
  subfolderlist <- c(paste0("C3_diff/",modeltype))
  resulttypes <- c("smooth.csv")
  
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
  
  
  save(smooth, file = paste0("results/C3_diff_",modeltype,".RData"))
  
}


##### plots ###########

f1temp <- function(x){(1.5*sin(pi*x)-2*x)}
f1 <- function(x) f1temp(x)-f1temp(0)
f2temp <- function(x){(1.5*cos(2*pi*(x+1/4)))}
f2 <- function(x) f2temp(x)-f2temp(0)
x <- seq(-1,1,length = 100)

K_choose <- 100
Nk_choose <- 10

theme_set(theme_bw())
theme_replace(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
GGPLOTTEXTSIZE <- 25
GGPLOTSTSIZE <- 24
GGPLOTLEGENDSIZE <- 23

PLOTWIDTH <- 12
PLOTHEIGHT <- 9

cols <- c("bmam"="coral",
          "bgamm"="deepskyblue",
          "f" = "black")

ltys <- c("bmam"="solid",
          "bgamm"="solid",
          "f" = "dotdash")

labs <- c("Marginal",
          "Conditional", 
          "True function")



filename <- "results/C3_diff_marginal.RData"
load(filename)
smooth_id <- smooth
estimationtable_x2 <- smooth_id %>%
  filter(var=="x2", K == K_choose, Nk == Nk_choose, sigma0 == 1)     # choose x2
estimationtable_x2 <- estimationtable_x2 %>% group_by(fx2) %>%
  summarise(bmam = mean(bmamfitted),
            bgamm = mean(bgammfitted))
estimationtable_x2 <- estimationtable_x2[match(round(f2(x), 3), round(estimationtable_x2$fx2, 3)),]
estimationtable_x2$x <- x
estimationtable_x2$f <- f2(x)

df_plot_x2_small <- pivot_longer(estimationtable_x2, cols = c(bmam, bgamm, f),
                                 names_to = "models", values_to = "estimate")


p2_small <- ggplot(data = df_plot_x2_small) + 
  geom_line(aes(x = x, y = estimate, color = models, linetype = models), size = 2.5)+
  scale_color_manual(name = NULL, values = cols, labels = labs)+
  scale_linetype_manual(name = NULL, values = ltys, labels = labs) +
  # ggtitle(TeX('$f^{M}_{2}$ $(X)$')) + 
  ylim(-3,3) +
  theme(text = element_text(size = GGPLOTTEXTSIZE),
        legend.text=element_text(size=GGPLOTLEGENDSIZE),
        plot.subtitle=element_text(size=GGPLOTSTSIZE))+
  xlab("X") + ylab(NULL)

p2_small


estimationtable_x2 <- smooth_id %>%
  filter(var=="x2", K == K_choose, Nk == Nk_choose, sigma0 == 2)    # choose x2
estimationtable_x2 <- estimationtable_x2 %>% group_by(fx2) %>%
  summarise(bmam = mean(bmamfitted),
            bgamm = mean(bgammfitted))
estimationtable_x2 <- estimationtable_x2[match(round(f2(x), 3), round(estimationtable_x2$fx2, 3)),]
estimationtable_x2$x <- x
estimationtable_x2$f <- f2(x)
df_plot_x2_large <- pivot_longer(estimationtable_x2, cols = c(bmam, bgamm, f),
                                 names_to = "models", values_to = "estimate")
p2_large <- ggplot(data = df_plot_x2_large) + 
  geom_line(aes(x = x, y = estimate, color = models, linetype = models), size = 2.5)+
  scale_color_manual(name = NULL, values = cols, labels = labs)+
  scale_linetype_manual(name = NULL, values = ltys, labels = labs) +
  # ggtitle(TeX('$f^{M}_{2}$ $(X)$')) + 
  ylim(-3,3) +
  theme(text = element_text(size = GGPLOTTEXTSIZE),
        legend.text=element_text(size=GGPLOTLEGENDSIZE),
        plot.subtitle=element_text(size=GGPLOTSTSIZE))+
  xlab("X") + ylab(NULL)
p2_large



filepath <- paste0("figures/diff/")
if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)

ggsave(p2_small, filename = paste0(filepath, "diff_small.pdf"), width=PLOTWIDTH, height=PLOTHEIGHT)
ggsave(p2_large, filename = paste0(filepath, "diff_large.pdf"), width=PLOTWIDTH, height=PLOTHEIGHT)

