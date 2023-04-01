## Section 5.4: Application to Beaver Foraging Behaviour

library(brms)
library(bmam)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(gridExtra)
library(stringr)
library(latex2exp)

if(FALSE){
  ## load data
  beaver2 <- readxl::read_xlsx("data/DatasetAnalysis2.xlsx")
  beaver2$Sex0 <- as.numeric(beaver2$Sex=="M"); beaver2$Sex0[beaver2$Sex0==0] <- 0.0000001
  beaver2$Sex1 <- as.numeric(beaver2$Sex=="F"); beaver2$Sex1[beaver2$Sex1==0] <- 0.0000001
  beaver2$Julian_dayF <- beaver2$Julian_day * as.numeric(beaver2$Sex=="F")
  beaver2$Julian_dayM <- beaver2$Julian_day * as.numeric(beaver2$Sex=="M")
  beaver2$Year <- as.factor(beaver2$Year)
  ## brms models
  set.seed(17)
  Tree_brms <- brm(bf(Tree ~ Year + s(Julian_dayF) + s(Julian_dayM) + 
                        (Sex0-1|Beaver)+(Sex1-1|Beaver)),
                   data = beaver2, family = "bernoulli", cores = 4, seed = 17,
                   save_pars = save_pars(all = TRUE),
                   warmup = 1000, iter = 2000, 
                   chains = 4, backend = "cmdstanr")
  Tree_bmam <- bmam(Tree_brms, centered = T)

  Grass_herb_brms <- brm(bf(Grass_herb ~ Year + s(Julian_dayF) + s(Julian_dayM) + 
                              (Sex0-1|Beaver)+(Sex1-1|Beaver)),
                         data = beaver2, family = "bernoulli", cores = 4, seed = 17,
                         save_pars = save_pars(all = TRUE),
                         warmup = 1000, iter = 2000, 
                         chains = 4, backend = "cmdstanr")
  Grass_herb_bmam <- bmam(Grass_herb_brms, center = T)
  save.image(file = "data/application_beavers.RData")
}


load("data/application_beavers.RData")

### draw plot ##########
theme_set(theme_bw())
theme_replace(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
GGPLOTTEXTSIZE <- 20
PLOTWIDTH <- 8
PLOTHEIGHT <- 6


show_gg <- function(gg,title,ylim){
  title_save <- title
  title <- str_replace_all(title,"\\\\"," / ")
  filepath <- paste0("figures/application/")
  if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)
  
  g1 <- gg$Conditional$Comparison[[1]]
  g1 <- g1 + ggtitle(title,subtitle = "Female") +
    xlab("Day of Year") + ylab(TeX('$f_{female}$')) +
    theme(text = element_text(size = GGPLOTTEXTSIZE)) + 
    guides(fill=guide_legend(title="Model"),
           color=guide_legend(title="Model")) +
    # ylim(ylim)+
    scale_y_continuous(limits = ylim, breaks = seq(ylim[1], ylim[2], by = 2))
  ggsave(plot = g1, filename = paste0(filepath, title_save,"_female.pdf"), width=PLOTWIDTH, height=PLOTHEIGHT)
  
  g2 <- gg$Conditional$Comparison[[2]]
  g2 <- g2 + ggtitle(title,subtitle = "Male") + 
    xlab("Day of Year") + ylab(TeX('$f_{male}$')) +
    theme(text = element_text(size = GGPLOTTEXTSIZE)) + 
    guides(fill=guide_legend(title="Model"),
           color=guide_legend(title="Model")) + 
    # ylim(ylim)+
    scale_y_continuous(limits = ylim, breaks = seq(ylim[1], ylim[2], by = 2))

  ggsave(plot = g2, filename = paste0(filepath,title_save,"_male.pdf"), width=PLOTWIDTH, height=PLOTHEIGHT)
  grid.arrange(g1,g2, nrow = 1)
}

show_gg(plot(Grass_herb_bmam, display = FALSE),"Herbs\\Grasses",c(-10,4))
show_gg(plot(Tree_bmam, display = FALSE),"Trees\\Shrubs",c(-4,4))


### Tables ########## 
summary(Grass_herb_bmam,digits = 4)
print(VarCorr(Grass_herb_bmam$Conditional$Brms),digits = 4)

summary(Tree_bmam,digits = 4)
print(VarCorr(Tree_bmam$Conditional$Brms),digits = 4)
