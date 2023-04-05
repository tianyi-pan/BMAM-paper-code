## Section 5.4: Application to Beaver Foraging Behaviour

library(brms)
library(bmam)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(gridExtra)
library(stringr)
library(latex2exp)

if(TRUE){
  ## load data
  beaver2 <- readxl::read_xlsx("data/DatasetAnalysis2.xlsx")
  beaver2$Sex0 <- as.numeric(beaver2$Sex=="M"); beaver2$Sex0[beaver2$Sex0==0] <- 0.0000001
  beaver2$Sex1 <- as.numeric(beaver2$Sex=="F"); beaver2$Sex1[beaver2$Sex1==0] <- 0.0000001
  beaver2$Julian_dayF <- beaver2$Julian_day * as.numeric(beaver2$Sex=="F")
  beaver2$Julian_dayM <- beaver2$Julian_day * as.numeric(beaver2$Sex=="M")
  beaver2$Year <- as.factor(beaver2$Year)
  beaver2$Sex <- as.factor(beaver2$Sex)
  Preddat <- data.frame(Julian_dayF = c(seq(min(beaver2$Julian_dayF[which(beaver2$Sex=="F")]),max(beaver2$Julian_dayF), length.out = 100),
                                       rep(0,100)),
                        Julian_dayM = c(rep(0,100), 
                                        seq(min(beaver2$Julian_dayM[which(beaver2$Sex=="M")]),max(beaver2$Julian_dayM), length.out = 100)),
                        Julian_day = c(seq(min(beaver2$Julian_dayF[which(beaver2$Sex=="F")]),max(beaver2$Julian_dayF), length.out = 100),
                                       seq(min(beaver2$Julian_dayM[which(beaver2$Sex=="M")]),max(beaver2$Julian_dayM), length.out = 100)),
                        Sex = rep(c(beaver2$Sex[which(beaver2$Sex == "F")[1]],
                                    beaver2$Sex[which(beaver2$Sex == "M")[1]]),
                                  each = 100),
                        varname = rep(c("Julian_dayF", "Julian_dayM"), each = 100),
                        Year = "2006")
  
  set.seed(17)
  Tree_brms <- brm(bf(Tree ~ Year + Sex + s(Julian_day, by = Sex) + 
                        (Sex0-1|Beaver)+(Sex1-1|Beaver)),
                   data = beaver2, family = "bernoulli", cores = 4, seed = 17,
                   save_pars = save_pars(all = TRUE),
                   warmup = 1000, iter = 2000, 
                   chains = 4, backend = "cmdstanr")
  Tree_bmam <- bmam(Tree_brms, centered = T, preddat = Preddat)
  
  Grass_herb_brms <- brm(bf(Grass_herb ~ Year + Sex + s(Julian_day, by = Sex) + 
                              (Sex0-1|Beaver)+(Sex1-1|Beaver)),
                         data = beaver2, family = "bernoulli", cores = 4, seed = 17,
                         save_pars = save_pars(all = TRUE),
                         warmup = 1000, iter = 2000, 
                         chains = 4, backend = "cmdstanr")
  Grass_herb_bmam <- bmam(Grass_herb_brms, center = T, preddat = Preddat)

  save.image(file = "data/application_beavers_sex.RData")
}


load("data/application_beavers_sex.RData")

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
  ggsave(plot = g1, filename = paste0(filepath, title_save,"_sex_female.pdf"), width=PLOTWIDTH, height=PLOTHEIGHT)
  
  g2 <- gg$Conditional$Comparison[[2]]
  g2 <- g2 + ggtitle(title,subtitle = "Male") + 
    xlab("Day of Year") + ylab(TeX('$f_{male}$')) +
    theme(text = element_text(size = GGPLOTTEXTSIZE)) + 
    guides(fill=guide_legend(title="Model"),
           color=guide_legend(title="Model")) + 
    # ylim(ylim)+
    scale_y_continuous(limits = ylim, breaks = seq(ylim[1], ylim[2], by = 2))

  ggsave(plot = g2, filename = paste0(filepath,title_save,"_sex_male.pdf"), width=PLOTWIDTH, height=PLOTHEIGHT)
  grid.arrange(g1,g2, nrow = 1)
}

show_gg(plot(Grass_herb_bmam, display = FALSE),"Herbs\\Grasses",c(-14,4))
show_gg(plot(Tree_bmam, display = FALSE),"Trees\\Shrubs",c(-3,6))


### Tables ########## 
summary(Grass_herb_bmam,digits = 4)
print(VarCorr(Grass_herb_bmam$Conditional$Brms),digits = 4)

summary(Tree_bmam,digits = 4)
print(VarCorr(Tree_bmam$Conditional$Brms),digits = 4)
