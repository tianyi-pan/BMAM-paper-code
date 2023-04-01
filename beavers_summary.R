## Section 2: Case Study Sex-Specific Trends in Beaver Foraging Behaviour

### import data #########
library(dplyr)
beaver2 <- readxl::read_xlsx("data/DatasetAnalysis2.xlsx")
beaver2$week <- as.factor(strftime(beaver2$Date, format = "%V"))


library(ggplot2)
theme_set(theme_bw())
theme_replace(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
GGPLOTTEXTSIZE <- 25
GGPLOTSTSIZE <- 24
GGPLOTLEGENDSIZE <- 23

PLOTWIDTH <- 12
PLOTHEIGHT <- 9


### Tables ##########
## observations
food <- function(data){
  print("overall")
  print(nrow(data))
  print(round(nrow(data)/3311, digits = 2))
  print("tree")
  print(c(round(sum(data$Tree), digits = 0), round(mean(data$Tree), digits = 2)))
  print("grass")
  print(c(round(sum(data$Grass_herb), digits = 0),round(mean(data$Grass_herb), digits = 2)))
  print("aquatic")
  print(c(round(sum(data$Aquatic), digits = 0),round(mean(data$Aquatic), digits = 2)))
}
time <- function(data){
  print("*****2000")
  a <- filter(data, Year == "2000") %>% nrow()
  print(a)
  print(round(a/nrow(data), digits = 2))
  
  print("*****2006")
  a <- filter(data, Year == "2006")%>% nrow()
  print(a)
  print(round(a/nrow(data), digits = 2))
  
  
  print("*****2007")
  a <- filter(data, Year == "2007")%>% nrow()
  print(a)
  print(round(a/nrow(data), digits = 2))
  
  
  print("*****3,4")
  a <- filter(data, Month %in% c("3","4")) %>% nrow() 
  print(a)
  print(round(a/nrow(data), digits = 2))
  
  
  print("*****5,6")
  a <- filter(data, Month %in% c("5","6")) %>% nrow() 
  print(a)
  print(round(a/nrow(data), digits = 2))
  
  print("*****7,8")
  a <- filter(data, Month %in% c("7","8")) %>% nrow() 
  print(a)
  print(round(a/nrow(data), digits = 2))
  
}

food(beaver2)
time(beaver2)

filter(beaver2, beaver2$Sex == "F") %>% food()
filter(beaver2, beaver2$Sex == "F") %>% time()

filter(beaver2, beaver2$Sex == "M") %>% food()
filter(beaver2, beaver2$Sex == "M") %>% time()




### plot1 ##################
cols <- c("Trees/Shrubs"="grey10",
          "Herbs/Grasses"="coral",
          "Aquatic Vegetation" = "deepskyblue")

fills <- c("Trees/Shrubs"="grey",
          "Herbs/Grasses"="coral",
          "Aquatic Vegetation" = "deepskyblue")


beavers_M <- filter(beaver2, beaver2$Sex == "M")
data <- beavers_M  %>%
  group_by(Month, Forage_Sub_Cat) %>%
  summarise(n = sum(Tree+Aquatic+Grass_herb)) %>%
  mutate(percentage = n / sum(n)) %>% as_tibble()
data <- add_row(data, Month = 3, Forage_Sub_Cat = "Grass_herb", n = 0,percentage = 0)
data$Forage_Sub_Cat <- factor(data$Forage_Sub_Cat,
                              levels = c("Tree","Grass_herb","Aquatic"), 
                              labels = c("Trees/Shrubs", "Herbs/Grasses", "Aquatic Vegetation"))


p_M_percent <- ggplot(data, aes(x=Month, y=percentage, color=Forage_Sub_Cat)) + 
  geom_line(size = 2, alpha = 0.9) + 
  geom_point(size = 5) + 
  ggtitle(NULL) + 
  ylim(0,1) +
  theme(text = element_text(size = GGPLOTTEXTSIZE),
        legend.text=element_text(size=GGPLOTLEGENDSIZE),
        plot.subtitle=element_text(size=GGPLOTSTSIZE),
        legend.spacing.y = unit(0.3, 'cm') )+
  guides(fill = guide_legend(byrow = TRUE))+
  scale_color_manual(name = NULL, values = cols)+
  xlab("Month") + ylab("Percentage")



p_M_percent

beavers_F <- filter(beaver2, beaver2$Sex == "F")
data <- beavers_F  %>%
  group_by(Month, Forage_Sub_Cat) %>%
  summarise(n = sum(Tree+Aquatic+Grass_herb)) %>%
  mutate(percentage = n / sum(n)) %>% as_tibble()
data <- add_row(data, Month = 3, Forage_Sub_Cat = "Grass_herb", n = 0,percentage = 0)
data$Forage_Sub_Cat <- factor(data$Forage_Sub_Cat,
                              levels = c("Tree","Grass_herb","Aquatic"), 
                              labels = c("Trees/Shrubs", "Herbs/Grasses", "Aquatic Vegetation"))

p_F_percent <- ggplot(data, aes(x=Month, y=percentage, color=Forage_Sub_Cat)) + 
  geom_line(size = 2, alpha = 0.9) + 
  geom_point(size = 5) + 
  ggtitle(NULL) + 
  ylim(0,1) +
  theme(text = element_text(size = GGPLOTTEXTSIZE),
        legend.text=element_text(size=GGPLOTLEGENDSIZE),
        plot.subtitle=element_text(size=GGPLOTSTSIZE),
        legend.spacing.y = unit(0.3, 'cm') )+
  guides(fill = guide_legend(byrow = TRUE))+
  scale_color_manual(name = NULL, values = cols)+
  xlab("Month") + ylab("Percentage")
p_F_percent


ggsave(p_F_percent, filename = "figures/beaver_summary/beaver_summary_F_percent.pdf", width=PLOTWIDTH, height=PLOTHEIGHT)
ggsave(p_M_percent, filename = "figures/beaver_summary/beaver_summary_M_percent.pdf", width=PLOTWIDTH, height=PLOTHEIGHT)


