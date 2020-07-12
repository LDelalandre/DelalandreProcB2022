source("R/Common variables.R")
library(cowplot)
library("gridExtra")

# plot the final ecosystem property in each condition ####
# more precisely, it is the averaged biomass on 10 evenly-spaced points in time during the last 1000 years.
PLOT <- list()
i <- 0
for (site in SITE){
  i <- i+1
  for (measure in MEASURE){
    result <- read.table(paste0("data/processed/",measure,"_",site,"_with interval_median.txt"),header=T)
    plot <-
      ggplot(result,aes(x=simul-1,y=decreasing,color="Removing distinct species first")) +
      labs(x="Number of species removed",y=measure) +
      # geom_line(size=1)+
      geom_ribbon(aes(ymin=int_min, ymax=int_max),fill="grey60", alpha=0.5,colour="black") +
      geom_line(aes(x=simul-1,y=increasing, color="#8766D"),size=2) +
      geom_line(aes(x=simul-1,y=decreasing, color="#00BFC4"),size=2) +
      theme(legend.position = "bottom") +
      ggtitle(site) +
      scale_x_continuous(breaks = 2*c(1:15)) +
      theme(legend.title = element_blank())+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none" ) # virer tous les titres
      # ggsave(paste0("figures/",measure,"_",site,".png"))
  }
  PLOT[[i]] <- plot
}


# construct the legend ####
plot <-
  ggplot(result,aes(x=simul-1,y=decreasing,color="Removing distinct species first")) +
  labs(x="Number of species removed",y=measure) +
  # geom_line(size=1)+
  geom_ribbon(aes(ymin=int_min, ymax=int_max),fill="grey60", alpha=0.5,colour="black") +
  geom_line(aes(x=simul-1,y=increasing, colour="#8766D"),size=2) +
  geom_line(aes(x=simul-1,y=decreasing, colour="#00BFC4"),size=2) +
  # geom_line(aes(x=simul-1,y=mean, colour="black"),size=1) +
  # theme(legend.position = "bottom") +
  ggtitle(site) +
  scale_x_continuous(breaks = 2*c(1:15)) +
  scale_color_manual(labels = c("Distinct species lost first", "Common species lost first"), values = c("#F8766D", "#00BFC4")) +
  labs(color = "Order of species loss :")   +
  theme(legend.background = element_rect(fill="white",
                                           size=1, linetype="solid", 
                                           colour ="black"))


get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(plot)
grid.arrange(PLOT[[1]], legend, ncol=2, widths=c(2.3, 2.3))

plot_grid(PLOT[[1]],PLOT[[2]],PLOT[[3]],PLOT[[4]],PLOT[[5]],PLOT[[6]],PLOT[[7]],PLOT[[8]],PLOT[[9]],PLOT[[10]],PLOT[[11]],legend, ncol = 4, nrow = 3)


# areas = f(species)
sit <- "Bever"
orde <- "decreasing"
data <- read.table("data/processed/specific_biom_prod_complete.txt",header=T)
sub <- subset(data,site==sit&order==orde)
ggplot(sub, aes(x=simul, y=mixture.t.ha., fill = species)) +
  geom_area()

sub$species <- factor(sub$species,levels=)
ggplot(sub, aes(simul, mixture.t.ha., group = species)) +
  geom_area(aes(fill = species),position = "stack") #+
  geom_text(aes(label = species), position = position_stack(vjust = 0.5))
  
ggplot(sub, aes(simul, mixture.t.ha.)) +
  geom_area(aes(fill = species))



# Boxplots mono two groups of distinctiveness ####
# identify the most distinct and most common halves with SName (for mixtures, where I don't have distinctiveness in the data frame)
PLOT <- list()
pval <- c()
i <- 0
for(sit in ord_plots){ # ordonner pour regrouper les sites (warm-wet sites...)
  i <- i+1
  biom <- read.table(paste0("data/processed/biomass_monoculture_",sit,".txt"),header=T)
  biom <- biom %>% 
    mutate(med = as.factor(ifelse(Di > median(biom$Di), "Distinct","Common"))) %>% 
    mutate(color = ifelse(med=="Common","#8766D","#00BFC4") )
  # shapiro.test : pas de normalité dans ces classes
  
  plot <- ggplot(biom, aes(x=med, y=monoculture.t.ha.,fill=color , colour=color)) +
    geom_boxplot() +
    ggtitle(sit) +
    labs(y="Biomass (t/ha)") +
    theme(axis.title.x=element_blank(),legend.position = "none" )
  
  # plot <- biom %>%
  #   ggplot(aes(abundance, fill=color, colour=color)) +
  #   geom_density(alpha=0.25)+
  #   geom_rug()
  test <- wilcox.test(biom$monoculture.t.ha.~biom$med)
  pval <- c(pval,test$p.value)
  # NB: In wilcox.test.default(x = c(0, 0, 8.4, 10.1, 0, 18, 1.9, 0, 2429.9,  :
  # impossible de calculer la p-value exacte avec des ex-aequos. C'est grave ?
  
  PLOT[[i]] <- plot
}


plot_grid(PLOT[[1]],PLOT[[2]],PLOT[[3]],PLOT[[4]],PLOT[[5]],PLOT[[6]],PLOT[[7]],PLOT[[8]],PLOT[[9]],PLOT[[10]],PLOT[[11]], ncol = 4, nrow = 3)


# Boxplots mixtures two groups of distinctiveness ####
dist <- read.table("data/raw/distinctiveness of the species.txt",header=T) %>% 
  mutate(med = as.factor(ifelse(Di > median(dist$Di), "Distinct","Common"))) %>% 
  mutate(color = ifelse(med=="Common","#8766D","#00BFC4") )
distinct_half <- subset(dist,med=="Distinct")

PLOT <- list()
pval <- c()
i <- 0
for(sit in ord_plots){ # ordonner pour regrouper les sites (warm-wet sites...)
  i <- i+1
  biom <- TOTAL
  biom <- biom %>% 
    subset(simul==1 & order=="increasing"&site==sit) %>% 
    mutate(med = as.factor(ifelse(species %in% distinct_half$SName, "Distinct","Common"))) %>%  
    mutate(color = ifelse(med=="Common","#8766D","#00BFC4") )
  # shapiro.test : pas de normalité dans ces classes
  
  plot <- ggplot(biom, aes(x=med, y=mixture.t.ha.,fill=color , colour=color)) +
    geom_boxplot() +
    ggtitle(sit) +
    labs(y="Biomass (t/ha)") +
    theme(axis.title.x=element_blank(),legend.position = "none" )
  
  test <- wilcox.test(biom$mixture.t.ha.~biom$med)
  pval <- c(pval,test$p.value)
  # NB: In wilcox.test.default(x = c(0, 0, 8.4, 10.1, 0, 18, 1.9, 0, 2429.9,  :
  # impossible de calculer la p-value exacte avec des ex-aequos. C'est grave ?
  
  PLOT[[i]] <- plot
}


plot_grid(PLOT[[1]],PLOT[[2]],PLOT[[3]],PLOT[[4]],PLOT[[5]],PLOT[[6]],PLOT[[7]],PLOT[[8]],PLOT[[9]],PLOT[[10]],PLOT[[11]], ncol = 4, nrow = 3)
