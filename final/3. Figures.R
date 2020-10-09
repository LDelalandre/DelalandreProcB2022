source("R/Common variables.R")
source("R/Before simulations.R")
library(cowplot)
library("gridExtra")
library(ggsignif)

# plot the final ecosystem property in each condition ####
# more precisely, it is the averaged biomass on 10 evenly-spaced points in time during the last 1000 years.
PLOT <- list()
i <- 0
measure <- "productivity_tot"
for (site in ord_plots){
  i <- i+1
  # for (measure in MEASURE){
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
      theme(plot.title = element_text(size=24)) +
      scale_x_continuous(breaks = 5*c(1:6)) +
      theme(legend.title = element_blank())+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none" ) + # virer tous les titres
      theme(axis.text=element_text(size=20))
      # ggsave(paste0("figures/",measure,"_",site,".png"))
  # }
  PLOT[[i]] <- plot
}

plot_prod_cold <- plot_grid(PLOT[[1]],PLOT[[2]],PLOT[[3]], ncol = 5, nrow = 4)
save_plot("paper/dist sp loss_cold sites.png",plot_prod_cold,base_height = 15)

plot_prod_warm_wet <- plot_grid(PLOT[[9]],PLOT[[10]],PLOT[[11]], ncol = 5, nrow = 4)
save_plot("paper/dist sp loss_warm-wet sites.png",plot_prod_warm_wet,base_height = 15)

plot_prod_warm_dry <- plot_grid(PLOT[[4]],PLOT[[5]],PLOT[[6]],PLOT[[7]],PLOT[[8]], ncol = 5, nrow = 4)
save_plot("paper/dist sp loss_warm-dry sites.png",plot_prod_warm_dry,base_height = 15)

# # construct the legend ####
# plot <-
#   ggplot(result,aes(x=simul-1,y=decreasing,color="Removing distinct species first")) +
#   labs(x="Number of species removed",y=measure) +
#   # geom_line(size=1)+
#   geom_ribbon(aes(ymin=int_min, ymax=int_max),fill="grey60", alpha=0.5,colour="black") +
#   geom_line(aes(x=simul-1,y=increasing, colour="#8766D"),size=2) +
#   geom_line(aes(x=simul-1,y=decreasing, colour="#00BFC4"),size=2) +
#   # geom_line(aes(x=simul-1,y=mean, colour="black"),size=1) +
#   # theme(legend.position = "bottom") +
#   ggtitle(site) +
#   scale_x_continuous(breaks = 2*c(1:15)) +
#   scale_color_manual(labels = c("Distinct species lost first", "Common species lost first"), values = c("#F8766D", "#00BFC4")) +
#   labs(color = "Order of species loss :")   +
#   theme(legend.background = element_rect(fill="white",
#                                            size=1, linetype="solid", 
#                                            colour ="black"))+
#   theme(legend.text = element_text( size=20),legend.title = element_text( size=20))
# 
# 
# get_legend<-function(myggplot){
#   tmp <- ggplot_gtable(ggplot_build(myggplot))
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#   legend <- tmp$grobs[[leg]]
#   return(legend)
# }
# 
# legend <- get_legend(plot)
# # grid.arrange(PLOT[[1]], legend, ncol=2, widths=c(2.3, 2.3))

# plot_prod <- plot_grid(PLOT[[1]],PLOT[[2]],PLOT[[3]],PLOT[[4]],PLOT[[5]],PLOT[[6]],PLOT[[7]],PLOT[[8]],PLOT[[9]],PLOT[[10]],PLOT[[11]],legend, ncol = 4, nrow = 3)
# save_plot("paper/Figure 2.png",plot_prod,base_height = 15)

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


# Plot biomass removal order####
PLOT <- list()
i <- 0
for (site in ord_plots){ # ATTENTION RE FAIRE TOURNER COTTBUS
  i <- i+1
  PROD <- NULL
  order <- "abundance"
  for(number in c(1:30)){
    prod <- try(read.table(paste0("data/raw/output-cmd2_",site,"_",order,".txt/forceps.",site,".site_",number,"_productivityScene.txt")),silent=T)
    if (class(prod) != "try-error"){# sometimes, the files are empty, and it returns an error message
      colnames(prod)<-colnames_prod
      years_to_keep <- max(prod$date) - c(900,800,700,600,500,400,300,200,100,0)
      prod_to_keep <- subset(prod,date %in% years_to_keep)
      prod_to_keep$totProdBiomass_t_ha <- prod_to_keep$adultProdBiomass_t_ha + prod_to_keep$saplingBiomass_t_ha
      
      # specific productivities
      productivities <- aggregate(prod_to_keep$totProdBiomass_t_ha, list(prod_to_keep$speciesShortName), mean)
      colnames(productivities) <- c("species","mixture_t_ha")
      productivities$site <- site
      productivities$order <- order
      productivities$simul <- number
      productivities$mixture_relative <- productivities$mixture_t_ha/sum(productivities$mixture_t_ha)
      PROD <- rbind(PROD,productivities)
    }
    
  }
  PROD_sp <- PROD
  
  
  
  PROD_tot <- aggregate(PROD_sp$mixture_t_ha,list(site = PROD_sp$site,order = PROD_sp$order,simul = PROD_sp$simul),FUN=sum)# sum of the productivities of the species
  
  prod_per_order <- spread(PROD_tot,order,x) # data frame with each order in column
  
  measure <- "productivity_tot"
  result <- read.table(paste0("data/processed/",measure,"_",site,"_with interval_median.txt"),header=T)
  abundance <- prod_per_order$abundance
  result <- cbind(result,abundance)
  
  plot <- ggplot(result,aes(x=simul-1,y=abundance,color="green")) +
    scale_color_manual(values=c("green4")) +
    labs(x="Number of species removed",y=measure) +
    geom_line()+
    geom_ribbon(aes(ymin=int_min, ymax=int_max),fill="grey60", alpha=0.5,colour="black") +
    theme(legend.position = "bottom") +
    ggtitle(site) +
    theme(plot.title = element_text(size=18)) +
    scale_x_continuous(breaks = 5*c(1:6)) +
    theme(legend.title = element_blank())+
    theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none" ) + # virer tous les titres
    theme(axis.text=element_text(size=15))
  
  PLOT[[i]] <- plot
}
plot_ab <- plot_grid(PLOT[[1]],PLOT[[2]],PLOT[[3]],PLOT[[4]],PLOT[[5]],PLOT[[6]],PLOT[[7]],PLOT[[8]],PLOT[[9]],PLOT[[10]],PLOT[[11]], ncol = 4, nrow = 3)
save_plot("paper/Figure abundance.png",plot_ab,base_height = 10)



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
    theme(axis.title.x=element_blank(),legend.position = "none" ) +
    geom_signif(comparisons = list(c("Common", "Distinct")),  map_signif_level=F,textsize=10,color="black") +
    # ylim(0,500) +
    theme(axis.title.y=element_blank(),axis.text.x=element_blank()) +
    scale_y_sqrt(limits=c(0,600)) +
    theme(axis.text=element_text(size=25),plot.title = element_text(size=28))
  
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

plot_mono_cold <- plot_grid(PLOT[[1]],PLOT[[2]],PLOT[[3]], ncol = 5, nrow = 4)
save_plot("paper/Boxplot mono_cold sites.png",plot_mono_cold,base_height = 15)

plot_mono_warm_wet <- plot_grid(PLOT[[9]],PLOT[[10]],PLOT[[11]], ncol = 5, nrow = 4)
save_plot("paper/Boxplot mono_warm-wet sites.png",plot_mono_warm_wet,base_height = 15)

plot_mono_warm_dry <- plot_grid(PLOT[[4]],PLOT[[5]],PLOT[[6]],PLOT[[7]],PLOT[[8]], ncol = 5, nrow = 4)
save_plot("paper/Boxplot mono_warm-dry sites.png",plot_mono_warm_dry,base_height = 15)

# plot_mono <-  plot_grid(PLOT[[1]],PLOT[[2]],PLOT[[3]],PLOT[[4]],PLOT[[5]],PLOT[[6]],PLOT[[7]],PLOT[[8]],PLOT[[9]],PLOT[[10]],PLOT[[11]], ncol = 4, nrow = 3)
# save_plot("paper/Figure 3_a.png",plot_mono,base_height = 10)

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
    theme(axis.title.x=element_blank(),legend.position = "none" ) +
    geom_signif(comparisons = list(c("Common", "Distinct")),  map_signif_level=F,textsize=10,color="black") +
    # ylim(0,max(biom$mixture.t.ha.)+50) +
    scale_y_sqrt(limits=c(0,350)) +
    theme(axis.title.y=element_blank(),axis.text.x=element_blank()) +
    theme(axis.text=element_text(size=25),plot.title = element_text(size=28))
    # stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=8,show_guide = FALSE)
  
  # test <- wilcox.test(biom$mixture.t.ha.~biom$med)
  # pval <- c(pval,test$p.value)
  # NB: In wilcox.test.default(x = c(0, 0, 8.4, 10.1, 0, 18, 1.9, 0, 2429.9,  :
  # impossible de calculer la p-value exacte avec des ex-aequos. C'est grave ?
  
  PLOT[[i]] <- plot
}

plot_mixt_cold <- plot_grid(PLOT[[1]],PLOT[[2]],PLOT[[3]], ncol = 5, nrow = 4)
save_plot("paper/Boxplot mixt_cold sites.png",plot_mixt_cold,base_height = 15)

plot_mixt_warm_wet <- plot_grid(PLOT[[9]],PLOT[[10]],PLOT[[11]], ncol = 5, nrow = 4)
save_plot("paper/Boxplot mixt_warm-wet sites.png",plot_mixt_warm_wet,base_height = 15)

plot_mixt_warm_dry <- plot_grid(PLOT[[4]],PLOT[[5]],PLOT[[6]],PLOT[[7]],PLOT[[8]], ncol = 5, nrow = 4)
save_plot("paper/Boxplot mixt_warm-dry sites.png",plot_mixt_warm_dry,base_height = 15)

# plot_mixt <- plot_grid(PLOT[[1]],PLOT[[2]],PLOT[[3]],PLOT[[4]],PLOT[[5]],PLOT[[6]],PLOT[[7]],PLOT[[8]],PLOT[[9]],PLOT[[10]],PLOT[[11]], ncol = 4, nrow = 3)
# save_plot("paper/Figure 3_b.png",plot_mixt,base_height = 10)

library(ggsignif)
plot + stat_summary(fun.y=mean, colour="black", geom="point", 
                  shape=18, size=8,show_guide = FALSE)

# + geom_signif(comparisons = list(c("Distinct", "Common")), 
            # map_signif_level=TRUE) #,aes(color="black"))


library(tidyverse)
library(ggsignif)


# catégorisation des sites ####
sites <- read.table("data/Site description.txt",header=T)
# site_descr <- 
  
  ggplot(data=sites,aes(x=Temp_moy,y=Annual_ppt,label=site)) +
  geom_point() +
  labs(x="Mean temperature (°C)",y="Annual precipitation (mm)") +
  ggrepel::geom_text_repel(aes(label = Site),size = 8)+
  theme(axis.title=element_text(size=21),axis.text=element_text(size=15))

ggsave("paper/figure 1_a.png",plot=site_descr,dpi="print",width=20,units="cm")

# Figure 1.B: PCA on the traits ####
traits<-read.table("data/Traits of the species_complete.txt",header=T)
c1<-choice_traits_1(traits) # data.frame with the traits of the species
# cprime <- select(c1,-c(A1max,A2))
ACP1<-PCA(c1)
c1$pc1 <- ACP1$ind$coord[, 1] 
c1$pc2 <- ACP1$ind$coord[, 2] 
c1$SName <- rownames(c1)
distinct_tot <- read.table("data/raw/distinctiveness of the species.txt",header = T)
c1$Distinctiveness <- distinct_tot$Di

plot_pca <- ggplot(data=c1,aes(x=pc1,y=pc2)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = SName),size = 5) +
  aes(color=Distinctiveness) +
  scale_colour_gradient(low="#00BFC4",high="#F8766D") +
  # scale_colour_gradient(low = "#132B43",high = "#56B1F7" ) +
  labs(x="Dim 1 (28.8%)",y="Dim2 (22.9%)") +
  theme(axis.title=element_text(size=15),axis.text=element_text(size=12))
ggsave ("paper/figure 1_b.png",plot=plot_pca,dpi="print",width=20,units="cm")

# table 1: traits ####
traits<-read.table("data/Traits of the species_complete.txt",header=T)
tr<-choice_traits_1(traits) # data.frame with the traits of the species
snsn <- read.table("data/correspondence_SName_Id.txt",header=T)
Species_name <- snsn$Name
Short_name <- snsn$SName


# cvcv <- read.table("data/distinctiveness with or without A1 and A2.txt",header=T)
Distinctiveness <- distinct_tot$Di

table1 <- cbind(Species_name,Short_name,tr,Distinctiveness)
write.table(table1,"paper/table1.txt",row.names=F,sep="\t")

