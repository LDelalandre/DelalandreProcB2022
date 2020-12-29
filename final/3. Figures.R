source("R/Common variables.R")
source("R/Before simulations.R")
library(cowplot)
library("gridExtra")
library(ggsignif)

distinct_tot <- read.table("data/raw/distinctiveness of the species.txt",header = T)



# Fig 2: PCA on the traits ####
traits<-read.table("data/Traits of the species_complete.txt",header=T)
c1<-choice_traits_1(traits) # data.frame with the traits of the species
# cprime <- select(c1,-c(A1max,A2))
ACP1<-PCA(c1)
c1$pc1 <- ACP1$ind$coord[, 1] 
c1$pc2 <- ACP1$ind$coord[, 2] 
c1$SName <- rownames(c1)
c1$Distinctiveness <- distinct_tot$Di

plot_pca <- ggplot(data=c1,aes(x=pc1,y=pc2)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = SName),size = 5) +
  aes(color=Distinctiveness) +
  scale_colour_gradient(low="#00BFC4",high="#F8766D") +
  # scale_colour_gradient(low = "#132B43",high = "#56B1F7" ) +
  labs(x="Dim 1 (28.8%)",y="Dim2 (22.9%)") +
  theme(axis.title=element_text(size=15),axis.text=element_text(size=12))
ggsave ("paper2/PCA.png",plot=plot_pca,dpi="print",height=15,width=20,units="cm")


# Table 1: traits ####
traits<-read.table("data/Traits of the species_complete.txt",header=T)
tr<-choice_traits_1(traits) # data.frame with the traits of the species
snsn <- read.table("data/correspondence_SName_Id.txt",header=T)
Species_name <- snsn$Name
Short_name <- snsn$SName

Distinctiveness <- distinct_tot$Di
table1 <- cbind(Species_name,Short_name,tr,Distinctiveness)
write.table(table1,"paper2/table1.txt",row.names=F,sep="\t")


# Correlation biomass-distinctiveness in monoculture ####
MONO <- read.table("data/processed/biomass_mono_ALL sites.txt",header=T)

MONOCULTURES <- NULL
for (sit in SITE){
  MONOsite <- 
    MONO %>% 
    filter(site==sit) %>% 
    mutate(site=sit) 
    # filter(persists==T) # If I keep only the local pool (species above biomass threshold), I lose too much statistical power
  # consequently, I don't filter them 
  DIST <- 
    read.table("data/raw/distinctiveness of the species.txt",header=T) %>% 
    filter(SName %in% MONOsite$SName) %>% 
    arrange(factor(SName,levels=MONOsite$SName))
  
  MONOsite$Di <- DIST$Di
  MONOCULTURES <- rbind(MONOCULTURES,MONOsite)
}


plotcor <- ggplot(MONOCULTURES,aes(x=Di,y=monoculture.t.ha.,label=SName))+
  geom_point()+
  facet_wrap(~site)+
  geom_smooth(method=lm)+
  ggpubr::stat_cor(method="spearman")+
  geom_label()

ggsave(filename = paste0("paper_2/correlation biomass_Di.png"), 
       plot = plotcor,
       width = 50, 
       height = 40,
       units = "cm",
       dpi = 300)

