source("final/0. Packages.R")
source("R/Common variables.R")
source("R/Before simulations.R")
source("R/comp_fct_dist.R")


distinct_tot <- read.table("data/raw/distinctiveness of the species.txt",header = T)


#________________________________________________________________________________
# 1) Fig. 3: PCA and bootstrap ####
# i) PCA ####
traits<-read.table("data/raw/Traits of the species_complete.txt",header=T)
c1<-select_traits(traits) # data.frame with the traits of the species
ACP1<-PCA(c1,graph=F)
c1$Dim.1 <- ACP1$ind$coord[, 1] 
c1$Dim.2 <- ACP1$ind$coord[, 2] 
c1$SName <- rownames(c1)
c1$Name <- traits$Name
c1$Distinctiveness <- distinct_tot$Di
c1$Distinctiveness <- c1$Distinctiveness/max(c1$Distinctiveness) # to have a relative distinctiveeness (good idea?)

axis <- ACP1$var$coord[,c(1,2)] %>% 
  data.frame() %>% 
  rownames_to_column(var="varnames")

axis2 <- transform(axis,
          Dim.1 = 4.7 * Dim.1,
          Dim.2 = 4.7* Dim.2)

var.explain.dim1 <- round(ACP1$eig[1,2],digits=1)
var.explain.dim2 <- round(ACP1$eig[2,2],digits=1)

plot_pca <-
  ggplot(data=c1,aes(x=Dim.1,y=Dim.2)) + 
  geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
  geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
  coord_equal() + 
  geom_text(data=axis2, aes(x=Dim.1, Dim.2, label=varnames), size = 5, vjust=1, color="black")+
  geom_segment(data=axis2, aes(x=0, y=0, xend=Dim.1-0.2, yend=Dim.2-0.2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  ggrepel::geom_label_repel(aes(label = Name),size = 5) + # or SName
  aes(color=Distinctiveness) +
  scale_colour_gradient(low="#00BFC4",high="#F8766D") +
  # scale_colour_gradient(low = "#132B43",high = "#56B1F7" ) +
  labs(x=paste0("Dim 1 (",var.explain.dim1,"%)"),
       y=paste0("Dim 2 (",var.explain.dim2,"%)") ) +
    theme_classic() +
  theme(axis.title=element_text(size=15),axis.text=element_text(size=15))
  

ggsave("figures_tables/PCA.png",plot_pca,height=26,width=33,units="cm",dpi="print")

# ii) SENSITIVITY ANALYSIS: bootstrap ####
traits<-read.table("data/raw/Traits of the species_complete.txt",header=T)
selected_traits<-select_traits(traits) # data.frame with the traits of the species

# Initiate distinctiveness computation
init_dist <- 
  comp_fct_dist(selected_traits) %>% 
  pull(Di)

# Perform bootstrap (subsample traits with replacement, and correlate species'
# distinctiveness ranking whith Spearman's rho).
corstat <- c()
for (i in c(1:10000)){
  col <- sample(c(1:14),size=14,replace=T) # sample with replacement the traits
  boot_traits <- selected_traits[,col]
  boot_dist <- comp_fct_dist(boot_traits) %>% 
    pull(Di)
  boot_cor <- cor(init_dist,boot_dist,method="spearman")
  corstat <- c(corstat,boot_cor)
}
write.table(corstat,"figures_tables/bootstrap on Di correlations.txt",row.names=F,col.names=F)

# Histogram of Spearman's rho for the bootstraps
hist_bootstrap <- ggplot(data.frame(corstat),aes(x=corstat)) +
  geom_histogram(col="black",fill="grey",binwidth = 0.05) +
  labs(x="Rho",y="frequency") +
  theme_classic() +
  theme(axis.title=element_text(size=15),axis.text=element_text(size=15))

# iii) Combine the two plots ####
fig2 <- plot_pca +
  annotation_custom(ggplotGrob(hist_bootstrap), 
                    xmin = 2.6, xmax = 5.4, ymin = 1.5, ymax = 3.8) + 
  geom_text(  label="A", x=-4,y=3.9,size = 13,color = "black") +
  geom_text(  label="B", x=3.9,y=3.9,size = 13,color = "black")

ggsave("figures_tables/Fig.3_pca.png",fig2,height=26,width=29,units="cm",dpi="print")


# 2) Fig. Sx: Variance Explained ####
percent_var <- factoextra::fviz_eig(ACP1, addlabels = TRUE, ylim = c(0, 30))
ggsave(filename = "figures/Fig.Sx_variance.png",plot = percent_var)


# 3) Fig.Sx: Richness ####
LH_ALL <- read.csv("data/processed/Loreau-Hector coefficients_per species.csv")

per_mono <- LH_ALL %>% 
  filter(simul==1&order=="random_3") %>% # With 30 species
  group_by(site) %>% 
  filter(persists_mono==T) %>%
  summarize(RS_mono=n()) %>% 
  arrange(factor(site, levels = SITE))

per_mixt <- LH_ALL %>% 
  filter(simul==1&order=="random_3") %>% # With 30 species
  group_by(site) %>% 
  filter(persists_mixt==T) %>%
  summarize(RS_mixt=n()) %>% 
  arrange(factor(site, levels = SITE)) 

richness <- merge(per_mono,per_mixt,by="site") %>% 
  arrange(factor(site, levels = SITE)) %>% 
  mutate(site=factor(site, levels=site))
ggplot(richness,aes(x=site,y=RS_mono))+
  geom_histogram(stat="identity",alpha=0.2) +
  theme(axis.text.x = element_text(angle = 60)) +
  geom_point( aes(x=site,y=RS_mixt),stat="identity")+
  ylab("Richness") +
  ggsave("figures_tables/Richness.png")


# 4) Fig.Sx: temporal change in biomass ####
# mean biomass in time

MEAN <- NULL
for (site in SITE){
  mean <- read.table(paste0("data/raw/Output_ForCEEPS/",site,"/output-cmd2_",site,"_decreasing.txt/forceps.",site,".site_1_mean.txt"))
  colnames(mean) <- colnames_mean
  mean2 <- mean %>%
    select(date,nTrees..ha.,totalBiomass.t.ha.) %>% 
    mutate(date = date-1950)
  mean2$site <- site
  MEAN <- rbind(MEAN,mean2)
}

MEAN2 <- MEAN %>% 
  mutate(site=factor(site, levels = SITE))
ggplot(MEAN2,aes(x=date,y=totalBiomass.t.ha.)) +
  facet_wrap(~site)+
  geom_line() +
  xlab("Year") +
  ylab("Biomass (t/ha)") +
  ggsave("figures_tables/Temporal_evolution_biomass.png",height=9,width=11)
