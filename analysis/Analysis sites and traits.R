source("R/Common variables.R")
source("R/Analysis_data.R")
source("R/Monocultures_functions.R")
source("R/Before simulations.R")

TOTAL <- read.table("data/processed/specific_biom_prod_complete.txt",header=T)
# Traits et sites
traits <- read.table("data/traits of the species_complete.txt",header=T)
sites <- read.table("data/Site description.txt",header=T)
rownames(sites) <- sites$Site

orders <- read.table("data/removal order and sp names.txt",header=T)
# Separate the sites depending on when the highest drop of ecosystem function occurs when removing species. ####
measure <- MEASURE[2]

DEC <- data.frame(matrix(0,nrow=11,ncol=30))
rownames(DEC) <- SITE
INC <- DEC
colnames(DEC) <- orders$species_decr # gives the species that is removed when that drop occurs
colnames(INC) <- orders$species_incr

i=0
for (site in SITE){
  i <- i+1
  B <- read.table(here::here(paste0("data/processed/",measure,"_",site,".txt")),header=T)
  #decreasing
  diff <- B$decreasing - c(B$decreasing[2:30],0)
  DEC[i,] <- diff # decrease in biomass
  
  #increasing
  diff <- B$increasing - c(B$increasing[2:30],0)
  INC[i,] <- diff
}



drop_dec <- c()
drop_inc <- c()
for (i in 1:length(SITE)){
  dec_highest_drop <- which(DEC[i,]==max(DEC[i,]))
  drop_dec <- c(drop_dec,dec_highest_drop)
  
  inc_highest_drop <- which(INC[i,]==max(INC[i,]))
  drop_inc <- c(drop_inc,inc_highest_drop)
}
data.frame(SITE,drop_dec,drop_inc)
# orders[which(orders$nb_sp_lost == dec_highest_drop),]$species_decr

# To have qualitatively the quantiles when biomass increases or decreases
# chose the proportion you want to keep:
# loss <- 0.9 # keep the 10% species who generate the higest loss of biomass when removed
# gain <- 0.1 # keep the 10% species who generate the highest gain of biomass when removed

# i=0
# for (site in SITE){
#   i <- i+1
#   B <- read.table(paste0("data/processed/biomass_tot_",site,".txt"),header=T)
#   #decreasing
#   diff <- B$decreasing - c(B$decreasing[2:30],0)
#   DEC[i,which(diff>quantile(diff,loss))] <- -1 # decrease in biomass
#   DEC[i,which(diff<quantile(diff,gain) & diff<0)] <- 1 # increase in biomass
#   
#   #increasing
#   diff <- B$increasing - c(B$increasing[2:30],0)
#   INC[i,which(diff>quantile(diff,loss))] <- -1
#   INC[i,which(diff<quantile(diff,gain) & diff<0)] <- 1 # increase in biomass

# Look at that drop for the "aggregated monocultures" ####
measure <- MEASURE[2]

DEC <- data.frame(matrix(0,nrow=11,ncol=30))
rownames(DEC) <- SITE
INC <- DEC
colnames(DEC) <- orders$species_decr
colnames(INC) <- orders$species_incr
DEC_mixt <- DEC
INC_mixt <- INC

i=0
for (site in SITE){
  i <- i+1
  #decreasing
  if (measure == "productivity_tot"){
    B <- read.table(paste0("data/processed/prod_removal_mixt_mono_",site,"_decreasing.txt"),header=T)
  } else{
    B <- read.table(paste0("data/processed/removal_mixt_mono_",site,"_decreasing.txt"),header=T)
  }
  
  diff <- B$Monoculture_biomass - c(B$Monoculture_biomass[2:30],0)
  diff_mixt <- B$Mixture_biomass - c(B$Mixture_biomass[2:30],0)
  DEC[i,] <- diff # decrease in biomass
  DEC_mixt[i,] <- diff_mixt
  
  #increasing
  if (measure == "productivity_tot"){
    B <- read.table(paste0("data/processed/prod_removal_mixt_mono_",site,"_increasing.txt"),header=T)
  } else{
    B <- read.table(paste0("data/processed/removal_mixt_mono_",site,"_increasing.txt"),header=T)
  }
  
  diff <- B$Monoculture_biomass - c(B$Monoculture_biomass[2:30],0)
  diff_mixt <- B$Mixture_biomass - c(B$Mixture_biomass[2:30],0)
  INC[i,] <- diff # decrease in biomass
  INC_mixt[i,] <- diff_mixt
}

drop_dec <- c()
drop_dec_mixt <- c()
drop_inc <- c()
drop_inc_mixt <- c()
for (i in 1:length(SITE)){
  dec_highest_drop <- which(DEC[i,]==max(DEC[i,]))
  drop_dec <- c(drop_dec,dec_highest_drop)
  drop_dec_mixt <- c(drop_dec_mixt, which(DEC_mixt[i,]==max(DEC_mixt[i,])) )
  
  inc_highest_drop <- which(INC[i,]==max(INC[i,]))
  drop_inc <- c(drop_inc,inc_highest_drop)
  drop_inc_mixt <- c(drop_inc_mixt, which(INC_mixt[i,]==max(INC_mixt[i,])) )
}
DROP <- data.frame(SITE,drop_dec,drop_dec_mixt,drop_inc,drop_inc_mixt)

# Plot that drop
ggplot(DROP,aes(x=SITE,y=drop_dec_mixt)) +
  labs(x="Site",y="Highest drop in productivity") +
  geom_boxplot() + 
  # geom_boxplot(aes(x=SITE,y=drop_dec, color="Monocultures")) +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=18, face="bold"),
        axis.title.y = element_text(size=20, face="bold")
        ) +
  ggsave("figures/Drop_productivity_mixtures.png")


# Computes "functional distinctiveness" on the SITES ####
# first, add drought index
drought <- read.table("data/drought index.txt",header=T)
drought <- arrange(drought,factor(drought$SITE,levels = sites$Site ))
site_to_keep <- select(sites,c(Temp_moy, Temp_min, Annual_ppt, Na.kg.ha.an., ProdMax.T.ha.an.))
site_val <- data.frame(site_to_keep,drought[,2])
distinctiveness <- traits_dist(site_val)
distinctiveness[order(distinctiveness$Di),]

ACP1 <- PCA(site_val) # site_val are quantitative values
ACP2 <- PCA(select(site_val,c(Temp_moy,Temp_min, Annual_ppt, Na.kg.ha.an.,drought...2.)))
ACP <- ACP2 # choose here which of the PCA you want to kepp for the following analysis

# variance explained by the axis
fviz_eig(ACP, addlabels = TRUE, ylim = c(0, 30)) # I keep the first two domensions only

# position of the site on the two first axis
distACP <- ACP1$ind$coord %>% 
  as.data.frame() %>%
  select(Dim.1    ,   Dim.2     ) %>%
  traits_dist() # here we compute the functional distinctiveness

fviz_pca_biplot(ACP, repel = TRUE, # biplot
                col.var = "#2E9FDF", # Couleur des variables
                col.ind = "#696969"  # Couleur des individus 
)

site_dist <-  distACP
site_dist[order(site_dist$Di),]

site_dist$mean <- (site_dist$Dim.1+site_dist$Dim.2)/2 #average the two first dimensions of the PCA
site_dist[order(site_dist$mean),]

# # Try to average all the variables. But I should scale them before ?
# MOY <- mutate(site_val,mean=(Temp_moy+ Annual_ppt+ Na.kg.ha.an.+ ProdMax.T.ha.an.+ drought...2.)/5)
# rownames(MOY) <- rownames(site_val)
# MOY[order(MOY$mean),]

# Other option : rank the sites according to their maximal productivity only
site_val[order(site_val$ProdMax.T.ha.an.),]
# ... it almost isolates the sites dependent on distinct species (but Cottbus is found within this group)
site_val[order(site_val$Annual_ppt),]

# Predict functioning from trait values ####
# Select the traits we use for the analysis (the ones used for computing the distinctiveness)
traits2 <- choice_traits_1(traits)
traits2 <- select(traits2,-c(A1max,A2))# remove traits that say explicitely that a species is a gymnosperm

# Add columns  with the biomas and productivity of the species to predict them from the traits - for each site
sit <- SITE[1]
#  [1] "GrandeDixence" "Bever"         "Davos"         "Adelboden"     "Huttwil"       "Schwerin"      "Bern"         
# [8] "Cottbus"       "Basel"         "Schaffhausen"  "Sion"    
orde <- ORDER[1]
simu <- 1
SUB <- subset(TOTAL, site==sit & order==orde & simul ==simu)
SUB2 <- arrange(SUB,factor(SUB$species,levels = traits$SName )) # species in the same order as in the trait data frame
prod_mixture <- SUB2$prod_mixture
pooled <- data.frame(prod_mixture,traits2)
# pooled<-subset(pooled,!(prod_mixture<0)) # pour pouvoir transformer les variables, il me faut une variable réponse positive
 
# library(corrplot)
# corrplot.mixed(cor(pooled), order="hclust", tl.col="black")

# linear model (12 potential explanatory variables !!)
null <- lm(prod_mixture^0.1~1,data=pooled)
full <- lm(prod_mixture^0.1~ S +   HMax + AMax+  G  +  DDMin  +   WiTN +WiTX + DrTol    +        
                       NTol      +  Brow + Ly  ,data=pooled) # puissance 0.25 : pas mal !
par(mfrow=c(2,2))
plot(full)

step(null, scope = ~ S +   HMax + AMax+  G  +  DDMin  +   WiTN +WiTX + DrTol    +        
       NTol      +  Brow + Ly ,
     direction="both", criterion = "BIC")
step(full, scope = ~ 1 ,
     direction="both", criterion = "AIC")



library(car)
predictors <- read.table("data/processed/Prediction productivity~traits.txt",header=T)
pred <- t(predictors[,-1])
colnames(pred) <- predictors$SITE
par(mfrow=c(1,1))
corrplot.mixed(cor(pred,method="spearman"), order="hclust", tl.col="black")



par(mfrow=c(2,2)) ; plot(model.complete)
# model.nul=lm(LogNbNids~1,data=Chenilles.cr)
# select.variables=step(model.nul,scope=~Altitude+Pente+NbPins+Hauteur+Diametre+...
#                       ... Densite+Orient+HautMax+NbStrat+Melange,direction="both",data=Chenilles.cr)
# summary(select.variables)

anova(lm(mixture.t.ha. ~ WiTX,data=pooled))
anova(model.complete)     

pooled
cor(pooled)

# Correlate biomass and distinctiveness in monoculture ####
site <- SITE[1]
correlations <- c()
par(mfrow=c(3,4))
for (site in SITE){
  biom <- read.table(paste0("data/processed/biomass_monoculture_",site,".txt"),header=T)
  correlations <- c(correlations,cor(biom$Di,biom$monoculture.t.ha.) )
  plot(biom$Di,biom$monoculture.t.ha.,main=site)
}
cor <- data.frame(SITE,correlations)
par(mfrow=c(1,1))
plot(cor$cor~cor$SITE)        

# Look at the distinctiveness of the remaining species ####
sit <- "Bever"
current_site <- read.table(paste0("data/processed/biomass_monoculture_",sit,".txt"),header=T)
sp <- as.character(subset(bever,monoculture.t.ha.>0)$SName)
traits2 <- subset(traits,SName %in% sp)

# compute functional distinctiveness
comp_fct_dist <- function(traits){
  c1<-choice_traits_1(traits2) # data.frame with the traits of the species
  ACP1<-PCA(c1)
  distACP <- ACP1$ind$coord %>% 
    as.data.frame() %>%
    select(Dim.1    ,   Dim.2     ,   Dim.3      ,  Dim.4) %>%
    traits_dist() # here we compute the functional distinctiveness
  distACP$SName<-rownames(distACP)
  
  distinct_tot <-  distACP
  distinct_tot[order(distinct_tot$Di),]
}

comp_fct_dist(traits2)

# On dirait que, quelque soit le site, les espèces les plus distinctes 
#  dans le pool local sont celles qui sont le pool régional
