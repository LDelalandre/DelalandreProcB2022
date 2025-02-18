source("R/Common variables.R")
source("R/Analysis_data.R")
source("R/Monocultures_functions.R")
source("R/Before simulations.R")
library(stringr)

TOTAL <- read.table("data/processed/specific_biom_prod_complete.txt",header=T)
# Traits et sites
traits <- read.table("data/traits of the species_complete.txt",header=T)
sites <- read.table("data/Site description.txt",header=T)
rownames(sites) <- sites$Site

orders <- read.table("data/removal order and sp names.txt",header=T)

# Correlation between traits ####
c1<-choice_traits_1(traits) # data.frame with the traits of the species

library("ggcorrplot")
corr <- cor(c1)
plot_cor <- ggcorrplot(corr,
           hc.order = TRUE,
           type = "lower",
           lab = TRUE,digits=1)
ggsave ("paper/correlation between traits.png",plot=plot_cor,dpi="print",width=17,units="cm")

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

#  [1] "GrandeDixence" "Bever"         "Davos"         "Adelboden"     "Huttwil"       "Schwerin"      "Bern"         
# [8] "Cottbus"       "Basel"         "Schaffhausen"  "Sion"    
orde <- ORDER[1]
simu <- 1

sit <- SITE[1]
SUB <- subset(TOTAL, site==sit & order==orde & simul ==simu)
SUB2 <- arrange(SUB,factor(SUB$species,levels = traits$SName )) # species in the same order as in the trait data frame
prod_mixture <- SUB2$prod_mixture
biomass <- SUB2$mixture.t.ha.
pooled <- data.frame(biomass,prod_mixture,traits2)
colnames(pooled)[1:2] <- c("biomass","prod")
pooled$site <- sit
pooled$species <- rownames(pooled)
for (sit in SITE[2:11]){
  SUB <- subset(TOTAL, site==sit & order==orde & simul ==simu)
  SUB2 <- arrange(SUB,factor(SUB$species,levels = traits$SName )) # species in the same order as in the trait data frame
  prod_mixture <- SUB2$prod_mixture
  biomass <- SUB2$mixture.t.ha.
  pooled2 <- data.frame(biomass,prod_mixture,traits2)
  colnames(pooled2)[1:2] <- c("biomass","prod")
  pooled2$site <- sit
  pooled2$species <- rownames(pooled2)
  pooled <- rbind(pooled,pooled2)
}
write.table(pooled,"data/processed/data productivity and biomass~traits_mixture.txt",row.names=F)

# for monocultures ####
# create a data frame with all the monocultures
sit <- SITE[1]
prod_mono <- read.table(paste0("data/processed/productivity_monoculture_",sit,".txt"),header=T)
biom_mono <- read.table(paste0("data/processed/biomass_monoculture_",sit,".txt"),header=T)
pooled_mono <- data.frame(prod_mono$monoculture,biom_mono$monoculture.t.ha.,traits2)
colnames(pooled_mono)[1] <- "prod"
colnames(pooled_mono)[2] <- "biomass"
dat <- pooled_mono
dat$site <- sit
dat$species <- rownames(dat)
for (sit in SITE[2:11]){
  prod_mono <- read.table(paste0("data/processed/productivity_monoculture_",sit,".txt"),header=T)
  biom_mono <- read.table(paste0("data/processed/biomass_monoculture_",sit,".txt"),header=T)
  pooled_mono <- data.frame(prod_mono$monoculture,biom_mono$monoculture.t.ha.,traits2)
  colnames(pooled_mono)[1] <- "prod"
  colnames(pooled_mono)[2] <- "biomass"
  dat2 <- pooled_mono
  dat2$site <- sit
  dat2$species <- rownames(dat2)
  dat <- rbind(dat,dat2)
}

write.table(dat,"data/processed/data productivity~traits_monoculture.txt",row.names = F)
# pooled<-subset(pooled,!(prod_mixture<0)) # pour pouvoir transformer les variables, il me faut une variable réponse positive
 
# par(mfrow=c(1,1)) ; corrplot::corrplot.mixed(cor(pooled_mono), order="hclust", tl.col="black")
# I remove La (highly correlated with Ly), as well as NTol and WiTX (cor with DDmin > 0.5)

# on which data do I work
dat <- subset(pooled_mono,prod>0)
# dat <- pooled_mono

# linear model (12 potential explanatory variables !!)
# CAREFUL: in stressful sites, I sometimes have less that 12 species persisting... 
# ... and consequently more observations than variables.


null <- lm(prod~1,data=dat)
full <- lm(prod ~ S +   HMax + AMax+  G  +  DDMin  +   WiTN  + DrTol    +        
                        Brow + Ly  ,data=dat) # +WiTX + NTol      + La
par(mfrow=c(2,2)) ; plot(full)

step(null, scope = ~ S +   HMax + AMax+  G  +  DDMin  +   WiTN  + DrTol    +
          Brow + Ly ,
     direction="both", criterion = "BIC")
step(full, scope = ~ 1 ,
     direction="both", criterion = "AIC")

GrandeDixence <- lm(formula = prod ~ S + HMax + AMax + G + DDMin + WiTN + DrTol + 
                     Brow + Ly, data = dat)
Bever <- lm(formula = prod ~ S + G + DDMin + WiTN + DrTol + Brow, data = dat)
Davos <- lm(formula = prod ~ G + DDMin + DrTol + Brow + Ly, data = dat)
Adelboden <- lm(formula = prod ~ S + HMax + AMax + DrTol + Ly, data = dat)
Huttwil <- lm(formula = prod ~ S + HMax + AMax + DDMin + DrTol + Brow, data = dat)
Schwerin <- lm(formula = prod ~ S + HMax + G + DrTol + Brow, data = dat)
Bern <- lm(formula = prod ~ S + HMax + AMax + G + DrTol + Brow, data = dat)
Cottbus <- lm(formula = prod ~ S + HMax + G + DDMin + DrTol + Brow + Ly, 
              data = dat)
Basel <- lm(formula = prod ~ S + HMax + DrTol + Brow, data = dat)
Schaffhausen <- lm(formula = prod ~ S + HMax + AMax + DrTol + Brow, data = dat)
Sion <- lm(formula = prod ~ DrTol + Ly, data = dat)



car::Anova(Cottbus)
# MBESS::effect.size(Adelboden)

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
sit <- SITE[1]
for (sit in SITE){
  biom <- read.table(paste0("data/processed/biomass_monoculture_",sit,".txt"),header=T)
  ggplot(data=biom,
         aes(x=Di,y=monoculture.t.ha.,label=SName))+
    geom_point() +
    geom_text(aes(label=SName),hjust=0, vjust=0) +
    ggtitle(sit) #+
    # ggsave(paste0("figures/Correlation distinctiveness biomass/monocultures",sit,".png"))
  # mod <- lm(monoculture.t.ha.^0.01~Di,data=biom)
  # plot(mod)
}

# plot the correlations as a funciton of the sites
correlations <- c()
pval <- c()
for (site in SITE){
  print(site)
  biom <- read.table(paste0("data/processed/biomass_monoculture_",site,".txt"),header=T)
  # biom <- subset(biom,!(SName %in% c("PCem","LDec")))
  correlations <- c(correlations,cor(biom$Di,biom$monoculture.t.ha.,method="kendall") )
  test <- cor.test(biom$Di,biom$monoculture.t.ha.,alternative="two.sided",method = "kendall")
  pval <- c(pval,test$p.value)
} 
#  NB Impossible de calculer la p-value exacte avec des ex-aequos
#  --> Soit faire des permutations, soit se satisfaire de ces p-values
cor <- data.frame(SITE,correlations,pval)
temp_order <- c("GrandeDixence" ,"Bever"     ,    "Davos"      ,   "Adelboden"   ,  "Huttwil"    ,  
                "Schwerin"  ,    "Bern" ,         "Schaffhausen" , "Cottbus"    ,   "Basel"       ,
                "Sion" ) # sites ordered by increasing temperature
cor$SITE <- factor(cor$SITE,levels=temp_order)

ggplot(cor, aes(x=SITE, y=correlations)) +
  geom_boxplot() 


# Correlate biomass and distinctiveness in mixtures ####

sit <- SITE[1]
for (sit in SITE){
  BIOM <- subset(TOTAL, site==sit & order=="increasing" & simul==1)
  ggplot(data=BIOM,
         aes(x=dist,y=mixture.t.ha.,label=species))+
    geom_point() +
    geom_text(aes(label=species),hjust=0, vjust=0) +
    ggtitle(paste0(sit," mixture")) #+
    ggsave(paste0("figures/Correlation distinctiveness biomass/mixture",sit,".png"))
}

# Correlate biomass in mixture and mono ####
sit <- SITE[1]
for (sit in SITE){
  BIOM <- subset(TOTAL, site==sit & order=="increasing" & simul==1)
  ggplot(data=BIOM,
         aes(x=monoculture.t.ha.,y=mixture.t.ha.,label=species))+
    geom_point() +
    geom_text(aes(label=species),hjust=0, vjust=0) +
    ggtitle(paste0(sit)) +
  ggsave(paste0("figures/Correlation distinctiveness biomass/biomass mono_mixt_",sit,".png"))
}

# Correlate biomass and productivity in monocultures ####
sit <- SITE[1]
for (sit in SITE){
  BIOM <- subset(TOTAL, site==sit & order=="increasing" & simul==1)
  # mod <- lm(log(1+BIOM$prod_monoculture)~BIOM$monoculture.t.ha.)
  # plot(mod) # hypothèses non respectées.
  # summary(mod) # la biomasse explique une grande part de la variance dans la productivité
  ggplot(data=BIOM,
         aes(x=monoculture.t.ha.,y=prod_monoculture,label=species))+
    geom_point() +
    geom_text(aes(label=species),hjust=0, vjust=0) +
    ggtitle(paste0(sit)) +
    ggsave(paste0("figures/Correlation distinctiveness biomass/biomass_prod_mono_",sit,".png"))
}

# Correlate biomass and productivity in mixture ####
sit <- SITE[1]
for (sit in SITE){
  BIOM <- subset(TOTAL, site==sit & order=="increasing" & simul==1)
  ggplot(data=BIOM,
         aes(x=mixture.t.ha.,y=prod_mixture,label=species))+
    geom_point() +
    geom_text(aes(label=species),hjust=0, vjust=0) +
    ggtitle(paste0(sit)) +
    ggsave(paste0("figures/Correlation distinctiveness biomass/biomass_prod_mixt_",sit,".png"))
}


# Look at the distinctiveness of the remaining species ####

# compute functional distinctiveness
comp_fct_dist <- function(traits){
  c1<-choice_traits_1(traits) # data.frame with the traits of the species
  c1 <- select(c1,-c(A1max,A2))# remove traits that say explicitely that a species is a gymnosperm
  ACP1<-PCA(c1,graph=F)
  distACP <- ACP1$ind$coord %>% 
    as.data.frame() %>%
    select(Dim.1    ,   Dim.2     ,   Dim.3      ,  Dim.4) %>%
    traits_dist() # here we compute the functional distinctiveness
  distACP$SName<-rownames(distACP)
  
  distinct_tot <-  distACP
  distinct_tot[order(distinct_tot$Di),]
}


correlation <- c()
Di <- comp_fct_dist(traits)
orDi <- Di[order(Di$Di,decreasing=T),] # remove distinct species last
for (i in c(1:30)){
  subDi <- orDi[i:30,] # remaining species after having removed i-1 species 
  subtraits <- subset(traits,as.character(SName) %in% subDi$SName) # traits of these remaining species
  newDi <- comp_fct_dist(subtraits) # compute fct dist on these remaining species
  subDi <- subset(Di,SName %in% newDi$SName) # keep the original distinctiveness (calculated relatively to all the species) of these remaining species
  A <- newDi[str_order(newDi$SName),] # I need the species to be in the same order to compute spearman's correlation
  B <- subDi[str_order(subDi$SName),]
  cor <- cor(A$Di,B$Di,method="spearman") # compare how distinctiveness sorts the species in the two cases
  # remember that the order with all the species is the one we used for the removal experiments. Is it the same as that of a local distinctiveness ?
  # it is not exactly matthias' question, which would require taking the ACTUAL species remaining at the end of the simulations.
  correlation <- c(correlation,cor)
}
plot(c(1:length(correlation)),correlation)

# fct dist for remaining species ####
Di <- comp_fct_dist(traits)
orde <- "decreasing"
prod <- read.table("data/processed/productivity_specific_GrandeDixence_with monocultures.txt",header=T)
correlation <- c()
for (simu in c(3,5,7,9,11,13,15,17,19,21,23,25,27)){
  sub <- subset(prod,order==orde & simul==simu) # remaining species
  sub <- sub[which(sub$mixture_t_ha>0),]
  subtraits <- subset(traits,as.character(SName) %in% sub$species) # traits of these remaining species
  newDi <- comp_fct_dist(subtraits) # compute fct dist on these remaining species
  subDi <- subset(Di,SName %in% newDi$SName) # keep the original distinctiveness (calculated relatively to all the species) of these remaining species
  A <- newDi[str_order(newDi$SName),] # I need the species to be in the same order to compute spearman's correlation
  B <- subDi[str_order(subDi$SName),]
  cor <- cor(A$Di,B$Di,method="spearman") # compare how distinctiveness sorts the species in the two cases
  correlation <- c(correlation,cor)
}
plot(2*c(1:length(correlation)),correlation,xlab="Number of species removed",ylab="Correlation")

sub[order(sub$species),] # productivity of the species
newDi[order(newDi$SName),]

# look at the correlation between distinctiveness and productivity for these persisting species
NEW <- newDi[order(newDi$SName),]
NEW$prod <- sub[order(sub$species),]$mixture_t_ha
plot(NEW$Di,NEW$prod)


# Local pool_monocultures_recompute fct distinctiveness ####
sit <- "GrandeDixence"
sit <- "Bever"
sit <- "Davos"

Di <- comp_fct_dist(traits)
cor_site <- c()
nb_sp <- c()
for (sit in SITE){
  current_site <- read.table(paste0("data/processed/biomass_monoculture_",sit,".txt"),header=T)
  sp <- as.character(subset(current_site,monoculture.t.ha.>0)$SName)
  # if I want to do it for species that persist, and not for monocultures
  # current_site <- read.table(paste0("data/processed/productivity_specific_",sit,".txt"),header=T)
  # sp <- as.character(subset(current_site,mixture_t_ha>0 & order=="increasing" & simul==1)$species)
  
  traits3 <- subset(traits,SName %in% sp)
  newDi <- comp_fct_dist(traits3) # compute fct dist on these remaining species
  subDi <- subset(Di,SName %in% newDi$SName) # keep the original distinctiveness (calculated relatively to all the species) of these remaining species
  A <- newDi[str_order(newDi$SName),] # I need the species to be in the same order to compute spearman's correlation
  B <- subDi[str_order(subDi$SName),]
  coefcor <- cor(A$Di,B$Di,method="spearman") # compare how distinctiveness sorts the species in the two cases
  cor_site <- c(cor_site, coefcor)
  nb_sp <- c(nb_sp,length(sp))
}

data.frame(SITE,cor_site,nb_sp)

# pour ce qui suit, ne pas faire tourner la boucle, mais exécuter le code pour seulement un site
cor.test(A$Di,B$Di,alternative="two.sided",method = "spearman")

biom <- read.table(paste0("data/processed/biomass_monoculture_",sit,".txt"),header=T)
abdist <- data.frame(newDi$Di,newDi$SName) # se if the abundance/dist relationships holds locally (with monocultures)
colnames(abdist) <- c("Di","SName")
ab <- c()
for (i in c(1:dim(abdist)[1])){
  ab <- c(ab,biom[which(as.character(biom$SName)==as.character(abdist[i,]$SName)),]$monoculture.t.ha.)
}
abdist$abundance_mono <- ab
plot(abdist$Di,abdist$abundance)
abdist

# A Bever, la corrélation n'est que de 0.42... Donc pas super maintien de l'ordre de distinctiveness pour le pool local.
A1 <- A[1:floor(dim(A)[1]/2),] # se if the second half (most distinct species) still correlates or not
A2 <- A[(floor(dim(A)[1]/2)+1):dim(A)[1],]
B1 <- B[1:floor(dim(B)[1]/2),]
B2 <- B[(floor(dim(B)[1]/2)+1):dim(B)[1],]
cor(A1$Di,B1$Di,method="spearman") # 0.23
cor(A2$Di,B2$Di,method="spearman") # 0.81 à Bever
# It seems that the most distinct species' order does not change must when we compute distinctiveness on the local pool
# For unstressful sites, the correlation is perfect, because we keep all the species (they can all grow in monoculture)
# In stressful sites : either the correlation is not so nad (GrandeDixence, 0.65)
# or it is quite bad, but is it good for the distinct half of the species (cf. Bever, 0.4 but 0.81 for the distinct half)

# sensib_fct dist_traits ####
# traits2 is the data frame of traits
selected_traits<-choice_traits_1(traits) # data.frame with the traits of the species
# selected_traits <- select(selected_traits,-c(A1max,A2))# remove traits that say explicitely that a species is a gymnosperm

comp_dist_selected_traits <- function(selected_traits){
  # This function computes fct distinctinveness on a trait matrix containing only the 12 traits
  # that were used for ranking the species.
  # I use it for bootstraping on these 12 traits
  ACP1<-PCA(selected_traits,graph=F)
  distACP <- ACP1$ind$coord %>% 
    as.data.frame() %>%
    select(Dim.1    ,   Dim.2     ,   Dim.3      ,  Dim.4) %>%
    traits_dist() # here we compute the functional distinctiveness
  distACP$SName<-rownames(distACP)
  
  distinct_tot <-  distACP
  distinct_tot[order(distinct_tot$SName),]$Di
}


init_dist <- comp_dist_selected_traits(selected_traits)

corstat <- c()
for (i in c(1:1000)){
  col <- sample(c(1:12),size=12,replace=T) # randomly select 12 columns to keep (with replacement) from the trait matrix
  boot_traits <- selected_traits[,col]
  boot_dist <- comp_dist_selected_traits(boot_traits)
  boot_cor <- cor(init_dist,boot_dist,method="spearman")
  corstat <- c(corstat,boot_cor)
}
summary(corstat)
sd(corstat)

jpeg("paper/histogram bootstrap.jpg", width = 350, height = 350)
hist(corstat,main="Bootstrap distribution of rho", xlab="Rho")
dev.off()


# compute functional distinctiveness on traits not related to envt response (Annette's remark) ####
traits.simulations <- select(traits,Name,SName,S, HMax, AMax,   G, DDMin, WiTN, WiTX, DrTol, NTol, Brow,   Ly, La,A1max, A2)
dist.all.traits <- comp_fct_dist(traits.simulations)

traits.not.envt <- select(traits.simulations,Name,SName,S, HMax, AMax,   G, Brow,   Ly, La,A1max, A2)
dist.not.envt <- comp_fct_dist(traits.not.envt)

plot(dist.all.traits$Di,dist.not.envt$Di)
abline(0,1)

data.frame(dist.all.traits$SName,dist.not.envt$SName)
cor(dist.all.traits$SName,dist.not.envt$SName)

traits.simulations <- select(traits,S, HMax, AMax,   G, DDMin, WiTN, WiTX, DrTol, NTol, Brow,   Ly, La,A1max, A2)
traits.not.envt <- select(traits.simulations,S, HMax, AMax,   G, Brow,   Ly, La,A1max, A2)
a <- comp_dist_selected_traits(traits.simulations)
b <- comp_dist_selected_traits(traits.not.envt)
cor(a,b,method="spearman")

# prod added to PCA on traits ####
# I will add columns to traits 2. 
# In each column, I put the productivity of each species in monoculture in a given site
traits2[order(rownames(traits2)),]
for (sit in SITE){
  prod <- read.table("data/processed/productivity_specific_GrandeDixence_with monocultures.txt",header=T)
  PR <- read.table(paste0("data/processed/productivity_monoculture_",sit,".txt"),header=T)
  traits2 <- cbind(traits2,PR[order(PR$SName),]$monoculture)
}
colnames(traits2)[c(13:23)] <- SITE

# then I can make a PCA, on which plot the supplementary variables (for each site,  
# the productivity of the species)
res.pca <- PCA(traits2, ind.sup = NULL, 
               quanti.sup = 13:23, quali.sup = NULL, graph=FALSE)
fviz_pca_var(res.pca)


traits2.active <- traits2[,1:12]
library("ade4")
PCA(traits2[,c(1:12)],graph=T)
s.arrow(traits2[,c(13:23)])



