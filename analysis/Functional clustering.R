source("R/Analysis_data.R")
source("R/Common variables.R")
library("functClust")
library(dplyr)

# Table format functClust data ####
# Je veux :
#   - lire les propriétés écosystémiques finales (cf. mes autres codes R)
#   - savoir les espèces présentes (au-delà d'un seuil) finales
#   - mettre tout ça dans un tableau du format d'entrée de functclust
#   - je peux aussi intégrer toutes les monocultures dedans. Demander à Benoît.


table_functClust <- function(site){
  dat <- data.frame(matrix(nrow=0,ncol=33))
  colnames(dat) <- c("assemblage",as.character(read.table("data/raw/distinctiveness of the species.txt",header = T)$SName),"biomass","productivity")
  j <- 1 # number of the assemblage = position (line) in the data.frame dat
  table <- read.table(paste0("mass grave/data processed 2020_28_12/biomass_specific_",site,"_with monocultures.txt"),header=T)
  table_mono <- read.table(paste0("mass grave/data processed 2020_28_12/biomass_monoculture_",site,".txt"),header=T)
  prod <- read.table(paste0("mass grave/data processed 2020_28_12/productivity_specific_",site,"_with monocultures.txt"),header=T)
  prod_mono <- read.table(paste0("mass grave/data processed 2020_28_12/productivity_monoculture_",site,".txt"),header=T)
  for (order in ORDER){
    a <- order
    sub1 <- subset(table,order==a)
    pr1 <- subset(prod,order==a)
    for (i in unique(sub1$simul)){ # the removal experiments where there were species remaining at the end
      sub <- subset(sub1, simul==i) # biomass of each species in one site, one order, one simul
      pr <- subset(pr1,simul==i) # productivity of each species in one site, one order, one simul
      
      dat[j,]$assemblage <- paste0(site,"_",order,"_",i)
      dat[j,]$biomass <- sum(sub$mixture.t.ha)
      dat[j,]$productivity <- sum(pr$mixture_t_ha)
      dat[j,which(!(colnames(dat) %in% c(as.character(sub$species),"biomass","assemblage","productivity") ))] <- 0
      dat[j,which(colnames(dat) %in% sub$species)] <- 1
      
      j <- j+1 # current position in the data.frame dat.
    }
  }
  for (k in 1:dim(table_mono)[1]){
    dat[j,]$assemblage <- paste0("monoculture_",site,"_",k)
    dat[j,]$biomass <- table_mono[k,]$monoculture.t.ha.
    dat[j,]$productivity <- prod_mono[k,]$monoculture
    dat[j,]
    
    dat[j,which(!(colnames(dat) %in% c(as.character(table_mono[k,]$SName),"biomass","assemblage","productivity") ))] <- 0
    dat[j,which(colnames(dat) %in% table_mono[k,]$SName)] <- 1
    
    j <- j+1
  }
  dat
}

for (site in SITE){
  dat <- table_functClust(site)
  write.table(dat,paste0("data/processed/functClust/functClust_data_",site,".txt"),sep="\t",row.names=F)
}

# functClust analysis ####
site <- "Bever"

dat <- read.table(paste0("mass grave/data processed 2020_28_12/functclust/functClust_data_",site,".txt"),header=T)
nbElt <- 30
dat <- distinct(dat) # Romove duplicated rows.

dat2 <- dat[,-(which(colSums(dat[,-1])==0)+1)] # Remove absent species, and ecosystem properties without value (i.e. all values = 0)
nbElt2 <- nbElt-length(which(colSums(dat[,2:31])==0))# count how many species remain in the data frame
# dat2[which(is.na(rowSums(dat2[,-1])==0)),]

res1 <- fclust(dat,nbElt)
res2 <- fclust(dat2, nbElt2)

nbcl=3
pdf(paste0("figures/functClust/functClust_",site,"_nbcl=",nbcl,".pdf"),width=10)
plot_fclust(res1,nbcl=nbcl)
dev.off()

save(res1, file = paste0("data/processed/functclust/",site,".RData"))
load( paste0("data/processed/functclust/",site,".RData"))

# sans les monocultures: ne semble pas changer grand-chose
# res2 <- fclust(select(dat2,-c(biomass_mono)), nbElt2)
# pdf(paste0("figures/functClust/functClust_",site,"_nbcl=",nbcl,"_sans monocultures.pdf"),width=10)
# plot_fclust(res2,nbcl=nbcl)
# dev.off()

# Plot Distinctiveness f(functional cluster)

# Composition and interaction effects ####
mOccur <- as.matrix(select(dat,-c(assemblage,biomass)))
fobs <- as.vector(select(dat,biomass))
mult <- multiplicative_decomposition(fobs, mOccur, rm.mono = F)
part_effects <- cbind(mult$alpha,mult$beta)
colnames(part_effects) <- c("interaction","composition")
hist(part_effects$)