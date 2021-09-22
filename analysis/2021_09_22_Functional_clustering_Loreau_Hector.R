source("R/Analysis_data.R")
source("R/Common variables.R")
library("functClust")
library(tidyverse)

if(file.exists(file.path("data/processed/functClust"))==F){
  dir.create(file.path("data/processed/functClust"))  
}

LH_all_per_sp <- read.csv("data/processed/Loreau-Hector coefficients_per species.csv") %>% 
  mutate(site = factor(site,levels=SITE))

# Table format functClust data ####
# Je veux :
#   - lire les propriétés écosystémiques finales (cf. mes autres codes R)
#   - savoir les espèces présentes (au-delà d'un seuil) finales
#   - mettre tout ça dans un tableau du format d'entrée de functclust
#   - je peux aussi intégrer toutes les monocultures dedans. Demander à Benoît.


table_functClust <- function(LH_all_per_sp,sit){
  df_LH <- LH_all_per_sp %>% 
    filter(site == sit)
  
  dat <- data.frame(matrix(nrow=0,ncol=35))
  colnames(dat) <- c("assemblage",
                     as.character(read.table("data/raw/distinctiveness of the species.txt",header = T)$SName),
                     "productivity","DeltaY","Selection","Complementarity")
  j <- 1 # number of the assemblage = position (line) in the data.frame dat

  for (order in ORDER){
    a <- order
    df_LH_order <- subset(df_LH,order==a)

    for (i in df_LH_order$simul %>% unique() %>% sort() ){ # the removal experiments where there were species remaining at the end
      sub <- subset(df_LH_order, simul==i) 

      dat[j,]$assemblage <- paste0(site,"_",order,"_",i) # paste order and simul (an assemblage is a simul)
      
      dat[j,]$productivity <- sum(sub$YOi)
      dat[j,]$DeltaY <- unique(sub$DeltaY)
      dat[j,]$Selection <- unique(sub$Selection)
      dat[j,]$Complementarity <- unique(sub$Complementarity)
      
      dat[j,which(!(colnames(dat) %in% c(as.character(sub$SName),"assemblage","productivity","DeltaY","Selection","Complementarity") ))] <- 0
      dat[j,which(colnames(dat) %in% sub$SName)] <- 1
      
      j <- j+1 # current position in the data.frame dat.
    }
  }
  # Above was the code for monocultures - but there is no selection nor complementarity effect in monoculture. 
  # I thus remove it.
  
  # for (k in 1:dim(table_mono)[1]){
  #   dat[j,]$assemblage <- paste0("monoculture_",site,"_",k)
  #   dat[j,]$biomass <- table_mono[k,]$monoculture.t.ha.
  #   dat[j,]$productivity <- prod_mono[k,]$monoculture
  #   dat[j,]
  #   
  #   dat[j,which(!(colnames(dat) %in% c(as.character(table_mono[k,]$SName),"biomass","assemblage","productivity") ))] <- 0
  #   dat[j,which(colnames(dat) %in% table_mono[k,]$SName)] <- 1
  #   
  #   j <- j+1
  # }
  dat
}

for (site in SITE){
  dat <- table_functClust(LH_all_per_sp,site)
  write.table(dat,paste0("data/processed/functClust/functClust_data_",site,".txt"),sep="\t",row.names=F)
}

#_______________________________________________________________________________
# functClust analysis ####
PROPERTIES <- c("productivity","DeltaY","Selection","Complementarity")

for (site in SITE){
  dat <- read.table(paste0("data/processed/functclust/functClust_data_",site,".txt"),header=T)
  
  # Number of species persisting
  dat <- distinct(dat) # Remove duplicated rows.
  zeros <- which(colSums(dat[,-c(1,32,33,34,35)])==0)
  if (length(zeros)>0){
    dat2 <- dat[,-(zeros+1)] # Remove absent species, and ecosystem properties without value (i.e. all values = 0)
  } else {
    dat2 <- dat
  }
 
  nbElt <- dat2 %>% 
    select(-c("assemblage","productivity","DeltaY","Selection","Complementarity")) %>% 
    colnames() %>% 
    length()
  
  # Chose the ecosystem property to study
  for (property in PROPERTIES){
    focus <- which(PROPERTIES == property)
    properties_discarded <- PROPERTIES[-focus]
    # I remove the other 3 properties
    dat3 <- dat2 %>% 
      select(-rlang::UQ(sym(properties_discarded[1])),
             -rlang::UQ(sym(properties_discarded[2])),
             -rlang::UQ(sym(properties_discarded[3])) )
    
    res1 <- fclust(dat3,nbElt)
    
    save(res1, file = paste0("data/processed/functclust/",site,"_",property,".RData"))
  }
}





# PLot ####
for (site in SITE[4:11]){
  for (property in PROPERTIES){
    load( paste0("data/processed/functclust/",site,"_",property,".RData"))
    
    nbcl=2
    pdf(paste0("figures/functClust/functClust_",site,"_nbcl=",nbcl,"_",property ,".pdf"),width=10)
    fclust_plot(res1,nbcl=nbcl)
    dev.off()
  }
}





