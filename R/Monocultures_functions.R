source("R/Common variables.R")

specific_values <- function(site){
  # Extract the biomasse of each species in a monoculture and make a data.frame with the following columns:
  # Di SName biomass_monoculture     sd_biom Id relative_biomass
  
  # Rq: be careful of the way the folders and files are named.
  # NB: add productivity and sd(productivity) to it when I have the adequate simulations!
  biomass<-c()
  sd_biom<-c()
  for(i in c(1:30)){
    mean <- read.table(paste0("data/raw/output-cmd2_",site,"_monoculture.txt/forceps.",site,".site_",i,"_mean.txt"))
    # NB it is in the good order because we made the monocultures in the order of the species Id number.
    colnames(mean)<-colnames_mean
    # biomass averaged
    years_to_keep <- max(mean$date) - c(900,800,700,600,500,400,300,200,100,0)
    meanbiom <- mean(subset(mean,date %in% years_to_keep)$totalBiomass.t.ha.)
    biomass<-c(biomass,meanbiom)
    sd_biom<-c(sd_biom,sd(mean$totalBiomass.t.ha.))
  }
  
  # specific_values: a data frame with final biomass of each monoculture, etc. 
  specific_values <- read.table("data/raw/distinctiveness of the species.txt",header=T)
  specific_values$'monoculture(t/ha)' <- biomass
  specific_values$sd_biom <- sd_biom
  specific_values$Id <- c(0:29)
  specific_values$monoculture_relative <- specific_values$'monoculture(t/ha)' / sum(specific_values$'monoculture(t/ha)')
  specific_values
}


biomasses <- function(site,specific_val){
  # returns a data frame with colnames:
  # species mixture.t.ha. mixture_relative site      order simul monoculture(t/ha) monoculture_relative       Di
  # Each species that apears is a species that is present at the end of a mixture simulation. Consequently, we don't have all the species.
  # NB: 
  
  # Example for the variables: site="Bern" ; specific_val = specific_values(site) ; number=1 ; order="decreasing" ; Nbpatches = 10 (ou 50)
  # Number is the number of the simulation (1: we didn't remove any species, until 30: there is only one species left)
  biomasses <- read.table(paste0("data/processed/specific_biomass_final_",site,".txt"),header=T)
  
  # Add a column with absolute biomass of the same species in monoculture
  biom_mono <- c()
  dist <- c()
  for (i in 1:dim(biomasses)[1]){
    biom_mono <- c(biom_mono,specific_val[which(specific_val$SName == as.character(biomasses[i,]$species) ),]$'monoculture(t/ha)')
    dist <- c(dist,specific_val[which(specific_val$SName == as.character(biomasses[i,]$species) ),]$Di)
  }
  biomasses$`monoculture(t/ha)` <- biom_mono
  biomasses$dist <- dist
  
  # Add a column with relative biomass of the same species in a mix of monocultures with the same species
  mono_sum <- aggregate(biomasses$`monoculture(t/ha)`, list(order=biomasses$order,simul=biomasses$simul), sum) # mean per species
  biom_relat <- c()
  for (i in 1:dim(biomasses)[1]){
    summed_biomass_mono <- mono_sum[which(mono_sum$order==biomasses[i,]$order & mono_sum$simul==biomasses[i,]$simul),]$x
    biom_relat[i] <- biomasses[i,]$`monoculture(t/ha)`/summed_biomass_mono
  }
  biomasses$monoculture_relative <- biom_relat
  
  # vérification : somme des biomasses relatives = 1. Ok.
  # aggregate(biomasses$monoculture_relative,list(order=biomasses$order,simul=biomasses$simul),FUN=sum)
  # sum(subset(biomasses,simul==1&order=="decreasing")$monoculture_relative)
  biomasses
}

