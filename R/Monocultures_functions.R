source("R/Common variables.R")



add_zeros <- function(site,productivity_specif,measure){
  sit <- site
  # This function adds lines with species present at the beginning of a simulation (=regional pool) but
  # absent in the realized pool (=in mixture)
  # Useful for the Loreau-Hector analysis.
  data <- productivity_specif
  data$species <- as.character(data$species)
  data$site <- as.character(data$site)
  data$order <- as.character(data$order)
  # Include species whose biomass or productivity = 0
  
  # First, have a data frame with the Id of species in all the oders of removal
  SPord <- read.table("data/removal orders.txt",header=T)
  
  # Second, have the correspondence between species Id and species short name
  name_sname <- read.table("data/correspondence_SName_Id.txt",header=T)
  
  # Third, add lines in the file data with 
  for (sim in c(1:30)){# all the simul
    for (ord in ORDER){
      spId <- SPord[,which(colnames(SPord) == ord)] # list of sp Id ordered for a given order of removal
      pr_in <- spId[sim:length(spId)] # list of the Id of species that were present at the beginning of that simulation
      zeros <- name_sname$SName[which(name_sname$Id %in% pr_in)] # Short names of the same species
      sub <- subset(data,order==ord&simul==sim & site==sit)
      to_add <- zeros[which(!(zeros %in% sub$species))] # species that were present at the beginning of the simul, but not at the end.
      for(sp in to_add){
        data <- rbind(data,c(sp,0,0,0,sit,ord,sim))
      }
    }
    
  }
  data
}

biomasses <- function(specif.biomass,specific_val){
  # returns a data frame with colnames:
  # species mixture.t.ha. mixture_relative site      order simul monoculture(t/ha) monoculture_relative       Di
  # Each species that apears is a species that is present at the end of a mixture simulation. Consequently, we don't have all the species.
  # NB: 
  
  # Example for the variables: site="Bern" ; specific_val = specific_values(site) ; number=1 ; order="decreasing" ; Nbpatches = 10 (ou 50)
  # Number is the number of the simulation (1: we didn't remove any species, until 30: there is only one species left)
  biomasss <- specif.biomass
  
  biomasses <- add_zeros(site=site,biomasses=biomasss,measure="biomass_tot")
  # Add a column with absolute biomass of the same species in monoculture
  biom_mono <- c()
  dist <- c()
  ab <- c()
  for (i in 1:dim(biomasses)[1]){
    biom_mono <- c(biom_mono,specific_val[which(specific_val$SName == as.character(biomasses[i,]$species) ),]$'monoculture(t/ha)')
    dist <- c(dist,specific_val[which(specific_val$SName == as.character(biomasses[i,]$species) ),]$Di)
    ab <- c(ab,specific_val[which(specific_val$SName == as.character(biomasses[i,]$species) ),]$abundance)
  }
  biomasses$`monoculture(t/ha)` <- biom_mono
  biomasses$dist <- dist
  biomasses$mono_abundance <- ab
  
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




productivities <- function(site,specific_val){
  # Adds a column with productivity in monoculture to data frames with prod in mixture
  # and adds lines with species absent from the local pool but present in the regional (their prod is zero)
  # returns a data frame with colnames:
  # species productivity.t.ha. mixture_relative site      order simul monoculture(t/ha) monoculture_relative       Di
  # Each species that appears is a species that is present at the end of a mixture simulation. Consequently, we don't have all the species.
  
  # Example for the variables: site="Bern" ; specific_val = specific_values(site)
  # Number is the number of the simulation (1: we didn't remove any species, until 30: there is only one species left)
  productivityy <- read.table(paste0("data/processed/productivity_specific_",site,".txt"),header=T)
  productivity <- add_zeros(site=site,biomasses=productivityy,measure="productivity_tot")
  # Add a column with absolute biomass of the same species in monoculture
  biom_mono <- c()
  dist <- c()
  for (i in 1:dim(productivity)[1]){
    biom_mono <- c(biom_mono,specific_val[which(specific_val$SName == as.character(productivity[i,]$species) ),]$monoculture)
    dist <- c(dist,specific_val[which(specific_val$SName == as.character(productivity[i,]$species) ),]$Di)
  }
  productivity$monoculture <- biom_mono
  productivity$dist <- dist
  
  # Add a column with relative productivity of the same species in a mix of monocultures with the same species
  mono_sum <- aggregate(productivity$monoculture, list(order=productivity$order,simul=productivity$simul), sum) # mean per species
  biom_relat <- c()
  for (i in 1:dim(productivity)[1]){
    summed_biomass_mono <- mono_sum[which(mono_sum$order==productivity[i,]$order & mono_sum$simul==productivity[i,]$simul),]$x
    biom_relat[i] <- productivity[i,]$monoculture/summed_biomass_mono
  }
  productivity$monoculture_relative <- biom_relat
  
  # vérification : somme des biomasses relatives = 1. Ok.
  # aggregate(biomasses$monoculture_relative,list(order=biomasses$order,simul=biomasses$simul),FUN=sum)
  # sum(subset(biomasses,simul==1&order=="decreasing")$monoculture_relative)
  productivity
}
