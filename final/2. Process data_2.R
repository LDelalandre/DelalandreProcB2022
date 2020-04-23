source("R/Common variables.R")
source("R/Analysis_data.R")
source("R/Monocultures_functions.R")

# CAREFUL not to run all the code: it includes processing data and writing tables, which:
# - takes a lot of time
# - can overwrite previous data frames

# Biomass per species_removal experiments ####
for (site in SITE){ # Write a table with the specific final biomasses of all sets of simulations for each site
  # NB: Takes a long time to run!!
  BIOMASSES <- biomass_specific(site=site) 
  write.table(BIOMASSES,paste0("data/processed/biomass_specific_",site,".txt"),sep="\t",row.names=F)
}

# Productivity per species_removal experiment ####
for (site in SITE){
  PROD <- productivity_specific(site)
  write.table(PROD,paste0("data/processed/productivity_specific_",site,".txt"),sep="\t",row.names=F)
}

# Monocultures: biomass per species ####
for (site in SITE){
  specific_val <- specific_biomasses(site)
  write.table(specific_val,paste0("data/processed/biomass_monoculture_",site,".txt"),row.names=F,sep="\t")
  biomass <- biomasses(site = site,specific_val = specific_val)
  write.table(biomass,paste0("data/processed/biomass_specific_",site,"_with monocultures.txt"),row.names=F,sep="\t")
}

# Monocultures:  productivities per species ####
for (site in SITE){
  specific_val <- specific_productivities(site)
  write.table(specific_val,paste0("data/processed/productivity_monoculture_",site,".txt"),row.names=F,sep="\t")
  prod <- productivities(site = site,specific_val = specific_val)
  write.table(prod,paste0("data/processed/productivity_specific_",site,"_with monocultures.txt"),row.names=F,sep="\t")
}

# Pool all the data ####
# First: productivity: make a real subset of the species above biomass threshold
TOTAL <- c()
for (site in SITE[]){
  biom <- read.table(paste0("data/processed/biomass_specific_",site,"_with monocultures.txt"),header=T)
  prod <- read.table(paste0("data/processed/productivity_specific_",site,"_with monocultures.txt"),header=T)
  for (ord in as.character(unique(biom$order))){
    sub_biomass <- subset(biom,order==ord)
    for (nb in unique(sub_biomass$simul)){
      sub_biom <- subset(sub_biomass,simul==nb)
      sub_prod <- subset(prod,order==ord & simul == nb)
      to_keep <- subset(sub_prod,species %in% sub_biom$species)
      to_keep$monoculture_relative <- to_keep$monoculture/sum(to_keep$monoculture)
      to_keep$mixture_relative <- to_keep$mixture_t_ha/sum(to_keep$mixture_t_ha)
      sub_biom <- arrange(sub_biom,species)
      to_keep <- arrange(to_keep,species) # so that I am sure that they are in the same order
      sub_biom$prod_mixture <- to_keep$mixture_t_ha
      sub_biom$prod_mixture_relative <- to_keep$mixture_relative
      sub_biom$prod_monoculture <- to_keep$monoculture
      sub_biom$prod_monoculture_relative <- to_keep$monoculture_relative
      TOTAL <- rbind(TOTAL,sub_biom)
    }
  }
}
write.table(TOTAL,"data/processed/specific_biom_prod_complete.txt",row.names=F)
# this is a data frame with all the biomass and productivity for all the species
# then I need to group for all the simul, using +/- hte following functions:
# BIOMASSES_tot <- aggregate(BIOMASSES_sp$mixture.t.ha.,list(site = BIOMASSES_sp$site,order = BIOMASSES_sp$order,simul = BIOMASSES_sp$simul),FUN=sum)# sum of the biomasses of the species
# biomass_per_order <- spread(BIOMASSES_tot,order,x) # data frame with each order in column

# Data frames to plot species removal experiments ####
for (sit in SITE){ 
  BIOMASSES_sp <- subset(TOTAL,site==sit)
  for (measure in MEASURE){
    if (measure=="biomass_tot"){
      BIOMASSES_tot <- aggregate(BIOMASSES_sp$mixture.t.ha.,
                                 list(site = BIOMASSES_sp$site,order = BIOMASSES_sp$order,simul = BIOMASSES_sp$simul),
                                 FUN=sum)
    } else {
      BIOMASSES_tot <- aggregate(BIOMASSES_sp$prod_mixture,
                                 list(site = BIOMASSES_sp$site,order = BIOMASSES_sp$order,simul = BIOMASSES_sp$simul),
                                 FUN=sum)
    }
    biomass_per_order <- spread(BIOMASSES_tot,order,x)
    write.table(biomass_per_order,paste0("data/processed/",measure,"_",sit,".txt"),sep="\t",row.names=F)
    
    conf <- confidence_interval(sit,measure)
    write.table(conf,paste0("data/processed/",measure,"_",sit,"_with interval.txt"),sep="\t",row.names=F)
  }
}



