source("R/Common variables.R")
source("R/Analysis_data.R")
source("R/Monocultures_functions.R")

# CAREFUL not to run all the code: it includes processing data and writing tables, which:
# - takes a lot of time
# - can overwrite previous data frames

# Biomass per species_removal experiments ####
for (site in SITE){ # Write a table with the specific final biomasses of all sets of simulations for each site
  # NB: does not process data for temperature increase simulations
  # NB: Takes a long time to run!!
  BIOMASSES <- NULL
  for (order in ORDER){
  BIOMASSES <- rbind(BIOMASSES,biomass_specific(site=site,order=order))
  }
  write.table(BIOMASSES,paste0("data/processed/biomass_specific_",site,".txt"),sep="\t",row.names=F)
}

# PRODUCTIVITY ####
persist <- function(sp,persistent.sp){
  # returns true if and only if species sp is present in the vector persisten.sp
  any(grepl(sp,persistent.sp))
}

# Specific productivity ####
ALL <- 
  read.table("data/processed/biomass_specific_ALL sites.txt",header=T)

for (sit in SITE){
  PROD <- NULL
  for (orde in ORDER){
    for(number in c(1:30)){
      # Compute productivity per species
      prod <- productivity_specific(site=sit,order=orde,number=number)
      
      # Vector of species above biomass threshold
      persistent.sp <- ALL %>% 
        filter(site==sit & order==orde & simul==number) %>% 
        pull(speciesShortName)
      
      # Add a column indicating if the given species persists
      prod$persists <- purrr::map_lgl(prod$species,persist,persistent.sp)
      
      PROD <- rbind(PROD,prod)
    }
  }
  write.table(PROD,paste0("data/processed/productivity_specific_",sit,".txt"),sep="\t",row.names=F)
}

# Community-level productivity ####
for (sit in SITE){
  prod <- 
    read.table(paste0("data/processed/productivity_specific_",sit,".txt"),header=T) %>% 
    filter(persists==T) %>% # keep species above biomass threshold
    group_by(site,order,simul) %>% 
    summarize(productivity=sum(mixture_t_ha)) %>% # compute community-level productivity
    spread(order,productivity) # spread the data frame to plot it
  prod[is.na(prod)] <- 0
  
  prod2 <- median_conf_int(as.data.frame(prod))
  write.table(prod2,paste0("data/processed/",measure,"_",sit,"_with interval_median.txt"),row.names=F)
}

# MONOCULTURES ####
# Monocultures: biomass per species ####
for (site in SITE){
  # Biomass per species in monoculture
  specific_val <- specific_biomasses(site)
  write.table(specific_val,paste0("data/processed/biomass_monoculture_",site,".txt"),row.names=F,sep="\t")
  
  # Biomass per species in mixture and monoculture in the same data frame (and including zeros)
  specif.biomass <- read.table(paste0("data/processed/biomass_specific_",site,".txt"),header=T)
  biomass <- biomasses(specif.biomass=specif.biomass,specific_val = specific_val)
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
TOTAL <- read.table("data/processed/specific_biom_prod_complete.txt",header=T)
for (sit in SITE){ 
  BIOMASSES_sp <- subset(TOTAL,site==sit)
  for (measure in c("biomass_tot","productivity_tot")){
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
  }
}

# Add median confidence interval to biomass and productivity removal experiments ####
for (sit in SITE){
  for (measure in c("biomass_tot","productivity_tot")){
    table <- read.table(paste0("data/processed/",measure,"_",sit,".txt"),header=T)
    table3 <- median_conf_int(table)
    
    write.table(table3,paste0("data/processed/",measure,"_",sit,"_with interval_median.txt"))
    
  }
}


# TEMPORAL STABILITY ####

# Specific temporal stability ####
for (site in SITE){
  TS_prod_sp <- NULL
  for (order in ORDER){
    TS_prod_sp=rbind(TS_prod_sp,sd_productivity_specific(site=site,order=order))
  }
  write.table(TS_prod_sp, paste0("data/processed/TS_productivity_specific_",site,".txt"),sep="\t")
}

# Community-level temporal stability ####
for (site in SITE){ # Write TS and sd
  TS_prod_tot <- NULL
  for (order in ORDER){
    TS_prod_tot <- rbind(TS_prod_tot,sd_productivity_tot(site))
  }
  write.table(TS_prod_tot, paste0("data/processed/TS_sd_productivity_tot_",site,".txt"),sep="\t",row.names = F)
}

for (site in SITE){ # Write TS in a format ready to plot (spread the data)
  TS_prod_tot <- read.table(paste0("data/processed/TS_sd_productivity_tot_",site,".txt"),header=T)
  TempStab <- select(TS_prod_tot,TS,site,order,simul)
  write.table( spread(TempStab,order,TS) , paste0("data/processed/TS_productivity_tot_",site,".txt"),sep="\t",row.names = F)
  }

# Add median confidence interval to temporal stability data ####
for (site in SITE){
  table <- read.table(paste0("data/processed/TS_productivity_tot_",site,".txt"),header=T)
  table3 <- median_conf_int(table)
  
  write.table(table3,paste0("data/processed/TS_productivity_tot_",site,"_with interval_median.txt"))
}

