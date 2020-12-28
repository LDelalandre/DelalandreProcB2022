source("R/Common variables.R")
source("R/Analysis_data.R")
source("R/Monocultures_functions.R")

# CAREFUL not to run all the code: it includes processing data and writing tables, which:
# - takes a lot of time
# - can overwrite previous data frames

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

