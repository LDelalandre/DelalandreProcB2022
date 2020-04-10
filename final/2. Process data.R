source("R/Analysis_data.R")
source("R/Common variables.R")
# CAREFUL not to run all the code: it includes processing data and writing tables, which:
# - takes some time
# - can overwrite previous data frames

# Biomass removal experiments ####
for (site in SITE){ # Write a table with the specific final biomasses of all sets of simulations for each site
  # NB: Takes a long time to run!!
  BIOMASSES <- biomass_specific(site=site) 
  write.table(BIOMASSES,paste0("data/processed/biomass_specific_",site,".txt"),sep="\t",row.names=F)
}

for (site in SITE){
  # Total biomass of the simulation
  biomass_per_order <- biomass_tot(site)
  write.table(biomass_per_order,paste0("data/processed/biomass_tot_",site,".txt"),sep="\t",row.names=F)
}


# sd(biomass) removal experiments ####
# specific sd(biomass)
for (site in SITE){ # NB: quite long to run!
  sd_biom_sp <- sd_biomass_specific(site)
  write.table(sd_biom_sp, paste0("data/processed/sd_biomass_specific_",site,".txt"),sep="\t")
}

# total sd(biomass)
for (site in SITE){ # Takes a looooong time to run!!
  sd_biom_tot <- sd_biomass_tot(site)
  write.table(sd_biom_tot, paste0("data/processed/sd and CV_biomass_tot_",site,".txt"),sep="\t",row.names = F)
}

for (site in SITE){ # Spread the datasets so that plotting is straightforward.
  sd_biom_tot <- read.table(paste0("data/processed/sd and CV_biomass_tot_",site,".txt"),header=T)
  # CV of total biomass
  CV <- select(sd_biom_tot,CV,site,order,simul)
  write.table( spread(CV,order,CV) , paste0("data/processed/CV_biomass_tot_",site,".txt"),sep="\t",row.names = F)
  # sd of total biomass
  sigma <- select(sd_biom_tot,sd,site,order,simul)
  write.table( spread(sigma,order,sd) , paste0("data/processed/sd_biomass_tot_",site,".txt"),sep="\t",row.names = F)
}


# Productivity removal experiment ####

for (site in SITE){
  PROD <- productivity_specific(site)
  write.table(PROD,paste0("data/processed/productivity_specific_",site,".txt"),sep="\t",row.names=F)
  
  PROD_tot <- productivity_total(site)
  write.table(PROD_tot,paste0("data/processed/productivity_tot_",site,".txt"),sep="\t",row.names=F)
}

# sd(productivity) removal experiments ####
for (site in SITE){
  sd_prod_sp <- sd_productivity_specific(site)
  write.table(sd_prod_sp, paste0("data/processed/sd_productivity_specific_",site,".txt"),sep="\t")
}

# total sd(productivity)
for (site in SITE){ 
  sd_prod_tot <- sd_productivity_tot(site)
  write.table(sd_prod_tot, paste0("data/processed/sd and CV_productivity_tot_",site,".txt"),sep="\t",row.names = F)
}

for (site in SITE){ 
  sd_prod_tot <- read.table(paste0("data/processed/sd and CV_productivity_tot_",site,".txt"),header=T)
  # CV 
  CV <- select(sd_prod_tot,CV,site,order,simul)
  write.table( spread(CV,order,CV) , paste0("data/processed/CV_productivity_tot_",site,".txt"),sep="\t",row.names = F)
  # sd 
  sigma <- select(sd_biom_tot,sd,site,order,simul)
  write.table( spread(sigma,order,sd) , paste0("data/processed/sd_productivity_tot_",site,".txt"),sep="\t",row.names = F)
}

# Add confidence intervals ####
for (site in SITE){ # Add confidence intervals...
  for (measure in MEASURE){
    conf <- confidence_interval(site,measure)
    write.table(conf,paste0("data/processed/",measure,"_",site,"_with interval.txt"),sep="\t",row.names=F)
  }
}


