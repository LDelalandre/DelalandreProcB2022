source("R/Analysis_data.R")
source("R/Common variables.R")
# CAREFUL not to run all the code: it includes processuing data and writing tables, which:
# - takes some time
# - can overwrite previous data frames

# Make the treatment of the measures automatic (duces the amount of code)
# Biomass removal experiments ####
for (site in SITE){ # Write a table with the specific final biomasses of all sets of simulations for each site
  # NB: Takes up to a few minuts to run!
  BIOMASSES <- biomass_specific(site=site) 
  write.table(BIOMASSES,paste0("data/processed/biomass_specific_",site,".txt"),sep="\t",row.names=F)
}

measure <- MEASURE[1] # biomass
for (site in SITE){
  # Total biomass of the simulation
  biomass_per_order <- biomass_tot(site)
  write.table(biomass_per_order,paste0("data/processed/biomass_tot_",site,".txt"),sep="\t",row.names=F)
  # Add a confidence interval to it
  biomass_per_order2 <- confidence_interval(site,measure = MEASURE[1])
  write.table(biomass_per_order2,paste0("data/processed/",measure,"_",site,"_with interval.txt"),sep="\t",row.names=F)
}

# plot the final biomasses in each condition
# more precisely, it is the averaged biomass on 10 evenly-spaced points in time during the last 1000 years.
for (site in SITE){
  result <- read.table(paste0("data/processed/total_biomass_final_",site,"_with interval.txt"),header=T)
  ggplot(result,aes(x=simul-1,y=decreasing,color="Removing distinct species first")) +
    labs(x="Number of species removed",y="Total biomass (t/ha)") +
    geom_line()+
    geom_ribbon(aes(ymin=int_min, ymax=int_max),fill="grey60", alpha=0.5,colour="black") +
    geom_line(aes(x=simul-1,y=increasing, color="Removing distinct species last")) +
    theme(legend.position = "bottom") +
    ggtitle(site) +
    ggsave(paste0("figures/Biomass_",site,".png"))
}


# sd(biomass) removal experiments ####
for (site in SITE){ # specific sd(biomass)
  sd_biom_sp <- sd_biomass_specific(site)
  write.table(sd_biom_sp, paste0("data/processed/sd_biomass_specific_",site,".txt"),sep="\t")
}

for (site in SITE){
  sd_biom_tot <- sd_biomass_tot(site)
  # CV of total biomass
  CV <- select(sd_biom_tot,CV,site,order,number)
  write.table( spread(CV,order,CV) , paste0("data/processed/CV_biomass_tot_",site,".txt"),sep="\t",row.names = F)
  # sd of total biomass
  sigma <- select(sd_biom_tot,sd_biomass,site,order,number)
  write.table( spread(sigma,order,sd_biomass) , paste0("data/processed/sd_biomass_tot_",site,".txt"),sep="\t",row.names = F)
}

for (site in SITE){ # Add confidence intervals...
  # ... to sd(biomass)
  measure <- MEASURE[2] # sd_biomass
  biomass_per_order2 <- confidence_interval(site,measure = MEASURE[1])
  write.table(biomass_per_order2,paste0("data/processed/",measure,"_",site,"_with interval.txt"),sep="\t",row.names=F)
  
  # ... to CV(biomass)
  measure <- MEASURE[3] # CV_biomass
  biomass_per_order2 <- confidence_interval(site,measure = MEASURE[1])
  write.table(biomass_per_order2,paste0("data/processed/",measure,"_",site,"_with interval.txt"),sep="\t",row.names=F)
}


# plot sd(biomass)
for (site in SITE){
  result <- read.table(paste0("data/processed/sd_biomass_species removal experiments_",site,".txt"))
  ggplot(result,aes(x=nb_removed,y=sd_biomass_dec,color="Removing distinct species first")) +
    labs(x="Number of species removed",y="sd(Total biomass)") +
    geom_line()+
    geom_ribbon(aes(ymin=int_min, ymax=int_max),fill="grey60", alpha=0.5,colour="black") +
    geom_line(aes(x=nb_removed,y=sd_biomass_inc, color="Removing distinct species last")) +
    theme(legend.position = "bottom") +
    ggtitle(site) +
    ggsave(paste0("figures/sd_biomass_",site,".png"))
}


# Productivity removal experiment ####

for (site in SITE){
  PROD <- productivity_specific(site)
  write.table(PROD,paste0("data/processed/productivity_specific_",site,".txt"),sep="\t",row.names=F)
  
  PROD_tot <- productivity_total(site)
  write.table(PROD_tot,paste0("data/processed/productivity_tot_",site,".txt"),sep="\t",row.names=F)
}


# plot productivity
for (site in SITE2){
  result <- read.table(paste0("data/processed/Productivity_species removal experiments_",site,".txt"))
  ggplot(result,aes(x=nb_removed,y=prod_dec,color="Removing distinct species first")) +
    labs(x="Number of species removed",y="Productivity") +
    geom_line()+
    geom_ribbon(aes(ymin=int_min, ymax=int_max),fill="grey60", alpha=0.5,colour="black") +
    geom_line(aes(x=nb_removed,y=prod_inc, color="Removing distinct species last")) +
    theme(legend.position = "bottom") +
    ggtitle(site) +
    ggsave(paste0("figures/Productivity_",site,".png"))
}

# sd(productivity) removal experiments ####
