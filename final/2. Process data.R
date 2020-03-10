source("R/Analysis_data.R")
library(ggplot2)

SITE <- c("Bern","Bever","Cottbus","Huttwil")

# Biomass and sd(biomass) removal experiments ####
for (site in SITE){ # Write a table with the final biomasses of all sets of simulations for each site
  # NB: Takes up to a few minuts to run!
  removal_exp_final_biomasses(site=site) 
}


# plot the final biomasses in each condition
# more precisely, it is the averaged biomass on 10 evenly-spaced points in time during the last 1000 years.
for (site in SITE){
  result <- read.table(paste0("data/processed/Biomass_species removal experiments_",site,".txt"))
  ggplot(result,aes(x=nb_removed,y=biomass_dec,color="Removing distinct species first")) +
    labs(x="Number of species removed",y="Total biomass (T/ha)") +
    geom_line()+
    geom_ribbon(aes(ymin=int_min, ymax=int_max),fill="grey60", alpha=0.5,colour="black") +
    geom_line(aes(x=nb_removed,y=biomass_inc, color="Removing distinct species last")) +
    theme(legend.position = "bottom") +
    ggtitle(site) +
    ggsave(paste0("figures/Biomass_",site,".png"))
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
#NB: SO FAR, I HAVEN'T RUN THE SIMUL with productivity output for Bern and Bever
#I thus use, instead of SITE:
SITE2 <- c("Huttwil","Cottbus")
for (site in SITE2){ # Write a table with the productivity of all sets of simulations for each site
  removal_exp_productivity(site=site) 
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
