library(tidyverse)
biomass_per_species <- read.table("data/processed/biomass_specific_ALL sites.txt",header=T)

biomass_per_species_random <- biomass_per_species %>% 
  filter(!(order %in% c("increasing","decreasing")))

write.csv2(biomass_per_species_random,"data/processed/biomass_random_pour_cyrille.csv",row.names=F)
