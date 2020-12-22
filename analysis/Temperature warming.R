source("R/Common variables.R")
source("R/Analysis_data.R")
source("R/Monocultures_functions.R")


site <- "Davos"
order <- "3degrees_increasing"

# Biomass in mixture
biom_per_species <- biomass_specific(site=site,order=order)

# Biomass in monoculture
specific_val <- specific_biomasses(site)



# Pool the two
# ATTENTION à ce qui se passe ici, je crois que ça modifie mon fichier biom_per_species de départ...
biomass <- biomasses(specif.biomass = biom_per_species,specific_val = specific_val)

BIOMASSES_sp <- biomass
BIOMASSES_sp <- rename(BIOMASSES_sp,mixture.t.ha.='mixture(t/ha)' )
aggregate(as.numeric(BIOMASSES_sp$mixture.t.ha.),list(site = BIOMASSES_sp$site,order = BIOMASSES_sp$order,simul = BIOMASSES_sp$simul),FUN=sum)# sum of the biomasses of the species


  BIOMASSES_sp <- read.table(paste0("data/processed/biomass_specific_",site,".txt"),header=T)
  BIOMASSES_tot <- aggregate(BIOMASSES_sp$mixture.t.ha.,list(site = BIOMASSES_sp$site,order = BIOMASSES_sp$order,simul = BIOMASSES_sp$simul),FUN=sum)# sum of the biomasses of the species
  
  biomass_per_order <- spread(BIOMASSES_tot,order,x) # data frame with each order in column
  biomass_per_order
