source("R/Common variables.R")
source("R/Analysis_data.R")
source("R/Monocultures_functions.R")

persist <- function(sp,persistent.sp){
  # returns true if and only if species sp is present in the vector persisten.sp
  any(grepl(sp,persistent.sp))
}
# I should use this function in the BIOMASS section: instead of filter,
# I should add a column indicating if a species persists or not.
# Then, when using this file later, I should filter species which persist: filter (persists==T)

# CAREFUL not to run all the code: it includes processing data and writing tables, which:
# - takes a lot of time
# - can overwrite previous data frames


# BIOMASS ####

# NB: Approximately 2 hours to run the code in this BIOMASS section!!
ALL <- NULL
# Compute the biomass of each species in each site, order of removal and simulation.
# NB: I apply a biomass threshold, under which species are considered as absent from the community.

# The data frame ALL contains the list of species in each site, order and simulation,to be used for 
# subsetting productivity data.

for (site in SITE){
  for (order in ORDER){
    for (number in c(1:30)){
      res<-try(read.table(paste0("data/raw/output-cmd2_",site,"_",order,".txt/forceps.",site,".site_",number,"_complete.txt")),silent=T) 
      if (class(res) != "try-error"){# sometimes, the files are empty, and it returns an error message
        colnames(res) <- colnames_res
        # gives mean biomass (in t/ha) of each species on 10 years sampled every 100 years
        # (done for one species pool (simul = 1 : all the species, simul = 30 : 1 sp remaining), for one
        # order of removal, and in one site)
        res2 <- res %>% 
          group_by(date,speciesShortName) %>% 
          select(date,speciesShortName,biomass.kg.) %>% 
          summarise(biom_tot=sum(biomass.kg.)) %>% # sum of the biomass per species per date
          mutate(biom_tot=biom_tot/(1000*0.08*Nbpatches)) %>% # so that the unit becomes t/ha
          mutate(site=site,order=order,simul=number) 
        
        # Apply a biomass threshold and keep species whose biomass in the last year is above a given fraction (the variable "threshold")
        # of the total biomass of the community.
        # Filtered contains the filtered species:
        filtered <- res2 %>% 
          group_by(date) %>% 
          filter(date==3950) %>% 
          mutate(sum=sum(biom_tot)) %>% 
          filter(biom_tot>threshold*sum)
        
        # Filter species in res2, and then average their biomass on 10 years every 100 years.
        res3 <- res2 %>% 
          filter(speciesShortName %in% filtered$speciesShortName) %>% 
          group_by(site,order,simul,speciesShortName) %>% 
          summarise(mean_biom_tot=mean(biom_tot)) # average of biomass on 10 years every 100 years
      }
      ALL <- rbind(ALL,res3)
    }
  }
}

# write.table(ALL,"data/processed/biomass_specific_ALL sites.txt",sep="\t",row.names=F)

# PRODUCTIVITY ####
# Specific productivity ####
ALL <- read.table("data/processed/biomass_specific_ALL sites.txt",header=T)

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
  write.table(prod2,paste0("data/processed/productivity_tot_",sit,"_with interval_median.txt"),row.names=F)
}

# MONOCULTURES ####

# Biomass per species ####
# Generates a unique data frame with site, order, simul, the biomass of each species, 
# and whether it persists or not (above threshold or not).
MONO <- c()
for (sit in SITE){
  specific_val <- 
    specific_biomasses(sit) %>% 
    mutate(site=sit)
  persistent_sp <- specific_val %>% 
    filter(monoculture_relative>threshold) %>%
    pull(SName)
  
  specific_val$persists <- purrr::map_lgl(specific_val$SName,persist,persistent_sp)
  MONO <- rbind(MONO,specific_val)
}

write.table(MONO,"data/processed/biomass_mono_ALL sites.txt",sep="\t",row.names=F)


# Productivity per species ####
for (site in SITE){
  # 1) Compute specific productivity (and TS of productivity) per site in monoculture
  specific_val <- specific_productivities(site)
  write.table(specific_val,paste0("data/processed/productivity_monoculture_",site,".txt"),row.names=F,sep="\t")
  
  # 2) Add a column with productivity in monoculture to the files giving productivity in mixture
  # NB: For Loreau-Hector, it adds species that are absent, but were present in the regional pool 
  # (and attributes them a productivity value o zero)
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

