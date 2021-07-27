# Script to process raw data #

# NB: raw data weight ~ 120Go. They will be available as zipped files on Zenodo soon.
# Raw data curation (aggregate biomass and productivity from the individual level to species level) lasts ~2 to 3 hours on my computer.
# Intermediary data is available in the "processed" folder.
# If you want to avoid data curation step, jump directly to script 3. to 6. in the "final" folder. 

source("final/0. Packages.R")
source("R/Common variables.R")
source("R/Analysis_data.R")
source("R/Monocultures_functions.R")


if (file.exists(file.path("data/raw/Output_ForCEEPS"))==F){
  cat("If you want to continue data curation, please download raw data from Zenodo (cf. script \"1.2. Import raw data\") \n
    N.B. Raw data is heavy (120 Go) and can be long to download and unzip. \n
    An alternative option is to go on with the analysis from data aggregated at the species level. \n
    To do so, please go to scripts 3. to 6. in the folder \"final\".")
} else {
  
  # Write informative tables ####
  make_table_name_sname()
  
  # CAREFUL not to run all the code: it includes processing data and writing tables, which:
  # - takes a lot of time
  # - can overwrite previous data frames
  
  # MONOCULTURES ####
  
  # Biomass per species ####
  # Generates a unique data frame with site, order, simul, the biomass of each species, 
  # and whether it persists or not (above threshold or not).
  MONO <- c()
  for (sit in SITE){
    specific_val <- 
      specific_biomasses_mono(sit) %>% 
      mutate(site=sit)
    persistent_sp <- specific_val %>% 
      filter(monoculture_relative>threshold) %>%
      pull(SName)
    
    specific_val$persists <- purrr::map_lgl(specific_val$SName,persist,persistent_sp)
    MONO <- rbind(MONO,specific_val)
  }
  
  write.table(MONO,"data/processed/biomass_mono_ALL sites.txt",sep="\t",row.names=F)
  
  
  # Productivity per species ####
  
  # Filter species above biomass threshold
  ALL_mono <- 
    read.table("data/processed/biomass_mono_ALL sites.txt",header=T) %>% 
    filter(persists==T) 
  PROD_mono <- NULL
  
  for (sit in SITE){
    # Vector of species above biomass threshold
    persistent.sp <- ALL_mono %>% 
      filter(site==sit) %>% 
      pull(SName)
    
    # Compute productivity per species
    prod_mono <- specific_productivities_mono(sit)
    prod_mono$site <- sit
    
    # Add a column indicating if the given species persists
    prod_mono$persists <- purrr::map_lgl(prod_mono$SName,persist,persistent.sp)
    
    PROD_mono <- rbind(PROD_mono,prod_mono)
  }
  write.table(PROD_mono,paste0("data/processed/productivity_monoculture_ALL sites.txt"),row.names=F,sep="\t")
  
  
  # BIOMASS MIXTURES ####
  
  # NB: Approximately 2 hours to run the code in this BIOMASS section!!
  ALL <- NULL
  # Compute the biomass of each "persistent" species in each site, order of removal and simulation.
  # NB: To define "persistent": I apply a biomass threshold, under which species are considered as absent from the community.
  
  # The data frame ALL contains the list of species in each site, order and simulation,to be used for 
  # subsetting productivity data.
  
  for (site in SITE){
    for (order in ORDER){
      for (number in c(1:30)){
        res<-try(read.table(paste0("data/raw/Output_ForCEEPS/",site,"/output-cmd2_",site,"_",order,".txt/forceps.",site,".site_",number,"_complete.txt")),silent=T) 
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
            mutate(sum=sum(biom_tot)) %>% # sum of total biomass used to compare with the biomass of one species
            filter(biom_tot>threshold*sum) %>% 
            pull(speciesShortName)
          
          
          # Average species biomass on 10 years every 100 years, and indicate if species are in the realized pool (in which case persists = TRUE)
          res3 <- res2 %>% 
            group_by(site,order,simul,speciesShortName) %>% 
            summarise(mean_biom_tot=mean(biom_tot)) %>%  # average of biomass on 10 years every 100 years
            mutate(
              persists = case_when(
                speciesShortName %in% filtered ~ TRUE,
                !(speciesShortName %in% filtered) ~ FALSE
              )
            ) # Add a column "persists" indicating if species is above biomass threshold
          
        }
        ALL <- rbind(ALL,res3)
      }
    }
  }
  
  write.table(ALL,"data/processed/biomass_specific_ALL sites.txt",sep="\t",row.names=F)
  
  # PRODUCTIVITY MIXTURES ####
  # Specific productivity ####
  ALL <- read.table("data/processed/biomass_specific_ALL sites.txt",header=T) %>% 
    filter(persists==TRUE)
  # NB: biomass of PERSISTENT species in all sites
  
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
  PROD <- NULL
  for (sit in SITE){
    prod <- 
      read.table(paste0("data/processed/productivity_specific_",sit,".txt"),header=T) %>% 
      filter(persists==T) %>% # keep species above biomass threshold
      group_by(site,order,simul) %>% 
      summarize(productivity=sum(mixture_t_ha)) %>% # compute community-level productivity
      spread(order,productivity) # spread the data frame to plot it
    prod[is.na(prod)] <- 0
    
    prod2 <- median_conf_int(as.data.frame(prod))
    PROD <- rbind(PROD,prod2)
  }
  write.table(PROD,paste0("data/processed/productivity_tot_with interval_median.txt"),row.names=F)
  
}

