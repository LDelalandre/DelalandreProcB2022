source("R/Common variables.R")
source("R/Analysis_data.R")

# I recode properly biomass computation (simpler and easier to track errors)

# Monocultures ####
MONO <- c()
for (sit in SITE){
  mono <- read.table(paste0("data/processed/biomass_monoculture_",sit,".txt"),header=T) %>% 
    filter(monoculture_relative>threshold) %>% 
    select(-abundance) %>% 
    mutate(site=sit)
  MONO <- rbind(MONO,mono)
}

write.table(MONO,"data/processed/biomass_mono_ALL sites.txt",sep="\t",row.names=F)


# Mixtures ####

ALL <- NULL
# Compute the biomass of each species in each site, order of removal and simulation.
# NB: I apply a biomass threshold, under which species are considered as absent from the community.

# The data frame ALL contains the list of species in each site, order and simulation,to be used for 
# subsetting productivity data.

# NB: Approximately 2 hours to run the code.

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


# Correlation biomass~distinctiveness ####

# in mono with threshold in biomass
MONOCULTURES <- NULL
for (sit in SITE){
  MONOsite <- 
    read.table("data/processed/biomass_mono_ALL sites.txt",header=T) %>% 
    filter(site==sit)
  
  DIST <- 
    read.table("data/raw/distinctiveness of the species.txt",header=T) %>% 
    filter(SName %in% MONOsite$SName) %>% 
    arrange(factor(SName,levels=MONOsite$SName))
  
  MONOsite$Di <- DIST$Di
  
  MONOCULTURES <- rbind(MONOCULTURES,MONOsite)
}

# in mono without threshold
MONO <- read.table("data/processed/biomass_mono_ALL sites.txt",header=T)

MONOCULTURES <- NULL
for (sit in SITE){
  MONOsite <- 
    MONO %>% 
    filter(site==sit) %>% 
    mutate(site=sit)
  
  DIST <- 
    read.table("data/raw/distinctiveness of the species.txt",header=T) %>% 
    filter(SName %in% MONOsite$SName) %>% 
    arrange(factor(SName,levels=MONOsite$SName))
  
  MONOsite$Di <- DIST$Di
  MONOCULTURES <- rbind(MONOCULTURES,MONOsite)
}


plotcor <- ggplot(MONOCULTURES,aes(x=Di,y=monoculture.t.ha.,label=SName))+
  geom_point()+
  facet_wrap(~site)+
  geom_smooth(method=lm)+
  ggpubr::stat_cor(method="spearman")+
  geom_label()

ggsave(filename = paste0("figures/correlation biomass_Di.png"), 
       plot = plotcor,
       width = 50, 
       height = 40,
       units = "cm",
       dpi = 300)


