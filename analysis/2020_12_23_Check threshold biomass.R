source("R/Common variables.R")
source("R/Analysis_data.R")

MONO <- c()
for (sit in SITE){
  mono <- read.table(paste0("data/processed/biomass_monoculture_",sit,".txt"),header=T) %>% 
    filter(monoculture_relative>0.01)
  MONO <- c(MONO,dim(mono)[1])
}
data.frame(SITE,MONO)



MIXT <- c()
for (sit in SITE){
  mixt <- read.table(paste0("data/processed/biomass_specific_",sit,"_with monocultures.txt"),header=T) %>% 
    filter(order=="increasing"&simul=="1") %>% 
    filter(mixture_relative>0.0001)
    # filter(mixture.t.ha.>0.0001*sum(mixture.t.ha.))
  MIXT <- c(MIXT,dim(mixt)[1])
}
data.frame(SITE,MIXT)


#############################
# I recode properly biomass computation (simpler and easier to track errors)

# Mixtures ####
ALL <- NULL
# I want the biomass of each species in a simul with all the species
order <- "decreasing"
number <- 1

# for (site in SITE){
  # for (order in ORDER){
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
          mutate(site=site,order=order,simul=number) %>% 
          group_by(site,order,simul,speciesShortName) %>% 
          summarise(mean_biom_tot=mean(biom_tot)) # average of biomass on 10 years every 100 years
      }
      ALL <- rbind(ALL,res2)
    # }
  # }
}

# Measure the number of species above a threshold = Chauvet's realized species richness.
ALL %>% 
  filter(order==order,number==number) %>% # keep the first simul (=all species) of a given order
  mutate(sum_biom=sum(mean_biom_tot)) %>% 
  filter(mean_biom_tot > 0.01*sum_biom) %>% 
  summarise(n=n()) %>% 
  arrange(factor(site,levels=SITE))

ALL2 <- ALL %>% 
  filter(order==order,number==number) %>% # keep the first simul (=all species) of a given order
  mutate(sum_biom=sum(mean_biom_tot)) %>% 
  filter(mean_biom_tot > 0.01*sum_biom)

# Compare old and new computations
biomBeverOld <- biomass_specific("Bever","decreasing")

ALL

# Plus petite brique (un ordre seulement, avec une sous-fonction du vieux code)
res<-try(read.table(paste0("data/raw/output-cmd2_",site,"_",order,".txt/forceps.",site,".site_",number,"_complete.txt")),silent=T) 
colnames(res) <- colnames_res
int <- temporal_plot(res)

int2 <- int %>% 
  group_by(species,date) %>% 
  summarise(biom_tot=sum(biomass)) %>% 
  summarise(mean_biom_tot=mean(biom_tot))
int2$mean_biom_tot <- int2$mean_biom_tot/(1000*0.08*Nbpatches)

ALL2 <- ALL %>% 
  filter(site==site,order==order,simul==number) %>%  # keep the first simul (=all species) of a given order
  arrange(factor(speciesShortName,levels=int2$species))

cbind(int2,ALL2) # same columns

int3 <- temporal_plot_threshold(int)
int3$biomass <- int3$biomass/(1000*0.08*Nbpatches)
int4 <- int3 %>% 
  group_by(species,date) %>% 
  summarise(biom_tot=sum(biomass)) %>% 
  summarise(mean_biom_tot=mean(biom_tot))

# Regarder sur la dernière année seulement (comme Morin 2011)

