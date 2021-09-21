source("R/Common variables.R")
library(tidyverse)

species <- read.table("data/processed/correspondence_SName_Id.txt",header=T)

MONOCULTURES <- read.table("data/processed/productivity_monoculture_ALL sites.txt",header=T) %>% 
  mutate(persists_mono = persists) %>% 
  select(-c(persists))

LH_all <- NULL
sit <- "Bern"
  
MONO_site <- MONOCULTURES %>% 
  filter(site == sit)

MIXTURES <- read.table(paste0("data/processed/productivity_specific_",sit,".txt"),header=T) %>% 
  mutate(persists_mixt = persists,
         SName = species) %>% 
  select(-c(persists,species))

merged <- merge(MIXTURES,MONO_site, by = "SName" ) %>% 
  group_by(order,simul) %>% 
  rename(site = site.x)

#__________________________________________________________________________________________________
merged_filtered <- merged %>% 
  filter(order=="random_10"&simul==16)


merged_2 <- merged_filtered[1:3,] %>% 
  rename(YOi=mixture_t_ha,Mi=monoculture) %>% 
  filter(Mi>0) %>% 
  group_by(site,order,simul) %>% 
  mutate(nb_sp_realized_pool=sum(persists_mixt)) %>% # nb of species in a community (which is a 2000-year-simul defined by a given regional pool of species)
  mutate(nb_sp_regional_pool=n()) # For simul 30, N=31-30=1 sp. For simul 1, N=31-1=30 sp, etc.


LHinfo <- merged_2 %>% 
  mutate(YO=sum(YOi)) %>%
  mutate(RYEi=1/nb_sp_regional_pool) %>% 
  # I divide by the regional pool, which represents the grains "seeed" (See Loreau & Hector, 2001, Nature): 
  # "RYEi = expected relative yield of species i in the mixture, which is simply its proportion seeded or planted"  
  mutate(RYOi=YOi/Mi) %>%
  mutate(YEi=RYEi*Mi) %>%
  mutate(YE=sum(YEi)) %>%
  mutate(DeltaY=YO-YE) %>%
  
  mutate(DeltaRYi=RYOi-RYEi) %>%
  mutate(Mavg = mean(Mi)) %>%
  mutate(DeltaRYavg = mean(DeltaRYi)) %>%
  
  mutate(Cpltarity = nb_sp_regional_pool*DeltaRYavg*Mavg) %>%
  mutate(Selection = DeltaY - Cpltarity)%>%
  mutate(Selection2=nb_sp_regional_pool*cov(DeltaRYi,Mi)) %>% 
  
  group_by(site,order,simul) %>%
  summarize(DeltaY=mean(DeltaY),Cpltarity=mean(Cpltarity),Selection=mean(Selection),Selection2=mean(Selection2))
