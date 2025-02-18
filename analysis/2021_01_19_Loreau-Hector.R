source("R/Common variables.R")
library("dplyr")
# NOTE: I don't keep this analysis, I pool on MONO

# functions ####
get_mono_prod <- function(spec,MONO){ # extracts the productivity in monoculture of one species from the MONO data frame.
  position <- grep(spec,MONO$SName) # where is the species spec in the MONO data frame
  MONO$monoculture[position] # Get its productivity in monoculture
}

# Add productivity in monoculture to mixture data frames ####
MONO1 <- read.table("data/processed/productivity_monoculture_ALL sites.txt",header=T)
POOLED_MIXT <- NULL
for (sit in SITE){
  MONO <-  MONO1 %>% 
    filter(site==sit) %>% 
    arrange(SName) 
  
  MIXT <- read.table(paste0("data/processed/productivity_specific_",sit,".txt"),header=T) %>% 
    arrange(order,simul,species)%>% 
    filter(persists==T) %>%  # select the species that persist
    mutate(monoculture_t_ha=purrr::map_dbl(species,get_mono_prod,MONO)) %>%  # add a column with species productivity in monoculture
    filter(monoculture_t_ha>0) # because we divide by monoculture productivity in loreau-hector, so it can't be divided by zero
  POOLED_MIXT <- rbind(POOLED_MIXT,MIXT)
}
write.table(POOLED_MIXT,"data/processed/mixt_mono_for LH_productivity_specific.txt",row.names=F,sep="\t")

# LH ####
POOLED_MIXT2 <- POOLED_MIXT %>% 
  filter(site==sit) %>% 
  rename(YOi=mixture_t_ha,Mi=monoculture_t_ha) %>% 
  select(-c(mixture_relative,persists)) %>% 
  group_by(site,order,simul) %>% 
  mutate(nb_sp_realized_pool=n()) %>% # nb of species in a community (which is a 2000-year-simul defined by a given regional pool of species)
  mutate(nb_sp_regional_pool=31-simul)# For simul 30, N=31-30=1 sp. For simul 1, N=31-1=30 sp, etc.
  
eee <- POOLED_MIXT2 %>% 
  mutate(YO=sum(YOi)) %>%
  mutate(RYEi=1/nb_sp_realized_pool) %>% # /!\ I'm not sure which pool I should take here! 
  # "N = number of species in the mixture" (Loreau and Hector, 2001).
  # And: "RYEi = expected relative yield of species i in the mixture, which is simply its proportion seeded or planted"
  # So it should be 1/regional pool, I guess.
  mutate(RYOi=YOi/Mi) %>%
  mutate(YEi=RYEi*Mi) %>%
  mutate(YE=sum(YEi)) %>%
  mutate(DeltaY=YO-YE) %>%
  
  mutate(DeltaRYi=RYOi-RYEi) %>%
  mutate(Mavg = mean(Mi)) %>%
  mutate(DeltaRYavg = mean(DeltaRYi)) %>%
  
  mutate(Cpltarity = nb_sp_realized_pool*DeltaRYavg*Mavg) %>%
  mutate(Selection = DeltaY - Cpltarity)%>%
  
  mutate(Selection2=nb_sp_realized_pool*cov(DeltaRYi,Mi))


