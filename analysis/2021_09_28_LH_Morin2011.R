library(tidyverse)
source("R/Common variables.R")

LH_all_per_sp <- read.csv("data/processed/Loreau-Hector coefficients_per species.csv") %>%
  mutate(site = factor(site,levels=SITE))

LH_all <- read.csv("data/processed/Loreau-Hector coefficients.csv") %>% 
  mutate(site = factor(site,levels=SITE))

LH_all2 <- LH_all %>% 
  filter(simul==1)

ggplot(LH_all2,aes(x=site,y=DeltaY)) +
  geom_boxplot()


# Reproduce figures from Morin 2011



# Figure 1
LH_all_xav <- filter(LH_all, simul %in% c(16:30))
ggplot(LH_all_xav,aes(x=30-simul+1,y=YO)) +
  geom_point()




ggplot(LH_all,aes(x=30-simul+1,y=YO)) +
  geom_point()+
  geom_smooth(method="lm")

# Non, /!\ lui c'est la realized sp richness, et moi c'est celle du départ !
# 1) pb dans la richesse réalisée
LH_all_per_sp %>% 
  filter(!(nb_sp_regional_pool == nb_sp_regional_pool2))

# pour faire ça :
LH_for_morin <- LH_all_per_sp %>% 
  tibble() %>% 
  filter(persists_mixt == T,
         persists_mono == T) %>% 
  group_by(site,order,simul) %>% 
  mutate(realized_sr = n())

LH_rich_test <- LH_for_morin %>% 
  filter(!(realized_sr == nb_sp_realized_pool)) %>% 
  select(site,order,simul,SName,nb_sp_realized_pool,realized_sr)
  
# pb dans les persistantes

LH_all_per_sp %>% 
  filter(persists_mixt == T) %>% 
  # filter(simul==1) %>% 
  select(persists_mono) %>% 
  group_by(persists_mono) %>% 
  summarize(n=n())
# --> Il y a des espèces qui persistent en mélange mais pas en mono.

rarecases <- LH_all_per_sp %>% 
  filter(persists_mixt == T & persists_mono == F)
rarecases %>% 
  pull(site) %>% 
  unique()

LH_all_per_sp %>% 
  filter(site=="GrandeDixence"&order=="random_21"&simul==19) %>% 
  filter(persists_mixt==T)


