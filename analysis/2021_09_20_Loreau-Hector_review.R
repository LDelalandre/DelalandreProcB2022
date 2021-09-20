source("R/Common variables.R")
library(tidyverse)


species <- read.table("data/processed/correspondence_SName_Id.txt",header=T)

MONOCULTURES <- read.table("data/processed/productivity_monoculture_ALL sites.txt",header=T) %>% 
  mutate(persists_mono = persists) %>% 
  select(-c(persists))

LH_all <- NULL
for (sit in SITE){
  MONO_site <- MONOCULTURES %>% 
    filter(site == sit)
  
  MIXTURES <- read.table(paste0("data/processed/productivity_specific_",sit,".txt"),header=T) %>% 
    mutate(persists_mixt = persists,
           SName = species) %>% 
    select(-c(persists,species))
  
  
  
  merged <- merge(MIXTURES,MONO_site, by = "SName" ) %>% 
    group_by(order,simul) %>% 
    rename(site = site.x)
  
  merged_2 <- merged %>% 
    rename(YOi=mixture_t_ha,Mi=monoculture) %>% 
    filter(Mi>0) %>% 
    group_by(site,order,simul) %>% 
    mutate(nb_sp_realized_pool=sum(persists_mixt)) %>% # nb of species in a community (which is a 2000-year-simul defined by a given regional pool of species)
    mutate(nb_sp_regional_pool=31-simul)# For simul 30, N=31-30=1 sp. For simul 1, N=31-1=30 sp, etc.
  
  
  LHinfo <- merged_2 %>% 
    mutate(YO=sum(YOi)) %>%
    mutate(RYEi=1/nb_sp_regional_pool) %>% # /!\ I'm not sure which pool I should take here! 
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
    
    mutate(Selection2=nb_sp_realized_pool*cov(DeltaRYi,Mi)) %>% 
    mutate(Relative_DeltaY = DeltaY/YE) %>% #/YO
    mutate(Relative_Selection = Selection/YE) %>% #/YO
    mutate(Relative_Complementarity = Cpltarity/YE) %>%  #/YO
    group_by(site,order,simul) %>%
    summarize(DeltaY=mean(DeltaY),Cpltarity=mean(Cpltarity),Selection=mean(Selection))
  
  LH_all <- rbind(LH_all, LHinfo)
}
write.csv(LH_all,"data/processed/Loreau-Hector coefficients.csv",row.names = F)


#_______________________________________________________________________________
# PLots ####
LH_all <- read.csv("data/processed/Loreau-Hector coefficients.csv")


# 15 species either Di or not ####
group_of_site <- function(site){
  if (site %in% c("GrandeDixence","Bever","Davos")){
    "cold"
  } else if (site %in% c("Adelboden","Huttwil")){
    "warm-wet"
  } else if (site %in% c("Bern","Schaffhausen","Basel")){
    "warm"
  } else if (site %in% c("Schwerin","Cottbus","Sion")){
    "warm-dry"
  } else {
    "problem"
  }
}

# how many species 
nbspecies <- 15

for (nbspecies in c(1:20)){
  # the most common
  commonx <- LH_all %>% 
    filter(order == "decreasing" & simul == nbspecies+1) %>% 
    mutate(categ_site = map_chr(site,group_of_site))
  # the most distinct
  distinctx <- LH_all %>% 
    filter(order == "increasing" & simul == nbspecies+1) %>% 
    mutate(categ_site = map_chr(site,group_of_site))
  
  grouped_commonx <- commonx %>% 
    ungroup() %>% 
    group_by(categ_site) %>%
    summarize_at(.vars = vars("DeltaY","Cpltarity","Selection"), .funs = ~ mean(.),na.rm = T)
  
  grouped_distinctx <- distinctx%>% 
    ungroup() %>% 
    group_by(categ_site) %>%
    summarize_at(.vars = vars("DeltaY","Cpltarity","Selection"), .funs = ~ mean(.),na.rm = T)
  
  
  
  # grouped_commonx %>%
  #   gather(DeltaY:Selection,key="effect",value="value") %>%
  #   filter(effect != "DeltaY") %>%
  #   ggplot(aes(x=categ_site,y=value,fill=effect)) +
  #   geom_col(width=0.6,position = "stack") + # or "dodge"
  #   scale_fill_grey(start = 0.6, end = 0.3)
  
  
  a <- grouped_commonx %>% 
    gather(DeltaY:Selection,key="effect",value="value") %>% 
    mutate(SpPool = "Common")
  b <- grouped_distinctx %>% 
    gather(DeltaY:Selection,key="effect",value="value") %>% 
    mutate(SpPool = "Distinct")
  
  plot <- rbind(a,b) %>% 
    mutate(effect = factor(effect, levels = c("DeltaY","Cpltarity","Selection"))) %>% 
    mutate(categ_site = factor( categ_site, levels = c("cold","warm-wet","warm","warm-dry"))) %>% 
    filter(effect!="DeltaY") %>% 
    ggplot(aes(x=categ_site,y=value,fill=effect)) +
    geom_col(width=0.6,position = "dodge") +
    scale_fill_brewer(palette="Dark2")+
    # scale_fill_grey(start = 0.7, end = 0.2)+
    facet_wrap(~SpPool) +
    theme(axis.text.x=element_text(angle = 35, vjust = 0.5))
  
  ggsave(paste0("figures/2021_09_Loreau-Hector/Sel_Cpltarity_groups of sites_",nbspecies," species.png"),plot)
  
  
}



