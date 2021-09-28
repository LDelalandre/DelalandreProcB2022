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
  
  LH1 <- merged %>% 
    rename(YOi=mixture_t_ha,Mi=monoculture) %>% 
    filter(Mi>0) %>% 
    group_by(site,order,simul) %>% 
    mutate(nb_sp_realized_pool=sum(persists_mixt)) %>% # nb of species in a community (which is a 2000-year-simul defined by a given regional pool of species)
    mutate(nb_sp_regional_pool=31-simul) %>% # For simul 30, N=31-30=1 sp. For simul 1, N=31-1=30 sp, etc.
    mutate(nb_sp_local_pool = n()) # number of rows. cf. Chauvet 2017.
  
  LH2 <- LH1 %>% 
    mutate(YO=sum(YOi)) %>%
    mutate(RYEi=1/nb_sp_regional_pool) %>% 
    # I divide by the regional pool, which represents the grains "seeed" (See Loreau & Hector, 2001, Nature): 
    # "RYEi = expected relative yield of species i in the mixture, which is simply its proportion seeded or planted"  
    mutate(RYOi=YOi/Mi) %>%
    mutate(YEi=RYEi*Mi) %>%
    mutate(YE=sum(YEi)) %>%
    mutate(DeltaY=YO-YE) %>%
    mutate(DeltaRYi=RYOi-RYEi)
    
  LH3 <- LH2 %>% 
    mutate(Selection = nb_sp_regional_pool * ( mean(DeltaRYi * Mi) - mean(DeltaRYi) * mean(Mi) ),
              Complementarity = nb_sp_regional_pool * mean(DeltaRYi)*mean(Mi),
              DeltaY = nb_sp_regional_pool * mean(DeltaRYi * Mi)) %>%  # idem as YO - YE
    mutate(sum = Selection + Complementarity) %>% 
    # NB: N is the number elements on which we sum (cf. Loreau & Hector, 2001), i.e. the "number of species seeded or planted"
    
    group_by(site,order,simul)
  
  LH4 <- LH3 %>% 
    select(SName,site,order,simul,persists_mono,Mi,persists_mixt,YOi,YE,DeltaY,Selection,Complementarity)
  
  LH_all <- rbind(LH_all, LH3)
}

write.csv(LH_all,"data/processed/Loreau-Hector coefficients_per species.csv",row.names = F)

LH_all2 <- LH_all %>% 
  group_by(site,order,simul) %>% 
  summarize(YO = mean(YO),YE = mean(YE),DeltaY=mean(DeltaY),Complementarity=mean(Complementarity),Selection=mean(Selection),sum=mean(sum),)

write.csv(LH_all2,"data/processed/Loreau-Hector coefficients.csv",row.names = F)


#_______________________________________________________________________________
# PLots ####
LH_all <- read.csv("data/processed/Loreau-Hector coefficients.csv") %>% 
  mutate(site = factor(site,levels=SITE))


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

for (nbspecies in c(1:30)){
  # the most common
  commonx <- LH_all %>% 
    filter(order == "decreasing" & simul == 30 - nbspecies+1) %>% 
    mutate(categ_site = map_chr(site,group_of_site))
  # the most distinct
  distinctx <- LH_all %>% 
    filter(order == "increasing" & simul == 30 - nbspecies+1) %>% 
    mutate(categ_site = map_chr(site,group_of_site))
  
  grouped_commonx <- commonx %>% 
    ungroup() %>% 
    group_by(categ_site) %>%
    summarize_at(.vars = vars("DeltaY","Complementarity","Selection"), .funs = ~ mean(.),na.rm = T)
  
  grouped_distinctx <- distinctx%>% 
    ungroup() %>% 
    group_by(categ_site) %>%
    summarize_at(.vars = vars("DeltaY","Complementarity","Selection"), .funs = ~ mean(.),na.rm = T)
  
  
  
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
    mutate(effect = factor(effect, levels = c("DeltaY","Complementarity","Selection"))) %>% 
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



# Search possible analyses ####
sim = 11
metrics <- "DeltaY"
for (sim in c(30:1)){
  toplot <- LH_all %>% 
    filter(simul==sim & order %in% c("decreasing","increasing")) %>% 
    mutate(order = factor(order,levels=c("increasing","decreasing")))
  
  ggplot(toplot,aes_string(x="order",y=metrics, group = 1, colour = "order")) +
    geom_line(col = "black") +
    geom_point(size=3) +
    facet_wrap(~site) +
    scale_x_discrete(labels= c("Distinct lost first","Common lost first")) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank()) +
    scale_colour_discrete(labels = c('Distinct species', 'Common species'),
                          name="") +
    ggtitle(paste(metrics,"; number of species = ",30-sim+1)) +
    ggsave(paste0("figures/2021_09_21_Loreau-Hector/",metrics,"_N=",30-sim+1,".png"))
}



metrics <- "Selection"
toplot <- LH_all %>% 
  filter(simul %in% c(10:20)) %>% 
  filter(order %in% c("decreasing","increasing")) %>% 
  mutate(order = factor(order,levels=c("increasing","decreasing")))

ggplot(toplot,aes_string(x="order",y=metrics, group = 1, colour = "order")) +
  # geom_line(col = "black") +
  geom_point(size=3) +
  facet_wrap(~site) +
  scale_x_discrete(labels= c("Distinct lost first","Common lost first")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  scale_colour_discrete(labels = c('Distinct species', 'Common species'),
                        name="") +
  ggtitle(paste(metrics,"; number of species = ",30-sim+1))
