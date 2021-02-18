source("R/Before simulations.R")
source("R/Common variables.R")
source("R/comp_fct_dist.R")

traits<-read.table("data/Traits of the species_complete.txt",header=T)
selected_traits<-choice_traits_1(traits) # data.frame with the traits of the species

# I add the group of species established with the PCA : 0: Common; the rest : the 30% (9/30) most distinct.
# 1: top-right; 2: top-left; 3: bottom-left
group_of_distinct <- function(species){
  if (species %in% c("LDec","PCem","PMon","PAbi")){
    "top-left"
  } else if (species %in% c("BPen","AVir")){
    "bottom-left"
  } else if (species %in% c("AAlb","TBac","QRob")){
    "top-right"
  } else {
    "common"
  }
}



prod_per_site_mixt <- function(ord,sim,mono){
  # chose the order of removal and the simul, and it plots abundance of each species present
  # mono is a boolean: give biomass of productivity
  MIXT <- read.table("data/processed/LH_productivity_specific_every condition.txt",header=T)
  
  MIXTURES <- NULL
  for (sit in SITE){
    MIXTsite <- 
      MIXT %>% 
      filter(site==sit) %>% 
      mutate(site=sit) %>% 
      filter(simul==sim & order==ord) %>% 
      arrange(species)
    # filter(persists==T) # If I keep only the local pool (species above biomass threshold), I lose too much statistical power
    # consequently, I don't filter them 
    DIST <- 
      read.table("data/raw/distinctiveness of the species.txt",header=T) %>% 
      filter(SName %in% MIXTsite$species) %>% 
      arrange(SName)
    
    MIXTsite$Di <- DIST$Di
    MIXTURES <- rbind(MIXTURES,MIXTsite)
  }
  
  # plot ####
  MIXTURES2 <- MIXTURES %>% mutate(group = map_chr(species,group_of_distinct))
  
  if (mono){
    plotcor <- ggplot(MIXTURES2,aes(x=Di,y=monoculture_t_ha,label=species,color=group))+
      geom_point()+
      facet_wrap(~site,ncol=3)+
      # geom_smooth(method=lm)+
      # ggpubr::stat_cor(method="spearman")+
      geom_label()
    ggsave(filename = paste0("figures/2021_02_18_correlation prod di/mono_ord=",ord,"_sp_lost=",nb_sp_lost,".png"), 
           plot = plotcor,
           width = 35, 
           height = 40,
           units = "cm",
           dpi = 150) # normally dpi=300 to print it
    
    plotcor
  } else{
    plotcor <- ggplot(MIXTURES2,aes(x=Di,y=mixture_t_ha,label=species,color=group))+
      geom_point()+
      facet_wrap(~site,ncol=3)+
      # geom_smooth(method=lm)+
      # ggpubr::stat_cor(method="spearman")+
      geom_label() +
      ylim(0,1.5)
    ggsave(filename = paste0("figures/2021_02_18_correlation prod di/mixt_ord=",ord,"_sp_lost=",nb_sp_lost,".png"), 
           plot = plotcor,
           width = 35, 
           height = 40,
           units = "cm",
           dpi = 150) # normally dpi=300 to print it
    
    # plotcor
  }

}

ord <- "increasing"
nb_sp_lost <-0

for (nb_sp_lost in c(0:29)){
  prod_per_site_mixt(ord,nb_sp_lost+1,mono=F)
}




# save ####

