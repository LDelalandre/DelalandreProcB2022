source("final/0. Packages.R")
source("R/Analysis_data.R")
source("R/Common variables.R")


# I) Loreau-hector's partitioning of biodiversity effects ####
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
    # I divide by the regional pool, which represents the grains "seeded" (See Loreau & Hector, 2001, Nature): 
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
    # NB: N is the number elements on which we sum (cf. Loreau & Hector, 2001), i.e. "the number of species seeded or planted"
    
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



# II) Identification of key species for productivity, selection, and complementarity ####

# Script adapted from Eva Maire (eva.maire@umontpellier.fr).                                  
# Analyses were used in Maire, E., Villéger, S., Graham, N., Hoey, A., Cinner, J.,    
# Ferse, S., Aliaume, C, Booth, D., Feary, D., Kulbicki, M., Sandin, S., Vigliola, L. & Mouillot, D. (2018).    
# Community-wide scan identifies fish species associated to coral reef services globally. Proceedings  B. 
                                                                               
# Step 1: The goal is to model a given ecosystem service (Y) according to species richness (R); 
# to check its relevance according to its explanatory power and to save its 
# Akaike Information Criterion (AIC_M0) as a reference for the next step     
#                                                                                                              
# Step2: The goal is to identify species key for the studied ecosystem service (Y) adding each candidate       
# species (presence-absence) as an additional explanatory variable to M0 to compute model M1 and its           
# associated AIC (AIC_M1). Finally, a species is declared as a key potential contributor to the ecosystem      
# service if ΔAIC (AIC_M0-AIC_M1) > 4 and if its partial effect is positive (positive coefficient in the model)
source("R/df_maire.R")
source("R/infer_key_species.R")

# II.1) Generate data ####
LH_all_per_sp <- read.csv("data/processed/Loreau-Hector coefficients_per species.csv") %>% 
  mutate(site = factor(site,levels=SITE)) %>% 
  filter(persists_mixt==T)

for (site in SITE){
  for (nb_sp in c(10:20) ){
    simul_to_keep <- 30 - nb_sp +1
    dat <- df_maire(LH_all_per_sp,site,simul_to_keep)
    write.table(dat,paste0("data/processed/maire/occurrence_ppty",site,"_",30-simul_to_keep+1,"species.txt"),sep="\t",row.names=F)
  }
}


# II.2) Inference ####
PROPERTIES <- c("productivity","DeltaY","Selection","Complementarity")
property <- "Complementarity"
site <- SITE[5]
nb_sp <- 15

SUMMARY_KEY_SPECIES <- NULL
for (site in SITE){
  for (property in PROPERTIES){
    OCCURRENCE_PPTY <- NULL
    for (nb_sp in c(10:20) ){
      occurrence_ppty <- read.table(paste0("data/processed/maire/occurrence_ppty",site,"_",nb_sp,"species.txt"),header=T)
      OCCURRENCE_PPTY <- rbind(OCCURRENCE_PPTY,occurrence_ppty)
    }
    OP <- OCCURRENCE_PPTY %>% 
      mutate(richness = 30 - simul + 1)
    
    Summary_key_species <- infer_key_species(OP) %>% 
      rownames_to_column("SName")
    if (dim(Summary_key_species)[1] > 0){
      Summary_key_species$site <- site
      Summary_key_species$property <- property
    }
    SUMMARY_KEY_SPECIES <- rbind(SUMMARY_KEY_SPECIES,Summary_key_species)
  }
}


# II.3) Plot nb and ppty key species ####
species <- read.table("data/raw/distinctiveness of the species.txt",header=T) %>% 
  pull(SName) 
distinct_sp <- species[1:10]
common_sp <- species[11:30]

SUMM2 <- SUMMARY_KEY_SPECIES %>%
  mutate(property = if_else(property=="productivity","Productivity",property)) %>% 
  mutate(status = if_else(SName %in% distinct_sp,"Distinct","Common")) %>% 
  arrange(factor(property, levels = c("Productivity","DeltaY","Selection","Complementarity"))) %>% 
  mutate(property = factor(property,levels=c("Productivity","DeltaY","Selection","Complementarity"))) %>% 
  arrange(factor(site, levels = SITE)) %>% 
  mutate(site = factor(site,levels=SITE)) %>% 
  filter(!(property == "DeltaY")) %>% 
  group_by(site,property,status) %>% 
  mutate(count = n()) %>% 
  mutate(status = if_else(status=="Common","Ordinary",status))

cols <- c("Distinct" = "#F8766D", "Ordinary" = "#00BFC4")
SUMM3 <- SUMM2 %>% 
  group_by(site,property,status) %>% 
  mutate(count = n()) %>% 
  summarize(coeff = mean(coeff_sp),count=mean(count))

ggplot(SUMM3,aes(x=property,y=coeff,fill=status))+
  scale_color_manual(values = cols, aesthetics = "fill") +
  geom_histogram(stat="identity",position = "dodge") +
  theme(axis.text.x = element_text(angle = 66,hjust=1)) +
  ylim(c(0,2))+
  facet_wrap(~site) +
  geom_text(aes(label=count), vjust=0,position = position_dodge(width = 1)) +
  labs(fill = "Category of species") +
  xlab("Property") +
  ylab("Mean effect of key species") +
  ggsave(paste0("figures_tables/mean_effect_key_sp.png"),height=7,width=7)

