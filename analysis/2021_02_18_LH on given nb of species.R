source("R/Common variables.R")
source("R/analysis_data.R")
library("tidyverse")

ALL_POOLED <- read.table("data/processed/LH_productivity_specific_every condition.txt",header=T)

ALL_POOLED2 <- ALL_POOLED %>% 
  # filter(site==sit) %>% 
  rename(YOi=mixture_t_ha,Mi=monoculture_t_ha) %>% 
  filter(Mi>0) %>% 
  # select(-c(mixture_relative,persists)) %>% 
  group_by(site,order,simul) %>% 
  mutate(nb_sp_realized_pool=n()) %>% # nb of species in a community (which is a 2000-year-simul defined by a given regional pool of species)
  mutate(nb_sp_regional_pool=31-simul)# For simul 30, N=31-30=1 sp. For simul 1, N=31-1=30 sp, etc.

LHtest <- ALL_POOLED2 %>% 
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
  mutate(Relative_Complementarity = Cpltarity/YE) #/YO

#_____________________________________________________
# Analysis ####
LHinfo <- LHtest %>% 
  group_by(site,order,simul) %>% 
  summarize(DeltaY=mean(DeltaY),Cpltarity=mean(Cpltarity),Selection=mean(Selection))

LHinfoRelative <- LHtest %>% 
  group_by(site,order,simul) %>% 
  summarize(DeltaY=mean(Relative_DeltaY),Cpltarity=mean(Relative_Complementarity),Selection=mean(Relative_Selection))

# 15 species either Di or not ####

group_of_site <- function(site){
  if (site %in% c("GrandeDixence","Bever","Davos")){
    "cold"
  } else if (site %in% c("Adelboden","Bern","Huttwil")){
    "warm-wet"
  } else {
    "warm-dry"
  }
}

# how many species 
nbspecies <- 10

# the most common
common15 <- LHinfo %>% 
  filter(order == "decreasing" & simul == nbspecies+1) %>% 
  mutate(categ_site = map_chr(site,group_of_site))
#♠ the most distinct
distinct15 <- LHinfo %>% 
  filter(order == "increasing" & simul == nbspecies+1) %>% 
  mutate(categ_site = map_chr(site,group_of_site))

grouped_common15 <- common15 %>% 
  ungroup() %>% 
  group_by(categ_site) %>% 
  summarize_at(c("DeltaY","Cpltarity","Selection"),mean,na.rm = T)

grouped_distinct15 <- distinct15%>% 
  ungroup() %>% 
  group_by(categ_site) %>% 
  summarize_at(c("DeltaY","Cpltarity","Selection"),mean,na.rm = T)



# grouped_common15 %>% 
#   gather(DeltaY:Selection,key="effect",value="value") %>% 
#   filter(effect != "DeltaY") %>% 
#   ggplot(aes(x=categ_site,y=value,fill=effect)) +
#   geom_col(width=0.6,position = "stack") + # or "dodge"
#   scale_fill_grey(start = 0.6, end = 0.3)


a <- grouped_common15 %>% 
  gather(DeltaY:Selection,key="effect",value="value") %>% 
  mutate(SpPool = "Common")
b <- grouped_distinct15 %>% 
  gather(DeltaY:Selection,key="effect",value="value") %>% 
  mutate(SpPool = "Distinct")

plot <- rbind(a,b) %>% 
  mutate(effect = factor(effect, levels = c("DeltaY","Cpltarity","Selection"))) %>% 
  filter(effect!="DeltaY") %>% 
  ggplot(aes(x=categ_site,y=value,fill=effect)) +
  geom_col(width=0.6,position = "dodge") +
  scale_fill_brewer(palette="Dark2")+
  # scale_fill_grey(start = 0.7, end = 0.2)+
  facet_wrap(~SpPool)

ggsave(plot,paste0("figures/2021_02_22_Loreau-Hector/Sel_Cpltarity_groups of sites_",nbspecies," species.png"))


# tests of differences between conditions 
sel_di <- with(distinct15,lm(Selection~categ_site))
cpl_di <- with(distinct15,lm(Cpltarity~categ_site))
sel_co <- with(common15,lm(Selection~categ_site))
cpl_co <- with(common15,lm(Cpltarity~categ_site))

plot(cpl_co)
anova(cpl_co)

# test of differences between groups of species
# selection
shapiro.test(distinct15$Selection) # not normally distributed
shapiro.test(common15$Selection)

wilcox.test(distinct15$Selection,common15$Selection)

# complementarity
shapiro.test(distinct15$Cpltarity) # not normally distributed
shapiro.test(common15$Cpltarity)

wilcox.test(distinct15$Cpltarity,common15$Cpltarity)

# test of the difference of cpltarity vs. selection effect
wilcox.test(distinct15$Cpltarity,distinct15$Selection)
wilcox.test(common15$Cpltarity,common15$Selection)


# LH manually ####
Davos15sp <- ALL_POOLED %>% 
  filter(order=="decreasing",simul==15,site=="Davos")
  
Davos15sp %>% 
  filter(monoculture_t_ha>0) %>% 
  mutate(RYOi=mixture_t_ha/monoculture_t_ha) %>% 
  mutate(RYEi=1/15) %>% 
  mutate(deltaRYi=RYOi-RYEi) %>% 
  select(species,monoculture_t_ha,deltaRYi,mixture_t_ha,RYOi) %>% 
  arrange(monoculture_t_ha) %>% 
  summarize(select=cov(deltaRYi,monoculture_t_ha))

Adelboden15sp <- ALL_POOLED %>% 
  filter(order=="decreasing",simul==15,site=="Adelboden")

Adelboden15sp %>% 
  filter(monoculture_t_ha>0) %>% 
  mutate(RYOi=mixture_t_ha/monoculture_t_ha) %>% 
  mutate(RYEi=1/15) %>% 
  mutate(deltaRYi=RYOi-RYEi) %>% 
  select(species,monoculture_t_ha,deltaRYi,mixture_t_ha,RYOi) %>% 
  arrange(monoculture_t_ha) 
  # summarize(select=cov(deltaRYi,monoculture_t_ha))



# x species from every order ####
x=15
subsetsp <- LHinfo %>% filter(!(order %in%c("decreasing","increasing")) & simul == x+1)


orders_remov <- read.table("data/Orders of removal of cmd files.txt")
dist <- read.table("data/raw/distinctiveness of the species.txt",header=T)
correspNameId <- read.table("data/correspondence_SName_Id.txt",header=T) %>% 
  mutate(Di=dist$Di)


get_mean_di <- function(rownb,nbspecies,desired_stat){
  sp <- as.numeric(orders_remov[rownb,c(1:nbspecies+1)]) # chose the nb of species. Are they distinct on average?
  spdi <- correspNameId %>% 
    filter(Id %in% sp) %>% 
    pull(Di)
  if(desired_stat == "mean"){
    mean(spdi)
  } else if (desired_stat == "median"){
    median(spdi)
  } else if (desired_stat == "sd"){
    sd(spdi)
  } else {
    "specify_good_stat"
  }
}

meanDi <- orders_remov %>%   
  rownames_to_column(var="rownb") %>% 
  mutate(meanDi = map_dbl(rownb,get_mean_di,x,"mean")) %>% 
  mutate(medianDi = map_dbl(rownb,get_mean_di,x,"median")) %>% 
  mutate(sdDi = map_dbl(rownb,get_mean_di,x,"sd")) %>% 
  select(V1,meanDi,medianDi,sdDi) %>% 
  rename(order=V1)

extract_meanDi <- function(orde,desired_stat){
  meanDi %>% 
    filter(order==orde) %>% 
    pull(paste0(desired_stat,"Di"))
}


subsetsp2 <- subsetsp %>% 
  mutate(meanDi = map_dbl(order,extract_meanDi,"mean")) %>% 
  mutate(medianDi = map_dbl(order,extract_meanDi,"median")) %>% 
  mutate(sdDi = map_dbl(order,extract_meanDi,"sd")) 


ggplot(subsetsp2,aes(x=sdDi,y=DeltaY))+
  geom_point() +
  facet_wrap(~site)

ggplot(subsetsp2,aes(x=sdDi,y=Cpltarity))+
  geom_point() +
  facet_wrap(~site)

ggplot(subsetsp2,aes(x=sdDi,y=Selection))+
  geom_point() +
  facet_wrap(~site)

hist(subsetsp2$DeltaY) # normalité


mod1 <- lm(DeltaY~sdDi + site,data=subsetsp2)
mod2 <- lm(Selection~sdDi + site,data=subsetsp2)
mod3 <- lm(Cpltarity~sdDi + site,data=subsetsp2)

library(car)
Anova(mod1)
Anova(mod2)
Anova(mod3)

plot(mod)
summary(mod)

mod <- lme4::lmer(DeltaY~meanDi + (meanDi|site),subsetsp2)
anova(mod,test="Chisq")

# all the simul, and nb of species as a fixed, continuous effect ####
LHinfo2 <- LHinfo %>% filter(!(order %in%c("decreasing","increasing")) & simul %in% c(7:23))

LHinfoRelative2 <- LHinfoRelative %>% filter(!(order %in%c("decreasing","increasing")) & simul %in% c(7:23))

LHinfo3 <- LHinfoRelative2 %>%  #LHinfo2 %>% 
  mutate(meanDi = map_dbl(order,extract_meanDi,"mean")) %>% 
  mutate(medianDi = map_dbl(order,extract_meanDi,"median")) %>% 
  mutate(sdDi = map_dbl(order,extract_meanDi,"sd")) %>% 
  mutate(nbsp = 31-simul)

hist(LHinfo2$DeltaY) # normalité

ggplot(LHinfo3,aes(x=sdDi,y=DeltaY))+
  geom_point() +
  geom_jitter() +
  geom_smooth(method="lm") +
  facet_wrap(~site)

ggplot(LHinfo3,aes(x=nbsp,y=DeltaY))+
  geom_point() +
  geom_jitter() +
  geom_smooth(method="lm") 
  facet_wrap(~site)


mod1 <- lm(DeltaY~sdDi + site + nbsp,data=LHinfo3)
mod2 <- lm(Selection~sdDi + site + nbsp,data=LHinfo3)
mod3 <- lm(Cpltarity~sdDi + site + nbsp,data=LHinfo3)

Anova(mod1)
summary(mod1)

Anova(mod2)
summary(mod2)

Anova(mod3)
summary(mod3)

LHinfo3 %>% 
  filter(sdDi<0.3) %>% 
  pull(order) %>% 
  unique()
#__________________________________________________________
# Separate into two halves (di or not)

commonhalf <- dist %>% 
  arrange(Di) %>% 
  filter(rownames(.)%in% c(1:15)) %>% 
  pull(SName)
  

# count the number of distinct
get_nb_di <- function(orde,nbspecies){
  # choose the nb of species. Are they distinct on average?
  sp <- orders_remov[,0:nbspecies+1] %>% 
    filter(V1 == orde)
  
  di_of_sp <- correspNameId %>% 
    filter(Id %in% sp) %>% 
    pull(SName)
  
  length(which(!(di_of_sp %in% commonhalf)))
}

subsetsp3 <- subsetsp %>% 
  mutate(nbDi = map_dbl(order,get_nb_di,5))


ggplot(subsetsp3,aes(x=nbDi,y=Cpltarity))+
  geom_point() +
  facet_wrap(~site)
