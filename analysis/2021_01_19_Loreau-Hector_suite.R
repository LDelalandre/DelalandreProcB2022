source("R/Common variables.R")
source("R/analysis_data.R")
library("dplyr")

ALL_POOLED <- read.table("data/processed/LH_productivity_specific_every condition.txt",header=T)

sit <- "Bern"

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
  mutate(Relative_DeltaY = DeltaY/YO)

# Generate columns with median confidence interval 
LHplot3 <- LHtest %>% 
  group_by(site,order,simul) %>% 
  summarize(Relative_DeltaY=mean(Relative_DeltaY)) %>% # compute community-level productivity
  spread(order,Relative_DeltaY) # spread the data frame to plot it
LHplot3[is.na(LHplot3)] <- 0

sit="Bever"
result <- as.data.frame(LHplot3) %>% 
  filter(site==sit) %>% 
  median_conf_int()

ggplot(result,aes(x=simul-1,y=decreasing)) +  
  gg_removal_productivity()+
  ggtitle(sit)


# plots
LHplot <- LHtest %>% 
  filter(order %in% c("decreasing","increasing"))

ggplot(LHplot,aes(x=simul,y=Relative_DeltaY,color=order))+
  geom_line()+
  facet_wrap(~site)+
  scale_color_manual(name="Species lost first",labels=c("Distinct", "Common"),values=c("#F8766D", "#00BFC4"))

sit <- "Bever"
LHplot2 <- LHtest %>% 
  filter(site==sit)
ggplot(LHplot2,aes(x=simul,y=Relative_DeltaY,color=order))+
  geom_line()+
  geom_point(data=filter(LHplot2,order=="decreasing"))+
  geom_point(data=filter(LHplot2,order=="increasing"))
  scale_color_manual(name="Species lost first",labels=c("Distinct", "Common"),values=c("#F8766D", "#00BFC4"))


#_______________________________________________________________________________
# PLOT as the other data ####
# NB: used the functions from"final/4.Main figur..."

# Environmental conditions
sites <- read.table("data/Site description.txt",header=T)

# Sites and coordinates used hereafter
COORD <- sites %>% 
  select(Site,Temp_moy,Annual_ppt) %>% 
  arrange(Temp_moy,Annual_ppt)


# datatoplot <- read.table(paste0("data/processed/",measure,"_with interval_median.txt"),header=T)
measure <- "Relative_DeltaY"
sit <- "Bern"

# generate graphs for each site
PLOT <- list()
i <- 0
for (sit in COORD$Site){
  i <- i+1
  result <- as.data.frame(LHplot3) %>% 
    filter(site==sit) %>% 
    median_conf_int()
  
  
  PLOT[[i]] <- 
    ggplot(result,aes(x=simul-1,y=decreasing)) +  
    gg_removal_productivity()+
    ggtitle(sit)
    # ggplot(result,aes(x=simul,y=Relative_DeltaY,color=order))+
    # gg_removal_productivity()
    # geom_line()+
    # scale_color_manual(name="Species lost first",labels=c("Distinct", "Common"),values=c("#F8766D", "#00BFC4"))+
    # ggtitle(sit)+
    # theme(legend.position = "none" ) 
}

# Add axes labels for one site (Bever chosen here)
result <- as.data.frame(LHplot3) %>% 
  filter(site=="Bever") %>% 
  median_conf_int()
PLOT[[2]] <- ggplot(result,aes(x=simul-1,y=decreasing)) + 
  gg_removal_productivity() +
  theme(axis.title=element_text()) +
  ggtitle("Bever") +
  ylab("Relative Delta Y")

# generate the final graph
environmental_plot(PLOT,COORD)

