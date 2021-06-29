library(tidyverse)

LH_ALL <- read.table("data/processed/LH_productivity_specific_every condition.txt",header=T)

per_mono <- LH_ALL %>% 
  filter(simul==1&order=="random_3") %>% # With 30 species
  group_by(site) %>% 
  filter(persists_mono==T) %>% 
  summarize(RS_mono=n())

per_mixt <- LH_ALL %>% 
  filter(simul==1&order=="random_3") %>% # With 30 species
  group_by(site) %>% 
  filter(persists_mixt==T) %>% 
  summarize(RS_mixt=n())

richness <- cbind(per_mono,per_mixt[,2])
ggplot(richness,aes(x=site,y=RS_mono))+
  geom_histogram(stat="identity",alpha=0.2) +
  theme(axis.text.x = element_text(angle = 60)) +
  geom_point( aes(x=site,y=RS_mixt),stat="identity")+
  ylab("Richness") +
  ggsave("paper_2/Richness.png")

