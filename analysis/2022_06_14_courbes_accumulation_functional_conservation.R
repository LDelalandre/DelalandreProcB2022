library(tidyverse)

data_biom <- read.table("data/processed/biomass_specific_ALL sites.txt",header=T)

data_biom_aggregated <- data_biom  %>% 
  filter(persists==T) %>% 
  filter(order == "increasing") %>% 
  group_by(site,simul) %>% 
  summarize(tot_biom = sum(mean_biom_tot))

plot_biom <- ggplot(data_biom_aggregated,aes(x=30-simul,y=tot_biom))+
  geom_line() +
  facet_wrap(~site) +
  xlab("number of species") +
  ylab("cumulated biomass (t/ha)") +
  ggtitle("Cumulated biomass of mixtures \n (Species added from distinct to ordinary)")

ggsave("figures/accumulation_curve_biomass.png",plot_biom)



data_prod_increasing <- read.table("data/processed/productivity_tot_with interval_median.txt",header=T) 

plot <- ggplot(data_prod,aes(x=30-simul,y=increasing))+
  geom_line() +
  facet_wrap(~site) +
  xlab("number of species") +
  ylab("cumulated productivity (t/ha)") +
  ggtitle("Cumulated productivity of mixtures \n (Species added from distinct to ordinary)")

ggsave("figures/accumulation_curve_productivity.png",plot)

