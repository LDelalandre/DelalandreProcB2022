library(tidyverse)

# mean biomass in time

MEAN <- NULL
for (site in SITE){
  mean <- read.table(paste0("data/raw/Output_ForCEEPS/",site,"/output-cmd2_",site,"_decreasing.txt/forceps.",site,".site_1_mean.txt"))
  colnames(mean) <- colnames_mean
  mean2 <- mean %>%
    select(date,nTrees..ha.,totalBiomass.t.ha.) %>% 
    mutate(date = date-1950)
  mean2$site <- site
  MEAN <- rbind(MEAN,mean2)
}

MEAN2 <- MEAN %>% 
  mutate(site=factor(site, levels = SITE))
ggplot(MEAN2,aes(x=date,y=totalBiomass.t.ha.)) +
  facet_wrap(~site)+
  geom_point()
