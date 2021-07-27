library(tidyverse)

prod_tot <- read.table("data/processed/productivity_tot_with interval_median.txt",header=T)

prod_tot %>% 
  group_by(site) %>% 
  summarize(AUC=sum(decreasing))
