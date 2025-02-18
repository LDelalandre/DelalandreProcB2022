library(tidyverse)
library(FactoMineR)
source("R/comp_fct_dist.R")

# Environmental conditions
sites <- read.table("data/Site description.txt",header=T) %>% 
  column_to_rownames("Site") %>% 
  select(Temp_moy,Annual_ppt,Na.kg.ha.an.,ProdMax.T.ha.an.)

traits_dist(sites) %>% 
  arrange(desc(Di)) %>% 
  rownames()



PCA <- PCA(sites)
coord_sites <- as.data.frame(PCA$ind$coord)
# Dim 1 : Productivity and Nutr availability
# Dim 2 : Temperature
# Precipitation on the two axis (possitive correlation with axis 1, negative with axis 2).