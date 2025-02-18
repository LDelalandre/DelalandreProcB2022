library(tidyverse)
library(priceTools)
# GitHub of priceTools: https://github.com/ctkremer/priceTools/blob/master/vignettes/priceTools-intro.Rmd


LH_all <- read.csv("data/processed/Loreau-Hector coefficients.csv") 
LH_sp <- read.csv("data/processed/Loreau-Hector coefficients_per species.csv")


sit <- "Bever"
orde <- "random_10"
simulX <- 18
simulY <- 14

dataX <- LH_sp %>% 
  filter(site==sit & order == orde & simul == simulX) %>% 
  select(SName,YOi)

dataY <- LH_sp %>% 
  filter(site==sit & order == orde & simul == simulY) %>% 
  select(SName,YOi)


comm <- data.setup(list(dataX,dataY))
head(comm)

  
pp <- price.part(comm)
pp[1:5]
