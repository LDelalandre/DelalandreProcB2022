source("R/Common variables.R")
library(tidyverse)

LH_all <- read.csv("data/processed/Loreau-Hector coefficients.csv") %>% 
  mutate(site = factor(site,levels=SITE))

species <- read.table("data/processed/correspondence_SName_Id.txt",header=T)
orders <- read.table("data/processed/Richness gradients orders.txt",header=T)


# Explore ####
sp15 <- orders %>% 
  filter(rownames(.)=="16") %>% # Species which will be removed when I have 30 - 16 + 1 = 15 species
  gather(key = order, value = Id15) 

sp15 %>% 
  mutate(Name15 = map(Id15, function(x) species %>% filter(Id==x) %>% pull(SName)  ))



