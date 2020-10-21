# Libraries ####
library("dplyr")
library("magrittr")
library(FactoMineR)
library("funrar")
library(ggplot2)
library(tidyr)

TOTAL <- read.table("data/processed/specific_biom_prod_complete.txt",header=T)

# Variables ####
colnames_mean<-colnames(read.table("data/colnames_mean.txt",header=T)) # idem
colnames_res<-colnames(read.table("data/colnames_res.txt", header=T))
colnames_prod <- colnames(read.table("data/colnames_productivityScene.txt",header=T))

Nbpatches <- 50
length <- 2000
yearstobejumped <- 999
timestep <- 100


SITE <- c("GrandeDixence","Bever","Davos","Adelboden","Huttwil","Schwerin","Bern","Cottbus","Basel","Schaffhausen"
          ,"Sion")
# ordered in growing temperature :
Site_descr <- read.table("data/Site description.txt",header=T)
ord_temp <- as.character(Site_descr[order(Site_descr$Temp_moy),]$Site)

ord_plots <- c("GrandeDixence","Bever","Davos","Basel","Cottbus","Schaffhausen","Schwerin","Sion",
               "Adelboden","Bern","Huttwil")


# SITE <- c("Bern","Bever","Cottbus","Huttwil","Adelboden","Basel","GrandeDixence")

ORDER <- c("increasing","decreasing","random_1","random_2","random_3","random_4",
           "random_5","random_6","random_7","random_8","random_9","random_10",
           "random_11","random_12","random_13","random_14",
           "random_15","random_16","random_17","random_18","random_19","random_20",
           "random_21","random_22","random_23","random_24",
           "random_25","random_26","random_27","random_28","random_29","random_30")

# MEASURE <- c("biomass_tot","sd_biomass_tot","CV_biomass_tot",
             # "productivity_tot","sd_productivity_tot","CV_productivity_tot")
MEASURE <- c("biomass_tot","productivity_tot","TS_productivity_tot")
