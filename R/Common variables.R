# Libraries ####
library("dplyr")
library("magrittr")
library(FactoMineR)
library("funrar")
library(ggplot2)
library(tidyr)

# Variables ####
colnames_mean<-colnames(read.table("data/colnames_mean.txt",header=T)) # idem
colnames_res<-colnames(read.table("data/colnames_res.txt", header=T))

Nbpatches <- 50
length <- 2000
yearstobejumped <- 999
timestep <- 100

SITE <- c("Bern","Bever","Cottbus","Huttwil")
ORDER <- c("increasing","decreasing","random_1","random_2","random_3","random_4",
           "random_5","random_6","random_7","random_8","random_9","random_10")