source("R/Analysis_data.R")
source("R/Common variables.R")
library("functClust")
library(tidyverse)

options(verbose = TRUE) # So that I see the advance of the computations


# res <- fclust_read(filename = paste0("figures/functClust/functclust/",site,"_",property))

# fclust_plot(res, main = paste(site,property), 
            # opt.tree = list("prd"),
            # nbcl=nbcl,
            # # opt.perf = list("stats_II"))
            # opt.perf = list("prd", "aov", pvalue = 0.01, "pub"))


# fclust_plot(res, paste(site,property), 
            # opt.motif = list("obs",  "aov", pvalue = 0.05))
# Ne marche pas !! (ou très très très très long...)


# Read the files printed by fclust_write.

site <- "Huttwil"
property <- "Selection"
nbCl <- 2 # number of clusters of species

inputs <- read.csv(paste0("figures/functClust/functClust/",site,"_",property,".inputs.csv")) %>% 
  rename(Ass = names.fobs)
matrices <- read.csv(paste0("figures/functClust/functClust/",site,"_",property,".matrices.csv"))
trees <- read.csv(paste0("figures/functClust/functClust/",site,"_",property,".trees.csv"))

# affectation of assemblages to different assembly motives 
affectAss <- matrices %>% 
  filter(mat == "mMotifs") %>% 
  filter(nbClu == nbCl) %>%
  select(-c(mat,nbClu)) %>% 
  gather(key = "Ass", value = "Motif") %>% 
  mutate(Motif = as.factor(Motif))

# which species in which cluster of species
# NB I name the cluster with the smaller number of species cluster a, and the other b.
affectSp <-
  trees %>% 
  filter(nbClu == nbCl) %>% 
  select(-c(nbOpt,nbClu)) %>% 
  gather(key = "Species", value = "Motif") %>% 
  filter(grepl(".II",Species)) %>% 
  filter(!grepl("R2.II",Species)) %>%  # remove R2
  mutate(Species = str_replace_all(Species, ".II", "") ) %>%
  group_by(Motif) %>% 
  mutate(nb_sp_in_motif = n()) %>% 
  ungroup() %>% 
  mutate(more_sp = max(nb_sp_in_motif)) %>% 
  mutate(more_sp = if_else (more_sp == nb_sp_in_motif, T, F)) %>% 
  mutate(Motif = if_else(more_sp,"b", "a"))

# which cluster of species in which assemblage
give_motif <- function(Assemblage){
  # example: Assemblage <- "Huttwil_increasing_29"
  sp_in_assemblage <- inputs %>% 
    filter(Ass == Assemblage) %>% 
    select(-c("Ass","fobs","names.xpr","xpr")) %>% 
    gather(key = Species, value = presence) %>% 
    filter(presence == 1) %>% 
    pull(Species)
  
  cluster_in_assembage <- affectSp %>% 
    filter(Species %in% sp_in_assemblage) %>% 
    pull(Motif) %>% 
    unique()
  
  if ("a" %in% cluster_in_assembage & "b" %in% cluster_in_assembage){
    "ab"
  } else if ("a" %in% cluster_in_assembage & !("b" %in% cluster_in_assembage)){
    "a"
  } else if (!("a" %in% cluster_in_assembage) & "b" %in% cluster_in_assembage){
    "b"
  }
}


motif_of_assemblages <- affectAss %>% 
  mutate(Motif_old = Motif) %>% 
  mutate(Motif = map_chr(Ass,give_motif))
  
  
# Observed function per assembly motif
fobs_per_motif <- merge(inputs,motif_of_assemblages,by = "Ass")

ggplot(fobs_per_motif,aes(x=Motif,y=fobs))+
  geom_boxplot()
