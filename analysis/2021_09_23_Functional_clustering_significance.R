source("R/Analysis_data.R")
source("R/Common variables.R")
library("functClust")
library(tidyverse)
library(ggsignif)

PROPERTIES <- c("productivity","DeltaY","Selection","Complementarity")
options(verbose = TRUE) # So that I see the advance of the computations

nbcl = 2
site <- "GrandeDixence"
property <- "DeltaY"
# res <- fclust_read(filename = paste0("figures/functClust/functclust/15_species/",site,"_",property))
# 
# fclust_plot(res, main = paste(site,property),
#             opt.tree = list("prd"),
#             nbcl=nbcl,
#             opt.motif = list("obs",  "aov", pvalue = 0.05))

res <- fclust_read(filename = paste0("figures/functClust/functclust/12_to_18_species/",site,"_",property))

fclust_plot(res, main = paste(site,property),
            opt.tree = list("prd"),
            nbcl=nbcl,
            opt.motif = list("obs",  "aov", pvalue = 0.05))

#____________________________________________
# Read the files printed by fclust_write
# and write plots of observed performance

site <- "Bern"
property <- "Complementarity"
nbCl <- 2 # number of clusters of species

nb_sp <- 15

folder <- if (nb_sp == "all"){
  "all_species_numbers"
} else{
  paste0(nb_sp,"_species")
}

for (site in SITE){
  for (property in PROPERTIES){
    inputs <- read.csv( paste0("figures/functClust/functClust/",folder,"/",site,"_",property,".inputs.csv") ) %>% 
      rename(Ass = names.fobs)
    matrices <- read.csv( paste0("figures/functClust/functClust/",folder,"/",site,"_",property,".matrices.csv") ) 
    trees <-read.csv( paste0("figures/functClust/functClust/",folder,"/",site,"_",property,".trees.csv") ) 
    
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
      gather(key = "Species", value = "Cluster") %>% 
      filter(grepl(".II",Species)) %>% 
      filter(!grepl("R2.II",Species)) %>%  # remove R2
      mutate(Species = str_replace_all(Species, ".II", "") ) %>%
      group_by(Cluster) %>% 
      mutate(nb_sp_in_cluster = n()) %>% 
      ungroup() %>% 
      mutate(more_sp = max(nb_sp_in_cluster)) %>% 
      mutate(more_sp = if_else (more_sp == nb_sp_in_cluster, T, F)) %>% 
      mutate(Cluster = if_else(more_sp,"b", "a"))
    
    
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
        pull(Cluster) %>% 
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
    
    #species in cluster a
    sp_in_a <- affectSp %>% 
      filter(Cluster=="a") %>% 
      pull(Species)
    
    sp_in_b <- affectSp %>% 
      filter(Cluster=="b") %>% 
      pull(Species)
    
    obsNb <- fobs_per_motif %>% 
      group_by(Motif) %>% 
      count()
    
    plot_obs_perf <- ggplot(fobs_per_motif,aes(x=Motif,y=fobs))+
      geom_boxplot() +
      labs(title = paste(site,property),
           subtitle = paste("Species in cluster a =",
                            gsub(", ",", ",toString(sp_in_a)), "- out of",length(sp_in_a) + length(sp_in_b),"species" ) ,
           caption= paste("n(", obsNb[1,1],") =",obsNb[1,2],"of",sum(pull(obsNb,n)), "assemblages"  ) ) +
      geom_signif(
        comparisons = list(c("ab", "b")),
        map_signif_level = TRUE
      )
    
    ggsave(plot = plot_obs_perf, filename = paste0("figures/functClust/",nb_sp,"_species/",site,"_",property ,".png") )
    
    
    
  }
}
