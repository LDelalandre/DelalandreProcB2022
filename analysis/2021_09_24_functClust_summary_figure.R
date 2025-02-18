source("R/Analysis_data.R")
source("R/Common variables.R")
library("functClust")
library(tidyverse)
library(ggsignif)

PROPERTIES <- c("productivity","DeltaY","Selection","Complementarity")
site <- "Bern"
property <- "Complementarity"
nbCl <- 2 # number of clusters of species
siteNb <- 0 # to save pots in the good order
nb_sp <- 30 # number of species 

for (site in SITE){
  siteNb <- siteNb + 1
  PLOT <- list()
  i <- 1
  for (property in PROPERTIES){

    
    folder <- if (nb_sp == "all"){
      "all_species_numbers"
    } else{
      paste0(nb_sp,"_species")
    }
    
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
    distinct_sp <- read.table("data/raw/distinctiveness of the species.txt",header=T) %>% 
      arrange(desc(Di)) %>% 
      slice(1:10) %>% 
      pull(SName)
    
    common_sp <- read.table("data/raw/distinctiveness of the species.txt",header=T) %>% 
      arrange(desc(Di)) %>% 
      slice(11:30) %>% 
      pull(SName)
    
    
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
      mutate(Cluster = if_else(more_sp,"b", "a")) %>% 
      mutate(Dist = if_else(Species %in% distinct_sp,"Di","Com") )
    
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
    fobs_per_assemblage <- merge(inputs,motif_of_assemblages,by = "Ass")
    
    #species in cluster a
    sp_in_a <- affectSp %>% 
      filter(Cluster=="a") %>% 
      pull(Species)
    
    sp_in_b <- affectSp %>% 
      filter(Cluster=="b") %>% 
      pull(Species)
    
    nb_dist_per_cluster <- affectSp %>% 
      group_by(Cluster) %>% 
      mutate(n=n(),
             dist = sum(Dist=="Di"),
             common = sum(Dist=="Com")) %>% 
      mutate(prop_di = dist/n) %>% 
      mutate(Dist = factor(Dist,levels = c("Com","Di")))
    
    
    
    # Summarize the info ####
    # Nb of observation per motif
    obsNb <- fobs_per_assemblage %>% 
      group_by(Motif) %>% 
      count()
    
    # Nb of species per cluster
    spNb <- affectSp %>% 
      group_by(Cluster) %>% 
      summarize(Cluster=Cluster,nb_sp = nb_sp_in_cluster) %>% 
      distinct()
    
    # Function observed per assemblage
    fobs_per_motif <- fobs_per_assemblage %>% 
      select(Ass,fobs,names.xpr,Motif) %>% 
      group_by(Motif) %>% 
      select(-Ass) %>% 
      summarize(meanfobs = mean(fobs),property=names.xpr) %>% 
      summarize(n = n()) # number of observations per motif
    
    plot_sp_per_cluster <- ggplot(nb_dist_per_cluster,aes(x=Cluster,fill=Dist))+
      geom_bar() +
      scale_fill_manual(name = "",
                        values = c("#00BFC4","#F8766D"),
                        labels = c("Common species", "Distinct species")) +
      labs(title = paste(site,property)) 
    
    # legend
    leg <- ggpubr::get_legend(plot_sp_per_cluster)
    ggpubr::as_ggplot(leg)
    
    plot_sp_per_cluster_no_legend <- plot_sp_per_cluster +
      theme(legend.position = "none" )
    
    
    # ggsave(plot= plot_sp_per_cluster,
    #        filename = paste0("figures/functClust/",nb_sp,"_species/nb Di per cluster/nb_di_",site,"_",property ,".png") )
    # 
    
    #____________________________________________________________________
    # plot
    
    fobs_per_motif2 <- merge(inputs,motif_of_assemblages,by = "Ass")
    
    plot_obs_perf <- ggplot(fobs_per_motif2,aes(x=Motif,y=fobs))+
      geom_boxplot() +
      labs(title = paste(site,property),
           subtitle = paste("Species in cluster a =",
                            gsub(", ",", ",toString(sp_in_a)), "- out of",length(sp_in_a) + length(sp_in_b),"species" ) ,
           caption= paste("n(", obsNb[1,1],") =",obsNb[1,2],"of",sum(pull(obsNb,n)), "assemblages"  ) ) +
      geom_signif(
        comparisons = list(c("ab", "b")),
        map_signif_level = TRUE
      )
    
    # ggsave(plot = plot_obs_perf, filename = paste0("figures/functClust/",nb_sp,"_species/",site,"_",property ,".png") )
    # 
    
    
    PLOT[[i]] <- plot_obs_perf
    i <- i+1
    PLOT[[i]] <- plot_sp_per_cluster_no_legend
    i <- i+1
  }
  
  
  
  combined_plots <- grid.arrange( PLOT[[1]],PLOT[[2]], PLOT[[3]],
                                  PLOT[[4]], PLOT[[5]],PLOT[[6]],
                                  PLOT[[7]],PLOT[[8]],
                                  ncol = 4, nrow = 2, 
                                  layout_matrix = rbind( c(1,3,5,7),c(2,4,6,8)))
  
  ggsave(paste0("figures/functClust/",siteNb,"combined_plot_",nb_sp,"_sp",site,".png")
         ,combined_plots,width = 15, height = 10)
}


