source("R/Common variables.R")

measure <- MEASURE[2]#5
# plot the final biomasses in each condition
# more precisely, it is the averaged biomass on 10 evenly-spaced points in time during the last 1000 years.
for (site in SITE){
  for (measure in MEASURE){
    result <- read.table(paste0("data/processed/",measure,"_",site,"_with interval.txt"),header=T)
    ggplot(result,aes(x=simul-1,y=decreasing,color="Removing distinct species first")) +
      labs(x="Number of species removed",y=measure) +
      geom_line()+
      geom_ribbon(aes(ymin=int_min, ymax=int_max),fill="grey60", alpha=0.5,colour="black") +
      geom_line(aes(x=simul-1,y=increasing, color="Removing distinct species last")) +
      theme(legend.position = "bottom") +
      ggtitle(site) +
      ggsave(paste0("figures/",measure,"_",site,".png"))
  }
}


# plot sd(biomass)
for (site in SITE){
  result <- read.table(paste0("data/processed/sd_biomass_species removal experiments_",site,".txt"))
  ggplot(result,aes(x=nb_removed,y=sd_biomass_dec,color="Removing distinct species first")) +
    labs(x="Number of species removed",y="sd(Total biomass)") +
    geom_line()+
    geom_ribbon(aes(ymin=int_min, ymax=int_max),fill="grey60", alpha=0.5,colour="black") +
    geom_line(aes(x=nb_removed,y=sd_biomass_inc, color="Removing distinct species last")) +
    theme(legend.position = "bottom") +
    ggtitle(site) +
    ggsave(paste0("figures/sd_biomass_",site,".png"))
}



# plot productivity
for (site in SITE2){
  result <- read.table(paste0("data/processed/Productivity_species removal experiments_",site,".txt"))
  ggplot(result,aes(x=nb_removed,y=prod_dec,color="Removing distinct species first")) +
    labs(x="Number of species removed",y="Productivity") +
    geom_line()+
    geom_ribbon(aes(ymin=int_min, ymax=int_max),fill="grey60", alpha=0.5,colour="black") +
    geom_line(aes(x=nb_removed,y=prod_inc, color="Removing distinct species last")) +
    theme(legend.position = "bottom") +
    ggtitle(site) +
    ggsave(paste0("figures/Productivity_",site,".png"))
}