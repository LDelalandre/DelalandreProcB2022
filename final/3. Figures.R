source("R/Common variables.R")

# plot the final ecosystem property in each condition
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
