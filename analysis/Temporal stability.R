
measure <- "TS_productivity_tot"
for (site in SITE){
  sd_prod_tot <- sd_productivity_tot(site)
  TSsite <- subset(sd_prod_tot,site==site)
  TSsite <- select(TSsite,TS,site,order,simul)
  write.table(spread(TSsite,order,TS),paste0("data/processed/",measure,"_",site,".txt"))
  
  conf <- confidence_interval(site,measure)
  write.table(conf,paste0("data/processed/",measure,"_",site,"_with interval.txt"))
  
  result <- read.table(paste0("data/processed/",measure,"_",site,"_with interval.txt"),header=T)
  ggplot(result,aes(x=simul-1,y=decreasing,color="Removing distinct species first")) +
    labs(x="Number of species removed",y=measure) +
    geom_line()+
    geom_ribbon(aes(ymin=int_min, ymax=int_max),fill="grey60", alpha=0.5,colour="black") +
    geom_line(aes(x=simul-1,y=increasing, color="Removing distinct species last")) +
    theme(legend.position = "bottom") +
    ggtitle(site) +
    scale_x_continuous(breaks = 2*c(1:15)) +
    ggsave(paste0("figures/",measure,"_",site,".png"))+
    theme(legend.title = element_blank())
}

