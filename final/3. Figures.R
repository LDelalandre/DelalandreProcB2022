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
      scale_x_continuous(breaks = 2*c(1:15)) +
      ggsave(paste0("figures/",measure,"_",site,".png"))
  }
}

sit <- "Bever"
orde <- "decreasing"
data <- read.table("data/processed/specific_biom_prod_complete.txt",header=T)
sub <- subset(data,site==sit&order==orde)
ggplot(sub, aes(x=simul, y=mixture.t.ha., fill = species)) +
  geom_area()

sub$species <- factor(sub$species,levels=)
ggplot(sub, aes(simul, mixture.t.ha., group = species)) +
  geom_area(aes(fill = species),position = "stack") #+
  geom_text(aes(label = species), position = position_stack(vjust = 0.5))
  
ggplot(sub, aes(simul, mixture.t.ha.)) +
  geom_area(aes(fill = species))
