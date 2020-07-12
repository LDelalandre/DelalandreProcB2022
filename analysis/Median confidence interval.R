
pbinom(9, size=30, prob=.5) # On cherche k (ici k=9) tel que la probabilité d'avoir moins de k succès soit de 0.025
qbinom(.025, size=30, prob=.5) # en fait il faut chercher dans ce sens. On cherche le 2.5ème quantile d'une binomiale avec 30 tirages et une proba de succès de 0.5

for (sit in SITE){
  table <- read.table(paste0("data/processed/productivity_tot_",sit,"_with interval.txt"),header=T)
  table2 <- select(table,starts_with("random")) 
  
  int_min <- c()
  int_max <- c()
  median <- c()
  for(i in c(1:30)){
    int_min <- c(int_min, table2[i,order(table2[i,])[10]  ] ) # avec k = 10
    int_max <- c(int_max, table2[i, order(table2[i,])[21] ] ) # avec 21 = 30-10+1 = N-k+1
    median <- c(median,median(as.numeric(table2[i,]),na.rm=T) )
  }
  table$int_min <- int_min
  table$int_max <- int_max
  table$mean <- median
  
  write.table(table,paste0("data/processed/productivity_tot_",sit,"_with interval_median.txt"))
  
}


# plot ####
measure <- "productivity_tot"
result <- read.table(paste0("data/processed/",measure,"_",sit,"_with interval_median.txt"),header=T)
plot <-
  ggplot(result,aes(x=simul-1,y=decreasing,color="Removing distinct species first")) +
  labs(x="Number of species removed",y=measure) +
  # geom_line(size=1)+
  geom_ribbon(aes(ymin=int_min, ymax=int_max),fill="grey60", alpha=0.5,colour="black") +
  geom_line(aes(x=simul-1,y=increasing, color="#8766D"),size=2) +
  geom_line(aes(x=simul-1,y=decreasing, color="#00BFC4"),size=2) +
  theme(legend.position = "bottom") +
  ggtitle(sit) +
  scale_x_continuous(breaks = 2*c(1:15)) +
  theme(legend.title = element_blank())+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none" ) # virer tous les titres
# ggsave(paste0("figures/",measure,"_",site,".png"))
plot
