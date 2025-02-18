# Remove abundant species first ####
PLOT <- list()
i <- 0
for (site in ord_plots){ # ATTENTION RE FAIRE TOURNER COTTBUS
  i <- i+1
  PROD <- NULL
  order <- "abundance"
  for(number in c(1:30)){
    prod <- try(read.table(paste0("data/raw/output-cmd2_",site,"_",order,".txt/forceps.",site,".site_",number,"_productivityScene.txt")),silent=T)
    if (class(prod) != "try-error"){# sometimes, the files are empty, and it returns an error message
      colnames(prod)<-colnames_prod
      years_to_keep <- max(prod$date) - c(900,800,700,600,500,400,300,200,100,0)
      prod_to_keep <- subset(prod,date %in% years_to_keep)
      prod_to_keep$totProdBiomass_t_ha <- prod_to_keep$adultProdBiomass_t_ha + prod_to_keep$saplingBiomass_t_ha
      
      # specific productivities
      productivities <- aggregate(prod_to_keep$totProdBiomass_t_ha, list(prod_to_keep$speciesShortName), mean)
      colnames(productivities) <- c("species","mixture_t_ha")
      productivities$site <- site
      productivities$order <- order
      productivities$simul <- number
      productivities$mixture_relative <- productivities$mixture_t_ha/sum(productivities$mixture_t_ha)
      PROD <- rbind(PROD,productivities)
    }
    
  }
  PROD_sp <- PROD
  
  
  
  PROD_tot <- aggregate(PROD_sp$mixture_t_ha,list(site = PROD_sp$site,order = PROD_sp$order,simul = PROD_sp$simul),FUN=sum)# sum of the productivities of the species
  
  prod_per_order <- spread(PROD_tot,order,x) # data frame with each order in column
  
  measure <- "productivity_tot"
  result <- read.table(paste0("data/processed/",measure,"_",site,"_with interval_median.txt"),header=T)
  abundance <- prod_per_order$abundance
  result <- cbind(result,abundance)
  
  plot <- ggplot(result,aes(x=simul-1,y=abundance,color="green")) +
    scale_color_manual(values=c("green4")) +
    labs(x="Number of species removed",y=measure) +
    geom_line()+
    geom_ribbon(aes(ymin=int_min, ymax=int_max),fill="grey60", alpha=0.5,colour="black") +
    theme(legend.position = "bottom") +
    ggtitle(site) +
    theme(plot.title = element_text(size=18)) +
    scale_x_continuous(breaks = 5*c(1:6)) +
    theme(legend.title = element_blank())+
    theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none" ) + # virer tous les titres
    theme(axis.text=element_text(size=15))
  
  PLOT[[i]] <- plot
}
plot_ab <- plot_grid(PLOT[[1]],PLOT[[2]],PLOT[[3]],PLOT[[4]],PLOT[[5]],PLOT[[6]],PLOT[[7]],PLOT[[8]],PLOT[[9]],PLOT[[10]],PLOT[[11]], ncol = 4, nrow = 3)
save_plot("paper2/Abundant species first.png",plot_ab,base_height = 10)