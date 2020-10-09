source("R/Common variables.R")


# write.table(data,paste0("data/processed/specific_",measure,"_with_zeros.txt"),row.names=F)
# now we have all the data, with zero 
# évidemment, c'est méééééga long, vu que c'est codé avec les pieds...

####
# NB: ajouter une option: relatif/absolu
TOTAL <- read.table("data/processed/specific_biom_prod_complete.txt",header=T)
measure <- "biomass_tot"
relatif=F

LH <- function(measure,relatif){
  # relatif=T or F
  if (measure=="biomass_tot"){
    if (relatif == F){
      data <- select(TOTAL,species, mixture.t.ha.,monoculture.t.ha.,site, order, simul)
      names(data)[names(data) == "mixture.t.ha."] <- "YOi"
      names(data)[names(data) == "monoculture.t.ha."] <- "Mi"
    } else{
      data <- select(TOTAL,species, mixture_relative,monoculture_relative,site, order, simul)
      names(data)[names(data) == "mixture_relative"] <- "YOi"
      names(data)[names(data) == "monoculture_relative"] <- "Mi"
    }

  }else{
    if (relatif == F){
      data <- select(TOTAL,species, prod_mixture,prod_monoculture,site, order, simul)
      names(data)[names(data) == "prod_mixture"] <- "YOi"
      names(data)[names(data) == "prod_monoculture"] <- "Mi"
    }
    data <- select(TOTAL,species, prod_mixture_relative,prod_monoculture_relative,site, order, simul)
    names(data)[names(data) == "prod_mixture_relative"] <- "YOi"
    names(data)[names(data) == "prod_monoculture_relative"] <- "Mi"
  }
  
  # dat <- read.table("data/processed/specific_biom_prod_with_zeros.txt",header=T)
  data$un <- rep(1,dim(data)[1])
  data <- data %>% group_by(site,order,simul)
  
  
  # data2 <- data %>%
  #   mutate(YO=sum(YOi)) %>%
  #   mutate(RYEi=1/(31-simul)) %>%# 1/N, N being the initial nb of species. For simul 30, N=31-30=1 sp. For simul 1, N=31-1=30 sp, etc.
  #   mutate(RYOi=YOi/Mi) %>%
  #   mutate(YEi=RYEi*Mi) %>%
  #   mutate(YE=sum(YEi)) %>%
  #   mutate(DeltaY=YO-YE) %>%
  #   mutate(DeltaRYi=RYOi-RYEi) %>%
  #   mutate(Mavg = mean(Mi)) %>%
  #   mutate(DeltaRYavg = mean(DeltaRYi)) %>%
  #   
  #   mutate(Cpltarity = (31-simul)*DeltaRYavg*Mavg) %>%
  #   mutate(Selection = DeltaY - Cpltarity)%>%
  #   
  #   mutate(Selection2=(31-simul)*(DeltaRYavg-DeltaRYi)*(Mavg-Mi))
  
  data2 <- data%>%
    subset(Mi!=0) %>%
    mutate(nb_sp=sum(un)) %>%
    mutate(YO=sum(YOi)) %>%
    mutate(RYEi=1/nb_sp) %>%# 1/N, N being the initial nb of species. For simul 30, N=31-30=1 sp. For simul 1, N=31-1=30 sp, etc.
    mutate(RYOi=YOi/Mi) %>%
    mutate(YEi=RYEi*Mi) %>%
    mutate(YE=sum(YEi)) %>%
    mutate(DeltaY=YO-YE) %>%
    mutate(DeltaRYi=RYOi-RYEi) %>%
    mutate(Mavg = mean(Mi)) %>%
    mutate(DeltaRYavg = mean(DeltaRYi)) %>%
    
    mutate(Cpltarity = nb_sp*DeltaRYavg*Mavg) %>%
    mutate(Selection = DeltaY - Cpltarity)%>%
    
    mutate(Selection2=nb_sp*(DeltaRYavg-DeltaRYi)*(Mavg-Mi))
  data2
}

# Compute loreau-hector values for relative and absolute biomass and productivity
for (relatif in c(T,F)){
  for (measure in MEASURE){
    data2 <- LH(measure,relatif)
    data3 <- summarise(data2,DeltaY=mean(DeltaY),Cpltarity=mean(Cpltarity),Selection=mean(Selection),Selection2=mean(Selection2))
    write.table(data3,paste0("data/processed/Loreau-Hector_",measure,"_relatif=",relatif,".txt"),row.names = F,sep="\t")
    # I use mean, because all the values are already the same for one group for DeltaRY etc.
  }
}


# sit <- "Bever"
# orde <- ORDER[1]
# simu <- 1
# LH evolution with species removal within site
for (relatif in c(T,F)){
  for (measure in MEASURE){
    data3 <- read.table(paste0("data/processed/Loreau-Hector_",measure,"_relatif=",relatif,".txt"),header=T)
    for (sit in SITE){
      for (orde in ORDER[1:2]){
        within_site <- subset(data3,site==sit & order == orde)
        # plot(toplot$simul,toplot$Selection,type="l")
        # lines(toplot$simul,toplot$Cpltarity,type="l",col="2")#,ylim=c(-10000,300))
        
        p <- ggplot(within_site,aes(x=simul-1,y=Selection,color="Selection")) +
          labs(x="Number of species removed",y=paste("LH coefficients",measure,"relatif=",relatif)) +
          geom_line()+
          # geom_ribbon(aes(ymin=int_min, ymax=int_max),fill="grey60", alpha=0.5,colour="black") +
          geom_line(aes(x=simul-1,y=Cpltarity, color="Complementarity")) +
          theme(legend.position = "bottom")
        # scale_x_continuous(breaks = 2*c(1:15)) +
        
        if (orde=="decreasing"){
          p=p+ggtitle(paste(sit,"Removing distinct species first"))
        } else if (orde=="increasing") {
          p=p+ggtitle(paste(sit,"Removing distinct species last"))
        }
        p + ggsave(paste0("figures/Loreau-Hector/within_site_",measure,"_relatif=",relatif,"_",sit,"_",orde,".png"))
      }
    }
  }
}

# Sort sites according to their hydric stress index ####
drought_index <- c()
for (site in SITE){
  mean <- read.table(paste0("data/raw/output-cmd2_",site,"_decreasing.txt/forceps.",site,".site_1_mean.txt"))
  colnames(mean) <- colnames_mean
  to_keep <- subset(mean,date %in% (as.numeric(max(mean$date))- c(900,800,700,600,500,400,300,200,100,0))) 
  drought_index <- c(drought_index,mean(to_keep$droughtIndexAnnual))
}
DROUGHT <- data.frame(SITE,drought_index)
write.table(DROUGHT,"data/drought index.txt",row.names=F)

# LH values for the first simulation of each site ####
for (relatif in c(T,F)){
  for (measure in MEASURE){
    for (Sorted_by in c("Temperature","Precipitation","Drought_index")){
      data3 <- read.table(paste0("data/processed/Loreau-Hector_",measure,"_relatif=",relatif,".txt"),header=T)
      across_sites <- subset(data3,order=="increasing" & simul==1) # the first simul is the same whether the order is increasing or decreasing
      if (Sorted_by == "Temperature"){
        ord <- Site_descr[order(Site_descr$Temp_moy),]$Site # Sites ordered with increasing mean temperature
      } else if (Sorted_by=="Precipitation") {
        ord <- Site_descr[order(Site_descr$Annual_ppt),]$Site # Sites ordered with increasing mean temperature
      } else {
        ord <- DROUGHT[order(DROUGHT$drought_index),]$SITE
      }
      across_sites$site <- factor(across_sites$site,levels = ord )
      
      p2 <- ggplot(across_sites,aes(x=site,y=Selection,color="Selection")) +
        labs(x=paste("Site sorted by",Sorted_by),y=paste("LH coefficients",measure,"relatif=",relatif)) +
        geom_boxplot()+
        # geom_ribbon(aes(ymin=int_min, ymax=int_max),fill="grey60", alpha=0.5,colour="black") +
        geom_boxplot(aes(x=site,y=Cpltarity, color="Complementarity")) +
        theme(legend.position = "bottom") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      p2 + ggsave(paste0("figures/Loreau-Hector/across_sites_",measure,"_relatif=",relatif,"_sorted by_",Sorted_by,".png"))
    }
  }
}



plot(within_site$DeltaY~within_site$simul,type="l")

data3 <- summarise(data2,DeltaY=mean(DeltaY),YO=mean(YO),YE=mean(YE),Cpltarity=mean(Cpltarity),Selection=mean(Selection),Selection2=mean(Selection2))
plot(data3$YE~data3$YO)
abline(0,1)
