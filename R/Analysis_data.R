library("dplyr")

# Read and process raw data ####
temporal_plot<-function(res){ 
  # Returns the biomass and the abundance of each species, every year
  # Compute species biomass (proxy for abundance) for the different years exported from the simulations
  # res is the data frame of the results complete (e.g. "forceps.Bern.site_1_complete.txt")
  res$species_date<-paste(res$speciesShortName,res$date)
  temp_plot<-data.frame(row.names=rownames(tapply(res$biomass.kg.,res$species_date,sum)))
  temp_plot$biomass <- tapply(res$biomass.kg.,res$species_date,sum)
  temp_plot$species <- rownames(temp_plot)
  library(tidyr)
  temp_plot<-separate(temp_plot, col = species , sep = " ", into = c("species", "date"),remove = F)
  
  # compute species abundance (number of individuals) for the different years
  res$nb<-matrix(data=1,nrow = dim(res)[1],ncol=1) # each indiv of each sp in each year... is one individual.
  temp_plot$abundance<-tapply(res$nb,res$species_date,sum)
  
  temp_plot
}


temporal_plot_threshold<-function(temp_plot){ 
  # Biomass and abundance of each species above an abundance threshold, every year
  # Threshold: every year, we keep the species whose biomass is above 0.001 times the total biomass of that year.
  # temp_plot is the output of function temporal_plot
  biomass_year <- c()
  temp_plot_threshold<-temp_plot
  c=1
  for (i in unique(temp_plot$date)){
    biomass_year<-c(biomass_year, sum( filter(temp_plot,date==i)$biomass ) )
    temp_plot_threshold<-filter(temp_plot_threshold, (date !=i) | (date==i & biomass>=0.001*biomass_year[c]) )
    c=c+1
  }
  temp_plot_threshold
}



richness<-function(temp_plot){
  # Richness and its changes in time
  # Returns a list of the richness value of every year
  rich<-c()
  dates<-as.character(sort(unique(as.numeric(temp_plot$date)))) # So that the years are in growing order
  for (i in dates ){
    rich<-c(rich, dim(subset(temp_plot,date==i))[1] )
  }
  rich
}


# Biomass and sd(biomass) - species removal #####
removal_exp_final_biomasses<- function(site) {
  # Takes the output of the simulations (complete files). 
  # Extracts final total biomass of the community for each number of species and each condition (decreasing or increasing order of distinctiveness)
  # Writes a data frame with these informations in a .txt file
  
  # Data for having the names of the columns
  colnames_mean<-colnames(read.table("data/colnames_mean.txt",header=T)) # idem
  colnames_res<-colnames(read.table("data/colnames_res.txt", header=T))
  
  richness<-c(30:1)
  
  # Species removed with increasing functional distinctiveness
  biomass_inc<-c()
  sd_biom_inc<-c()
  for(i in c(1:30)){
    mean<-read.table(paste0("data/raw/output-cmd2_",site,"_increasing.txt/forceps.",site,".site_",i,"_mean.txt"))
    colnames(mean)<-colnames_mean
    # biomass averaged
    years_to_keep <- max(mean$date) - c(900,800,700,600,500,400,300,200,100,0)
    meanbiom <- mean(subset(mean,date %in% years_to_keep)$totalBiomass.t.ha.)
    biomass_inc<-c(biomass_inc,meanbiom)
    sd_biom_inc<-c(sd_biom_inc,sd(mean$totalBiomass.t.ha.))
  }
  
  
  
  # Species removed with decreasing functional distinctiveness
  biomass_dec<-c()
  sd_biom_dec<-c()
  for(i in c(1:30)){
    mean<-read.table(paste0("data/raw/output-cmd2_",site,"_decreasing.txt/forceps.",site,".site_",i,"_mean.txt"))
    colnames(mean)<-colnames_mean
    # biomass averaged
    years_to_keep <- max(mean$date) - c(900,800,700,600,500,400,300,200,100,0)
    meanbiom <- mean(subset(mean,date %in% years_to_keep)$totalBiomass.t.ha.)
    biomass_dec<-c(biomass_dec,meanbiom)
    sd_biom_dec<-c(sd_biom_dec,sd(mean$totalBiomass.t.ha.))
  }
  
  # Species removed in random order
  nb_rand <- 10 # number of random simul
  RAND <- data.frame(matrix(NA,nrow=30,ncol=nb_rand))
  sd_RAND <- RAND
  for (j in c(1:nb_rand)){
    biomass_ran<-c()
    sd_biomass_ran<-c()
    for(i in c(1:30)){
      mean<-read.table(paste0("data/raw/output-cmd2_",site,"_random_",j,".txt/forceps.",site,".site_",i,"_mean.txt"))
      colnames(mean)<-colnames_mean
      # biomass averaged
      years_to_keep <- max(mean$date) - c(900,800,700,600,500,400,300,200,100,0)
      meanbiom <- mean(subset(mean,date %in% years_to_keep)$totalBiomass.t.ha.)
      biomass_ran <- c(biomass_ran, meanbiom)
      sd_biomass_ran <- c(sd_biomass_ran, sd(mean$totalBiomass.t.ha.))
    }
    RAND[,j]<-biomass_ran
    sd_RAND[,j]<- sd_biomass_ran
  }
  
  
  # Je construis un intervalle de confiance comme si les observations suivaient une loi normale... A am?liorer.
  int_min <- c()
  int_max <- c()
  mean <- c()
  for (j in c(1:30)){
    int_min <- c(int_min, mean(as.numeric(RAND[j,])) + 1.96 * sd(as.numeric(RAND[j,]))/sqrt(10) )
    int_max <- c(int_max, mean(as.numeric(RAND[j,])) - 1.96 * sd(as.numeric(RAND[j,]))/sqrt(10) )
    mean <- c(mean, mean(as.numeric(RAND[j,])) )
  }
  RAND$int_min <- int_min
  RAND$int_max <- int_max
  RAND$mean <- mean
  RAND$biomass_dec <- biomass_dec
  RAND$biomass_inc <- biomass_inc
  RAND$nb_removed <- c(0:29)
  
  # Same thing for sd(biomass). NB: it is computed on the biomass along ALL the simulation, includind the succession phase.
  sd_int_min <- c()
  sd_int_max <- c()
  sd_mean <- c()
  for (j in c(1:30)){
    sd_int_min <- c(sd_int_min, mean(as.numeric(sd_RAND[j,])) + 1.96 * sd(as.numeric(sd_RAND[j,]))/sqrt(10) )
    sd_int_max <- c(sd_int_max, mean(as.numeric(sd_RAND[j,])) - 1.96 * sd(as.numeric(sd_RAND[j,]))/sqrt(10) )
    sd_mean <- c(sd_mean, mean(as.numeric(sd_RAND[j,])) )
  }
  sd_RAND$int_min <- sd_int_min
  sd_RAND$int_max <- sd_int_max
  sd_RAND$mean <- sd_mean
  sd_RAND$sd_biomass_dec <- sd_biom_dec
  sd_RAND$sd_biomass_inc <- sd_biom_inc
  sd_RAND$nb_removed <- c(0:29)
  
  write.table(RAND,paste0("data/processed/Biomass_species removal experiments_",site,".txt"))
  write.table(sd_RAND,paste0("data/processed/sd_biomass_species removal experiments_",site,".txt"))
}


# Productivity - removal experiments ####
removal_exp_productivity <- function(site){
  colnames_prod <- read.table("data/colnames_productivityScene.txt",header=T)
  
  richness<-c(30:1)
  
  # Species removed with increasing functional distinctiveness
  prod_inc<-c()
  for(i in c(1:30)){
    prod <- read.table(paste0("data/raw/output-cmd2_",site,"_increasing.txt/forceps.",site,".site_",i,"_productivityScene.txt"))
    colnames(prod)<-colnames(colnames_prod)
    # productivity averaged
    years_to_keep <- max(prod$date) - c(900,800,700,600,500,400,300,200,100,0)
    prod_to_keep <- subset(prod,date %in% years_to_keep)
    # Annual productivity = sum of species specific annual productivities:
    annual_prod <- c()
    for(k in years_to_keep){
      annual_prod <- c(annual_prod, sum(prod_to_keep[which(prod_to_keep$date==k),]$adultProdBiomass) + sum(prod_to_keep[which(prod_to_keep$date==k),]$saplingBiomass))
    }
    prod_inc<-c(prod_inc,mean(annual_prod))
  }
  
  # Species removed with decreasing functional distinctiveness
  prod_dec<-c()
  for(i in c(1:30)){
    prod <- read.table(paste0("data/raw/output-cmd2_",site,"_decreasing.txt/forceps.",site,".site_",i,"_productivityScene.txt"))
    colnames(prod)<-colnames(colnames_prod)
    # productivity averaged
    years_to_keep <- max(prod$date) - c(900,800,700,600,500,400,300,200,100,0)
    prod_to_keep <- subset(prod,date %in% years_to_keep)
    # Annual productivity = sum of species specific annual productivities:
    annual_prod <- c()
    for(k in years_to_keep){
      annual_prod <- c(annual_prod, sum(prod_to_keep[which(prod_to_keep$date==k),]$adultProdBiomass) + sum(prod_to_keep[which(prod_to_keep$date==k),]$saplingBiomass))
    }
    prod_dec<-c(prod_dec,mean(annual_prod))
  }

  
  # Species removed in random order
  nb_rand <- 10 # number of random simul
  RAND <- data.frame(matrix(NA,nrow=30,ncol=nb_rand))
  
  for (j in c(1:nb_rand)){
    prod_ran<-c()
    for(i in c(1:30)){
      prod <- read.table(paste0("data/raw/output-cmd2_",site,"_random_",j,".txt/forceps.",site,".site_",i,"_productivityScene.txt"))
      colnames(prod)<-colnames(colnames_prod)
      # productivity averaged
      years_to_keep <- max(prod$date) - c(900,800,700,600,500,400,300,200,100,0)
      prod_to_keep <- subset(prod,date %in% years_to_keep)
      # Annual productivity = sum of species specific annual productivities:
      annual_prod <- c()
      for(k in years_to_keep){
        annual_prod <- c(annual_prod, sum(prod_to_keep[which(prod_to_keep$date==k),]$adultProdBiomass)+ sum(prod_to_keep[which(prod_to_keep$date==k),]$saplingBiomass) )
      }
      prod_ran<-c(prod_ran,mean(annual_prod))
    }
    RAND[,j]<-prod_ran
  }
  
  # Confidence interval (NB: to improve!) and .txt export
  int_min <- c()
  int_max <- c()
  mean <- c()
  for (j in c(1:30)){
    int_min <- c(int_min, mean(as.numeric(RAND[j,])) + 1.96 * sd(as.numeric(RAND[j,]))/sqrt(10) )
    int_max <- c(int_max, mean(as.numeric(RAND[j,])) - 1.96 * sd(as.numeric(RAND[j,]))/sqrt(10) )
    mean <- c(mean, mean(as.numeric(RAND[j,])) )
  }
  RAND$int_min <- int_min
  RAND$int_max <- int_max
  RAND$mean <- mean
  RAND$prod_dec <- prod_dec
  RAND$prod_inc <- prod_inc
  RAND$nb_removed <- c(0:29)

  write.table(RAND,paste0("data/processed/Productivity_species removal experiments_",site,".txt"))
}
  
# Price format ####

table_threshold_price <- function(site,order){
  # Writes a table with the final biomass and abundance of each species, the number of the simulation, and the order of removal
  # Simulation 1: all the species are present in the regional pool.
  # simul = 30: juste one species remains in the regional pool.
  # Order is an element of c("increasing", "decreasing", "random_1" ,  "random_2")
  res<-read.table(paste0("data/raw/output-cmd2_",site,"_",order,".txt/forceps.Bern.site_1_complete.txt"),header=F) 
  colnames(res) <- colnames(read.table(here::here("data","colnames_res.txt"),header=T))
  colnames(res)<-colnames_res
  temp_plot<-filter(temporal_plot(res),date==max(unique(res$date)))
  temp_plot$simul <- rep(1,dim(temp_plot)[1])
  data<-temporal_plot_threshold(temp_plot)
  for(i in c(2:30)){
    res<-read.table(paste0("data/raw/output-cmd2_",site,"_",order,".txt/forceps.Bern.site_",i,"_complete.txt"),header=F) 
    colnames(res)<-colnames_res
    temp_plot<-filter(temporal_plot(res),date==max(unique(res$date)))
    temp_plot$simul <- rep(i,dim(temp_plot)[1])
    data<-rbind(data,temporal_plot_threshold(temp_plot))
  }
  data$order <- rep(order,dim(data)[1])
  write.table(data,paste0("data/processed/table_price_threshold_",order,".txt"))
}
