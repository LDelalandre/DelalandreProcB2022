source("R/Common variables.R")

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


# Biomass - species removal #####
specific_biomass_final <- function(site) {
  # 1. Takes the output of the simulations (complete files).
  # 2. Selects biomass above threshold
  # 3. Selects a subset of ten years and averages the biomass of each species on these years
  # 4. Returns data frame (e.g. specific_biomass_final_Bern.txt) whose columns are:
  # "species"	"mixture(t/ha)"	"mixture_relative"	"site"	"order"	"simul"
  BIOMASSES <- as.data.frame(matrix(nrow=0,ncol=4,dimnames = list(NULL,c("species", "mixture(t/ha)", "mixture_relative", "simul"))))
  for (order in ORDER){
    biomass_inc<-c()
    sd_biom_inc<-c()
    for(number in c(1:30)){
      res<-try(read.table(paste0("data/raw/output-cmd2_",site,"_",order,".txt/forceps.",site,".site_",number,"_complete.txt")),silent=T) 
      if (class(res) != "try-error"){# sometimes, the files are empty, and it returns an error message
        colnames(res) <- colnames_res
        int <- temporal_plot(res)
        temp_plot <- temporal_plot_threshold(int) # NB: if error here, don't forget to charge the package "dplyr"
        
        dates <- as.numeric(max(temp_plot$date))- c(900,800,700,600,500,400,300,200,100,0) # years on which we average the biomass
        years_to_keep <- subset(temp_plot,date %in%dates) # we compute a mean value of biomass on those years
        
        biomasses <- aggregate(years_to_keep$biomass, list(years_to_keep$species), mean) # mean per species
        colnames(biomasses) <- c("species","mixture(t/ha)")
        biomasses$'mixture(t/ha)' <- biomasses$'mixture(t/ha)'/(1000*0.08*Nbpatches) # so that the unit becomes t/ha
        biomasses$mixture_relative <- biomasses$'mixture(t/ha)'/sum(biomasses$'mixture(t/ha)') #to have relative biomass
        biomasses$site <- site
        biomasses$order <- order
        biomasses$simul <- number
        BIOMASSES <- rbind(BIOMASSES,biomasses)
      }
    }
  }
  BIOMASSES
}  

total_biomass_final <- function(site){
  # 1. Reads the data frame printed by specific_biomass_final
  # 2. Sums the specific biomasses for each simul
  # 3. Returns a data.frame (e.g. total_biomass_final_Bern.txt)whose colnames are:
  # "site"	"simul"	"decreasing"	"increasing"	"random_1"	"random_10"	"random_2"	"random_3"	etc.
  # which gives the biomasses (t/ha) for each simul in each order in a given site.
  BIOMASSES_sp <- read.table(paste0("data/processed/specific_biomass_final_",site,".txt"),header=T)
  BIOMASSES_tot <- aggregate(BIOMASSES_sp$mixture.t.ha.,list(site = BIOMASSES_sp$site,order = BIOMASSES_sp$order,simul = BIOMASSES_sp$simul),FUN=sum)# sum of the biomasses of the species
  
  biomass_per_order <- spread(BIOMASSES_tot,order,x) # data frame with each order in column
  biomass_per_order
}


confidence_interval_biomass <- function(site){
  # Computes the confidence interval on the table written by tetal_biomass_final
  biomass_per_order <- read.table(paste0("data/processed/total_biomass_final_",site,".txt"),header=T)
  RAND <- biomass_per_order[,5:14]
  int_min <- c()
  int_max <- c()
  mean <- c()
  for (j in c(1:30)){
    int_min <- c(int_min, mean(as.numeric(RAND[j,])) + 1.96 * sd(as.numeric(RAND[j,]))/sqrt(10) )
    int_max <- c(int_max, mean(as.numeric(RAND[j,])) - 1.96 * sd(as.numeric(RAND[j,]))/sqrt(10) )
    mean <- c(mean, mean(as.numeric(RAND[j,])) )
  }
  biomass_per_order$int_min <- int_min
  biomass_per_order$int_max <- int_max
  biomass_per_order$mean <- mean
  
  biomass_per_order
}

# sd(biomass) - removal experiments ####

# Productivity - removal experiments ####

# sd(Productivity) - removal experiments ####

# Former Productivity - removal experiments ####
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

  write.table(RAND,paste0("data/processed/Productivity_species removal experiments_",site,".txt"),sep="\t")
}

