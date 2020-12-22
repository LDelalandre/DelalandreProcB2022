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
biomass_specific <- function(site,order) {
  # 1. Takes the output of the simulations (complete files).
  # 2. Selects biomass above threshold
  # 3. Selects a subset of ten years and averages the biomass of each species on these years
  # 4. Returns data frame (e.g. specific_biomass_final_Bern.txt) whose columns are:
  # "species" "abundance"	"mixture(t/ha)"	"mixture_relative"	"site"	"order"	"simul"
  BIOMASSES <- as.data.frame(matrix(nrow=0,ncol=4,dimnames = list(NULL,c("species", "mixture(t/ha)", "mixture_relative", "simul"))))
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
        abundances <-  aggregate(years_to_keep$abundance, list(years_to_keep$species), mean)
        colnames(biomasses) <- c("species","mixture(t/ha)")
        biomasses$abundance <-  abundances$x
        biomasses$'mixture(t/ha)' <- biomasses$'mixture(t/ha)'/(1000*0.08*Nbpatches) # so that the unit becomes t/ha
        biomasses$mixture_relative <- biomasses$'mixture(t/ha)'/sum(biomasses$'mixture(t/ha)') #to have relative biomass
        biomasses$site <- site
        biomasses$order <- order
        biomasses$simul <- number
        BIOMASSES <- rbind(BIOMASSES,biomasses)
      }
    }
  BIOMASSES
}  

biomass_tot <- function(site){
  # 1. Reads the data frame printed by specific_biomass_final
  # 2. Sums the specific biomasses for each simul
  # 3. Returns a data.frame (e.g. total_biomass_final_Bern.txt)whose colnames are:
  # "site"	"simul"	"decreasing"	"increasing"	"random_1"	"random_10"	"random_2"	"random_3"	etc.
  # which gives the biomasses (t/ha) for each simul in each order in a given site.
  BIOMASSES_sp <- read.table(paste0("data/processed/biomass_specific_",site,".txt"),header=T)
  BIOMASSES_tot <- aggregate(BIOMASSES_sp$mixture.t.ha.,list(site = BIOMASSES_sp$site,order = BIOMASSES_sp$order,simul = BIOMASSES_sp$simul),FUN=sum)# sum of the biomasses of the species
  
  biomass_per_order <- spread(BIOMASSES_tot,order,x) # data frame with each order in column
  biomass_per_order
}


# Productivity - removal experiments ####
productivity_specific <- function(site,order){
    PROD <- NULL
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
  PROD
}

productivity_total <- function(site){
  # 1. Reads the data frame printed by productivity_specific
  # 2. Sums the specific biomasses for each simul
  # 3. Returns a data.frame (e.g. total_biomass_final_Bern.txt)whose colnames are:
  # "site"	"simul"	"decreasing"	"increasing"	"random_1"	"random_10"	"random_2"	"random_3"	etc.
  # which gives the biomasses (t/ha) for each simul in each order in a given site.
  PROD_sp <- read.table(paste0("data/processed/productivity_specific_",site,".txt"),header=T)
  PROD_tot <- aggregate(PROD_sp$mixture_t_ha,list(site = PROD_sp$site,order = PROD_sp$order,simul = PROD_sp$simul),FUN=sum)# sum of the productivities of the species
  
  prod_per_order <- spread(PROD_tot,order,x) # data frame with each order in column
  prod_per_order
}

# sd(Productivity) - removal experiments ####
sd_productivity_specific <- function(site,order){
  SIGMA <- NULL
    for(number in c(1:30)){
      prod <- try(read.table(paste0("data/raw/output-cmd2_",site,"_",order,".txt/forceps.",site,".site_",number,"_productivityScene.txt")),silent=T)
      if (class(prod) != "try-error"){# sometimes, the files are empty, and it returns an error message
        colnames(prod)<-colnames_prod
        prod$totProdBiomass_t_ha <- prod$adultProdBiomass_t_ha + prod$saplingBiomass_t_ha
        dates <- as.numeric(max(prod$date))- c(900,800,700,600,500,400,300,200,100,0)
        prod <- subset(prod,date%in%dates) # keep 10 years every 100-year

        sigma <- aggregate(prod$totProdBiomass_t_ha, list(prod$speciesShortName), sd)
        colnames(sigma) <- c("species","sd")
        sigma$mean <- aggregate(prod$totProdBiomass_t_ha, list(prod$speciesShortName), mean)$x
        sigma$TS <- sigma$mean/sigma$sd
        sigma$site <- site
        sigma$order <- order
        sigma$simul <- number
        
        SIGMA <- rbind(SIGMA,sigma)
      }
    }
  SIGMA
}

sd_productivity_tot <- function(site){
  # independent from sd_productivity_specific
  SIGMA <- NULL
    for(number in c(1:30)){
      prod <- try(read.table(paste0("data/raw/output-cmd2_",site,"_",order,".txt/forceps.",site,".site_",number,"_productivityScene.txt")),silent=T)
      if (class(prod) != "try-error"){# sometimes, the files are empty, and it returns an error message
        colnames(prod)<-colnames_prod
        prod$totProdBiomass_t_ha <- prod$adultProdBiomass_t_ha + prod$saplingBiomass_t_ha
        dates <- as.numeric(max(prod$date))- c(900,800,700,600,500,400,300,200,100,0)
        prod <- subset(prod,date%in%dates) # keep 10 years every 100-year
        
        summed <- aggregate(prod$totProdBiomass_t_ha, list(prod$date), sum)
        sigma <- sd(summed$x)
        mu <- mean(summed$x)
        TS <- mu/sigma
        sigma <- data.frame(sd=sigma,mean=mu,TS,site,order,simul=number)
        
        SIGMA <- rbind(SIGMA,sigma)
      }
    }
  SIGMA
}


# Median confidence interval ####
# pbinom(9, size=30, prob=.5) # On cherche k (ici k=9) tel que la probabilité d'avoir moins de k succès soit de 0.025
# qbinom(.025, size=30, prob=.5) # en fait il faut chercher dans ce sens. On cherche le 2.5ème quantile d'une binomiale avec 30 tirages et une proba de succès de 0.5

median_conf_int <- function(table){
  # table is a data frame with in column :  site simul decreasing increasing random_1, etc.
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
  
  table
}


