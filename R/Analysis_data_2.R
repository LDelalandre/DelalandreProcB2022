source("R/Common variables.R")

# Biomass - species removal #####
biomass_specific_aggreg_mono <- function(site) {
  # 1. Takes the output of the simulations (complete files).
  # 2. Selects biomass above threshold
  # 3. Selects a subset of ten years and averages the biomass of each species on these years
  # 4. Returns data frame (e.g. specific_biomass_final_Bern.txt) whose columns are:
  # "species" "abundance"	"mixture(t/ha)"	"mixture_relative"	"site"	"order"	"simul"
  RES <- read.table(paste0("data/processed/Aggregated monocultures_method 2/artificial monocultures_",sit,".txt"),header=T)
  
  BIOMASSES <- as.data.frame(matrix(nrow=0,ncol=4,dimnames = list(NULL,c("species", "mixture(t/ha)", "mixture_relative", "simul"))))
  for (order in ORDER){
    orde <- order
    biomass_inc<-c()
    sd_biom_inc<-c()
    for(number in c(1:30)){
      res <- subset(RES,order==orde & simul==number)
      # res<-try(read.table(paste0("data/raw/output-cmd2_",site,"_",order,".txt/forceps.",site,".site_",number,"_complete.txt")),silent=T) 
      if (class(res) != "try-error"){# sometimes, the files are empty, and it returns an error message
        # colnames(res) <- colnames_res
        int <- temporal_plot(res)
        temp_plot <- temporal_plot_threshold(int) # NB: if error here, don't forget to charge the package "dplyr"
        
        dates <- as.numeric(max(temp_plot$date))- c(900,800,700,600,500,400,300,200,100,0) # years on which we average the biomass
        years_to_keep <- subset(temp_plot,date %in%dates) # we compute a mean value of biomass on those years
        
        if (dim(years_to_keep)[1]>0){
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
    }
  }
  BIOMASSES
}  

biomass_tot_aggreg_mono <- function(site){
  # 1. Reads the data frame printed by specific_biomass_final
  # 2. Sums the specific biomasses for each simul
  # 3. Returns a data.frame (e.g. total_biomass_final_Bern.txt)whose colnames are:
  # "site"	"simul"	"decreasing"	"increasing"	"random_1"	"random_10"	"random_2"	"random_3"	etc.
  # which gives the biomasses (t/ha) for each simul in each order in a given site.
  BIOMASSES_sp <- read.table(paste0("data/processed/Aggregated monocultures_method 2/biomass_specific_",site,"_aggreg_mono.txt"),header=T)
  BIOMASSES_tot <- aggregate(BIOMASSES_sp$mixture.t.ha.,list(site = BIOMASSES_sp$site,order = BIOMASSES_sp$order,simul = BIOMASSES_sp$simul),FUN=sum)# sum of the biomasses of the species
  
  biomass_per_order <- spread(BIOMASSES_tot,order,x) # data frame with each order in column
  biomass_per_order
}

# sd(biomass) - removal experiments ####
sd_biomass_specific <- function(site){
  SIGMA <- NULL
  for (order in ORDER){
    for(number in c(1:30)){
      res<-try(read.table(paste0("data/raw/output-cmd2_",site,"_",order,".txt/forceps.",site,".site_",number,"_complete.txt")),silent=T) 
      if (class(res) != "try-error"){# sometimes, the files are empty, and it returns an error message
        colnames(res) <- colnames_res
        int <- temporal_plot(res)
        temp_plot <- temporal_plot_threshold(int) # NB: if error here, don't forget to charge the package "dplyr"
        temp_plot$biomass <- temp_plot$biomass/(1000*0.08*Nbpatches) # so that unit is t/ha
        sigma <- aggregate(temp_plot$biomass, list(temp_plot$species), sd)
        colnames(sigma) <- c("species","sd")
        sigma$mean <- aggregate(temp_plot$biomass, list(temp_plot$species), mean)$x
        sigma$CV <- sigma$sd/sigma$mean
        sigma$site <- site
        sigma$order <- order
        sigma$simul <- number

        SIGMA <- rbind(SIGMA,sigma)
      }
    }
  }
  SIGMA
}

sd_biomass_tot <- function(site){ # Takes a looong time
  SIGMA <- NULL
  for (order in ORDER){
    for(number in c(1:30)){
      res<-try(read.table(paste0("data/raw/output-cmd2_",site,"_",order,".txt/forceps.",site,".site_",number,"_complete.txt")),silent=T) 
      if (class(res) != "try-error"){# sometimes, the files are empty, and it returns an error message
        colnames(res) <- colnames_res
        int <- temporal_plot(res)
        temp_plot <- temporal_plot_threshold(int) # NB: if error here, don't forget to charge the package "dplyr"
        temp_plot$biomass <- temp_plot$biomass/(1000*0.08*Nbpatches) # so that unit is t/ha
        summed <- aggregate(temp_plot$biomass,list(temp_plot$date),sum)
        sigma <- sd(summed$x)
        mu <- mean(summed$x)
        CV <- sigma/mu
        sigma <- data.frame(sd=sigma,mean=mu,CV,site,order,simul=number)
        
        SIGMA <- rbind(SIGMA,sigma)
      }
    }
  }
  SIGMA
}


# Productivity - removal experiments ####
productivity_specific <- function(site){
  PROD <- NULL
  for (order in ORDER){
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
sd_productivity_specific <- function(site){
  SIGMA <- NULL
  for (order in ORDER){
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
  }
  SIGMA
}

sd_productivity_tot <- function(site){
  SIGMA <- NULL
  for (order in ORDER){
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
  }
  SIGMA
}


# Confidence interval ####
confidence_interval <- function(site,measure){
  # Computes the confidence interval on the table written by total_biomass_final
  # measure is either "biomass" or "productivity"
  data <- read.table(paste0("data/processed/",measure,"_",site,".txt"),header=T)
  RAND <- select(data,ORDER[3:12])# select only the random results

  int_min <- c()
  int_max <- c()
  mean <- c()
  for (j in c(1:30)){
    int_min <- c(int_min, mean(as.numeric(RAND[j,]),na.rm=TRUE) + 1.96 * sd(as.numeric(RAND[j,]),na.rm=TRUE)/sqrt(10) )
    int_max <- c(int_max, mean(as.numeric(RAND[j,]),na.rm=TRUE) - 1.96 * sd(as.numeric(RAND[j,]),na.rm=TRUE)/sqrt(10) )
    mean <- c(mean, mean(as.numeric(RAND[j,]),na.rm=TRUE) )
  }
  data$int_min <- int_min
  data$int_max <- int_max
  data$mean <- mean
  
  data
}


