source("R/Common variables.R")

# Productivity - removal experiments ####
productivity_specific <- function(site,order,number){
    PROD <- NULL
    
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
  PROD
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


