source("R/Common variables.R")

# find if a species name is present in a vector of species names
persist <- function(sp,persistent.sp){
  # returns true if and only if species sp is present in the vector persistent.sp
  any(grepl(sp,persistent.sp))
}
# I should use this function in the BIOMASS section (in "2. Process data_2.R"): instead of filter,
# I should add a column indicating if a species persists or not.
# Then, when using this file later, I should filter species which persist: filter (persists==T)

# Monocultures ####

specific_biomasses_mono <- function(site){
  # Extract the biomasse of each species in a monoculture and make a data frame with the following columns:
  #  Di       SName      monoculture(t/ha)         Id         monoculture_relative
  sp_id <- read.table("data/processed/correspondence_SName_Id.txt",header=T)
  biomass<-c()
  Id <- c()
  SName <- c()
  abundance <- c()
  mean_biom <- c()
  for(i in c(1:30)){
    mean <- read.table(paste0("data/raw/Output_ForCEEPS/",site,"/output-cmd2_",site,"_monoculture.txt/forceps.",site,".site_",i,"_mean.txt"))
    # NB it is in the good order because we made the monocultures in the order of the species Id number.
    colnames(mean)<-colnames_mean
    # biomass averaged on 10 years every 100 years
    years_to_keep <- max(mean$date) - c(900,800,700,600,500,400,300,200,100,0)
    meanbiom <- mean(subset(mean,date %in% years_to_keep)$totalBiomass.t.ha.) # average of the last years
    biomass<-c(biomass,meanbiom) 
    SName <- c(SName,sp_id[i,]$SName)
    Id <- c(Id,sp_id[i,]$Id)
  }
  
  sp_biomass <- data.frame('monoculture(t/ha)'= biomass,SName=SName,Id)
  # specific_values: a data frame with final biomass of each monoculture, etc. 
  sp_di <- read.table("data/raw/distinctiveness of the species.txt",header=T)
  specific_values <- merge(sp_biomass,sp_di,by="SName")
  # specific_values$mean_biom <- mean_biom
  specific_values$monoculture_relative <- specific_values$monoculture.t.ha. / sum(specific_values$monoculture.t.ha.)
  specific_values
}


specific_productivities_mono <- function(site){
  # Extract the productivity of each species in a monoculture and make a data.frame with the following columns:
  # Di SName biomass_monoculture     sd_biom Id relative_biomass
  
  # Rq: be careful of the way the folders and files are named.
  # NB: add productivity and sd(productivity) to it when I have the adequate simulations!
  productivity<-c()
  SName <- c()
  # sd<-c()
  # mean <- c()
  for(number in c(1:30)){
    prod <- try(read.table(paste0("data/raw/Output_ForCEEPS/",site,"/output-cmd2_",site,"_monoculture.txt/forceps.",site,".site_",number,"_productivityScene.txt")),silent=T)
    if (class(prod) != "try-error"){# sometimes, the files are empty, and it returns an error message
      colnames(prod)<-colnames_prod
      prod$totProdBiomass_t_ha <- prod$adultProdBiomass_t_ha + prod$saplingBiomass_t_ha
      years_to_keep <- max(prod$date) - c(900,800,700,600,500,400,300,200,100,0)
      prod_to_keep <- subset(prod,date %in% years_to_keep)
      meanpr <- mean(prod_to_keep$totProdBiomass_t_ha)
      productivity<-c(productivity,meanpr)
      # sd<-c(sd,sd(prod$totProdBiomass_t_ha)) 
      # mean <- c(mean,mean(prod$totProdBiomass_t_ha))
      SName <- c(SName,unique(prod$speciesShortName))
    }
  }
  sp_prod <- data.frame(productivity,SName)
  
  # specific_values: a data frame with final biomass of each monoculture, etc. 
  sp_di <- read.table("data/raw/distinctiveness of the species.txt",header=T)
  specific_values <- merge(sp_di,sp_prod,by="SName") %>% 
    rename(monoculture=productivity)
  sp_id <- read.table("data/processed/correspondence_SName_Id.txt",header=T)
  specific_values <- merge(specific_values,sp_id,by="SName")
  specific_values$monoculture_relative <- specific_values$monoculture/ sum(specific_values$monoculture)
  specific_values
}


# Productivity - removal experiments ####
productivity_specific <- function(site,order,number){
    PROD <- NULL
    
      prod <- try(read.table(paste0("data/raw/Output_ForCEEPS/",site,"/output-cmd2_",site,"_",order,".txt/forceps.",site,".site_",number,"_productivityScene.txt")),silent=T)
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


