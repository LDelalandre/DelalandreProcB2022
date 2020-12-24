source("R/Common variables.R")

specific_biomasses <- function(site){
  # Extract the biomasse of each species in a monoculture and make a data frame with the following columns:
  #  Di       SName      monoculture(t/ha)         Id         monoculture_relative
  
  biomass<-c()
  abundance <- c()
  mean_biom <- c()
  for(i in c(1:30)){
    mean <- read.table(paste0("data/raw/output-cmd2_",site,"_monoculture.txt/forceps.",site,".site_",i,"_mean.txt"))
    # NB it is in the good order because we made the monocultures in the order of the species Id number.
    colnames(mean)<-colnames_mean
    # biomass averaged on 10 years every 100 years
    years_to_keep <- max(mean$date) - c(900,800,700,600,500,400,300,200,100,0)
    meanbiom <- mean(subset(mean,date %in% years_to_keep)$totalBiomass.t.ha.) # average of the last years
    biomass<-c(biomass,meanbiom) 
  }
  
  # specific_values: a data frame with final biomass of each monoculture, etc. 
  specific_values <- read.table("data/raw/distinctiveness of the species.txt",header=T)
  specific_values$'monoculture(t/ha)' <- biomass
  # specific_values$mean_biom <- mean_biom
  specific_values$Id <- c(0:29)
  specific_values$monoculture_relative <- specific_values$'monoculture(t/ha)' / sum(specific_values$'monoculture(t/ha)')
  specific_values
}

add_zeros <- function(site,biomasses,measure){
  sit <- site
  # This function adds lines with species present at the beginning of a simulation but
  # whose biomass is null a the end to the output of specific_biomasses.
  # Useful for the Loreau-Hector analysis.
  data <- biomasses
  data$species <- as.character(data$species)
  data$site <- as.character(data$site)
  data$order <- as.character(data$order)
  # Include species whose biomass or productivity = 0
  # First, have a data frame with the Id of species in all the oders of removal
  distinct_tot <- read.table("data/raw/distinctiveness of the species.txt",header=T)
  distinct_tot$Id <- c(0:29)
  SPord <- data.frame(matrix(data = NA, nrow = 30, ncol = 0))
  SPord$incr <-distinct_tot[order(distinct_tot$Di,decreasing=FALSE),]$Id
  SPord$decr <-distinct_tot[order(distinct_tot$Di,decreasing=TRUE),]$Id
  for (j in 1:30){
    set.seed(j)
    SPord[j+2] <- sample(c(0:29))
  }
  colnames(SPord) <- ORDER
  
  # Second, have the correspondence between species Id and species short name
  name_sname <- read.table("data/Traits of the species_complete.txt",header=T)
  name_sname <- name_sname[,c(1,2)]
  name_sname[,1] <- c(0:29) 
  colnames(name_sname) <- c("Id","SName")
  name_sname$SName <- as.character(name_sname$SName)
  
  # Third, add lines in the file data with 
  for (sim in c(1:30)){# all the simul
    for (ord in ORDER){
      spId <- SPord[,which(colnames(SPord) == ord)] # list of sp Id ordered for a given order of removal
      pr_in <- spId[sim:length(spId)] # list of the Id of species that were present at the beginning of that simulation
      zeros <- name_sname$SName[which(name_sname$Id %in% pr_in)] # Short names of the same species
      sub <- subset(data,order==ord&simul==sim & site==sit)
      to_add <- zeros[which(!(zeros %in% sub$species))] # species that were present at the beginning of the simul, but not at the end.
      for(sp in to_add){
        data <- rbind(data,c(sp,0,0,0,sit,ord,sim))
      }
    }
    
  }
  data
}

biomasses <- function(specif.biomass,specific_val){
  # returns a data frame with colnames:
  # species mixture.t.ha. mixture_relative site      order simul monoculture(t/ha) monoculture_relative       Di
  # Each species that apears is a species that is present at the end of a mixture simulation. Consequently, we don't have all the species.
  # NB: 
  
  # Example for the variables: site="Bern" ; specific_val = specific_values(site) ; number=1 ; order="decreasing" ; Nbpatches = 10 (ou 50)
  # Number is the number of the simulation (1: we didn't remove any species, until 30: there is only one species left)
  biomasss <- specif.biomass
  
  biomasses <- add_zeros(site=site,biomasses=biomasss,measure="biomass_tot")
  # Add a column with absolute biomass of the same species in monoculture
  biom_mono <- c()
  dist <- c()
  ab <- c()
  for (i in 1:dim(biomasses)[1]){
    biom_mono <- c(biom_mono,specific_val[which(specific_val$SName == as.character(biomasses[i,]$species) ),]$'monoculture(t/ha)')
    dist <- c(dist,specific_val[which(specific_val$SName == as.character(biomasses[i,]$species) ),]$Di)
    ab <- c(ab,specific_val[which(specific_val$SName == as.character(biomasses[i,]$species) ),]$abundance)
  }
  biomasses$`monoculture(t/ha)` <- biom_mono
  biomasses$dist <- dist
  biomasses$mono_abundance <- ab
  
  # Add a column with relative biomass of the same species in a mix of monocultures with the same species
  mono_sum <- aggregate(biomasses$`monoculture(t/ha)`, list(order=biomasses$order,simul=biomasses$simul), sum) # mean per species
  biom_relat <- c()
  for (i in 1:dim(biomasses)[1]){
    summed_biomass_mono <- mono_sum[which(mono_sum$order==biomasses[i,]$order & mono_sum$simul==biomasses[i,]$simul),]$x
    biom_relat[i] <- biomasses[i,]$`monoculture(t/ha)`/summed_biomass_mono
  }
  biomasses$monoculture_relative <- biom_relat
  
  # vérification : somme des biomasses relatives = 1. Ok.
  # aggregate(biomasses$monoculture_relative,list(order=biomasses$order,simul=biomasses$simul),FUN=sum)
  # sum(subset(biomasses,simul==1&order=="decreasing")$monoculture_relative)
  biomasses
}


specific_productivities <- function(site){
  # Extract the productivity of each species in a monoculture and make a data.frame with the following columns:
  # Di SName biomass_monoculture     sd_biom Id relative_biomass
  
  # Rq: be careful of the way the folders and files are named.
  # NB: add productivity and sd(productivity) to it when I have the adequate simulations!
  productivity<-c()
  sd<-c()
  mean <- c()
  CV <- c()
  for(number in c(1:30)){
    prod <- try(read.table(paste0("data/raw/output-cmd2_",site,"_monoculture.txt/forceps.",site,".site_",number,"_productivityScene.txt")),silent=T)
    if (class(prod) != "try-error"){# sometimes, the files are empty, and it returns an error message
      colnames(prod)<-colnames_prod
      prod$totProdBiomass_t_ha <- prod$adultProdBiomass_t_ha + prod$saplingBiomass_t_ha
      years_to_keep <- max(prod$date) - c(900,800,700,600,500,400,300,200,100,0)
      prod_to_keep <- subset(prod,date %in% years_to_keep)
      meanpr <- mean(prod_to_keep$totProdBiomass_t_ha)
      productivity<-c(productivity,meanpr)
      sd<-c(sd,sd(prod$totProdBiomass_t_ha)) 
      mean <- c(mean,mean(prod$totProdBiomass_t_ha))
    }
  }
  
  # specific_values: a data frame with final biomass of each monoculture, etc. 
  specific_values <- read.table("data/raw/distinctiveness of the species.txt",header=T)
  specific_values$monoculture <- productivity
  specific_values$sd <- sd
  specific_values$mean <- mean
  specific_values$CV <- sd/mean
  specific_values$Id <- c(0:29)
  specific_values$monoculture_relative <- specific_values$monoculture/ sum(specific_values$monoculture)
  specific_values
}


productivities <- function(site,specific_val){
  # returns a data frame with colnames:
  # species productivity.t.ha. mixture_relative site      order simul monoculture(t/ha) monoculture_relative       Di
  # Each species that apears is a species that is present at the end of a mixture simulation. Consequently, we don't have all the species.
  
  # Example for the variables: site="Bern" ; specific_val = specific_values(site) ; number=1 ; order="decreasing" ; Nbpatches = 10 (ou 50)
  # Number is the number of the simulation (1: we didn't remove any species, until 30: there is only one species left)
  productivityy <- read.table(paste0("data/processed/productivity_specific_",site,".txt"),header=T)
  productivity <- add_zeros(site=site,biomasses=productivityy,measure="productivity_tot")
  # Add a column with absolute biomass of the same species in monoculture
  biom_mono <- c()
  dist <- c()
  for (i in 1:dim(productivity)[1]){
    biom_mono <- c(biom_mono,specific_val[which(specific_val$SName == as.character(productivity[i,]$species) ),]$monoculture)
    dist <- c(dist,specific_val[which(specific_val$SName == as.character(productivity[i,]$species) ),]$Di)
  }
  productivity$monoculture <- biom_mono
  productivity$dist <- dist
  
  # Add a column with relative productivity of the same species in a mix of monocultures with the same species
  mono_sum <- aggregate(productivity$monoculture, list(order=productivity$order,simul=productivity$simul), sum) # mean per species
  biom_relat <- c()
  for (i in 1:dim(productivity)[1]){
    summed_biomass_mono <- mono_sum[which(mono_sum$order==productivity[i,]$order & mono_sum$simul==productivity[i,]$simul),]$x
    biom_relat[i] <- productivity[i,]$monoculture/summed_biomass_mono
  }
  productivity$monoculture_relative <- biom_relat
  
  # vérification : somme des biomasses relatives = 1. Ok.
  # aggregate(biomasses$monoculture_relative,list(order=biomasses$order,simul=biomasses$simul),FUN=sum)
  # sum(subset(biomasses,simul==1&order=="decreasing")$monoculture_relative)
  productivity
}
