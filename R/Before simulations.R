library("dplyr")
library("funrar")

# Chose traits to use to compute the distinctiveness ####
select_traits<-function(traits){ # start from the table traits, and keep the variable below for computing the distinctiveness
  rownames(traits)<-traits[,2]
  traits[,which((colnames(traits) %in% 
                   c(    "S"  ,                       
                         "HMax"    ,                   "AMax"   ,                    "G"  ,                        "DDMin"  ,                   
                         "WiTN"      ,                 "WiTX"       ,                   
                         "DrTol"    ,                  "NTol"       ,                "Brow"         ,              "Ly"      ,                  
                         "La"      ,    "A1max" , "A2"          )) )] 
}



# Compute distinctiveness #####
traits_dist<-function(traits){
  # Computes functional distinctiveness on all the species and adds a new columns with it on the input data frame (i.e. traits, here).
  dist_tot<-compute_dist_matrix(traits,metric='euclidean')
  tidy<-as.data.frame(rownames(traits)) # tidy format for computing distinctiveness in the fonction below
  colnames(tidy)<-"SName"
  distinct_tot<-distinctiveness_com(com_df=tidy,
                                    sp_col=colnames(tidy),abund=NULL,
                                    dist_matrix=dist_tot,relative=F)
  traits$Di<-distinct_tot$Di
  traits
}


traits_dist_gower<-function(traits){ # using both quantitative and qualitative traits --> gower's distance (and not euclidean)
  # Computes functional distinctiveness on all the species and adds a new columns with it on the input data frame (i.e. traits, here).
  dist_tot<-compute_dist_matrix(traits,metric='gower')
  tidy<-as.data.frame(rownames(traits)) # tidy format for computing distinctiveness in the fonction below
  colnames(tidy)<-"SName"
  distinct_tot<-distinctiveness_com(com_df=tidy,
                                    sp_col=colnames(tidy),abund=NULL,
                                    dist_matrix=dist_tot,relative=F)
  traits$Di<-distinct_tot$Di
  traits
}

# Write command files for ForCEEPS ####

# For these 3 functions:
# distinct_tot is a dataframe with the distinctiveness and the name of the species, ordered from 0 to 29.
# length is the number of years of the simulation
# site is a string of characters. e.g. "Bern"
# e.g. Cmd_decr(distACP,length = 2000, yearstobejumped = 999, timestep=100,site="Bern")

Cmd_decr<-function(distinct_tot,length,yearstobejumped,timestep,site){ # distinct_tot is the output of function traits_dist (above)
  # Writes Cmd file with species removed from the most to the least distinctive
  distinct_tot$Id <- c(0:29)
  ord_decr<-distinct_tot[order(distinct_tot$Di,decreasing=TRUE),]$Id
  
  sink(paste0("data/code_ForCEEPS_simulations/Cmd files/cmd2_",site,"_decreasing.txt"))
  cat("# Forceps script command file, format SimulationCommandReader2 (Xavier Morin)")
  cat("\n")
  cat("\n")
  cat(paste0("setupFileName = forceps.setup"))
  cat("\n")
  cat(paste0("numberOfYearsToBeJumped = ",yearstobejumped))
  cat("\n")
  cat(paste0("exportTimeStep = ",timestep))
  cat("\n")
  cat("\n")
  cat("#siteFileName	climateFileName	numberOfYears	potentialSpeciesList")
  cat("\n")
  for (i in 1:length(ord_decr)){
    cat(paste0("forceps.",site,".site	forceps.climate.",site,".txt	",length,"	"),ord_decr[i:length(ord_decr)])
    cat("\n")
  }
  sink()
}

Cmd_incr<-function(distinct_tot,length,yearstobejumped,timestep,site){  
  # Write Cmd file with species removed from the least to the most distinctive
  distinct_tot$Id <- c(0:29)
  ord_incr<-distinct_tot[order(distinct_tot$Di,decreasing=FALSE),]$Id
  
  sink(paste0("data/code_ForCEEPS_simulations/Cmd files/cmd2_",site,"_increasing.txt"))
  cat("# Forceps script command file, format SimulationCommandReader2 (Xavier Morin)")
  cat("\n")
  cat("\n")
  cat(paste0("setupFileName = forceps.setup"))
  cat("\n")
  cat(paste0("numberOfYearsToBeJumped = ",yearstobejumped))
  cat("\n")
  cat(paste0("exportTimeStep = ",timestep))
  cat("\n")
  cat("\n")
  cat("#siteFileName	climateFileName	numberOfYears	potentialSpeciesList")
  cat("\n")
  for (i in 1:length(ord_incr)){
    cat(paste0("forceps.",site,".site	forceps.climate.",site,".txt	",length,"	"),ord_incr[i:length(ord_incr)])
    cat("\n")
  }
  sink()
}

Cmd_rand<-function(distinct_tot,length,yearstobejumped,timestep,site){
  # Write Cmd files with species removed in random order
  # NB: here, we generate 300 dynamics (10 random removals of 30 species)
  for (j in 1:30){
    set.seed(j)
    ord_random<-sample(c(0:29))
    
    sink(paste0("data/code_ForCEEPS_simulations/Cmd files/cmd2_",site,"_random_",j,".txt"))
    cat("# Forceps script command file, format SimulationCommandReader2 (Xavier Morin)")
    cat("\n")
    cat("\n")
    cat(paste0("setupFileName = forceps.setup"))
    cat("\n")
    cat(paste0("numberOfYearsToBeJumped = ",yearstobejumped))
    cat("\n")
    cat(paste0("exportTimeStep = ",timestep))
    cat("\n")
    cat("\n")
    cat("#siteFileName	climateFileName	numberOfYears	potentialSpeciesList")
    cat("\n")
    for (i in 1:length(ord_random)){
      cat(paste0("forceps.",site,".site	forceps.climate.",site,".txt	",length,"	"),ord_random[i:length(ord_random)])
      cat("\n")
    }
    sink()
  }
}

Cmd_mono<-function(distinct_tot,length,yearstobejumped,timestep,site){  
  # Write Cmd file with species removed from the least to the most distinctive
  sink(paste0("data/code_ForCEEPS_simulations/Cmd files/cmd2_",site,"_monoculture.txt"))
  cat("# Forceps script command file, format SimulationCommandReader2 (Xavier Morin)")
  cat("\n")
  cat("\n")
  cat(paste0("setupFileName = forceps.setup"))
  cat("\n")
  cat(paste0("numberOfYearsToBeJumped = ",yearstobejumped))
  cat("\n")
  cat(paste0("exportTimeStep = ",timestep))
  cat("\n")
  cat("\n")
  cat("#siteFileName	climateFileName	numberOfYears	potentialSpeciesList")
  cat("\n")
  for (i in 0:29){
    cat(paste0("forceps.",site,".site	forceps.climate.",site,".txt	",length,"	"),i)
    cat("\n")
  }
  sink()
}

# Id of species in all the orders of removal ####
# It does not give the same random orders as the code I used for my simulations, 
# even though I used "set.seed"... Anyway, it is just another random order.
make_table_name_sname <- function(){
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
  write.table(SPord,"data/processed/Richness gradients orders.txt",row.names=F)
  
  # Correspondence between species Id and species short name
  name_sname <- 
    read.table("data/raw/Traits of the species_complete.txt",header=T) %>% 
    select(Name,SName)
  
  name_sname$Id <- c(0:29) 
  name_sname$SName <- as.character(name_sname$SName)
  write.table(name_sname,"data/processed/correspondence_SName_Id.txt",row.names=F)
}

