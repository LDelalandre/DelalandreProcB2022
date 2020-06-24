source("R/Common variables.R")
for (site in SITE){
  # species ordered according to their biomass
  biomass <- read.table(paste0("data/processed/biomass_specific_",site,".txt"),header=T)
  biosub <- subset(biomass,order=="decreasing" & simul==1)
  bio <- biosub[order(biosub$mixture.t.ha.,decreasing=T),]$species
  
  # correspondence between species name and species Id
  traits <- read.table("data/Traits of the species_complete.txt",header=T)
  traits$Id <- c(0:29)
  corresp <- select(traits,Id,Name,SName)
  order <- c()
  for (i in 1:length(bio)){
    order <- c(order,which(as.character(corresp$SName)==as.character(bio[i])))
  }
  order <- order-1 # because id = position in corresponding data frame - 1
  absent <- which(!(c(0:29)%in% order))-1 # species that were absent from the mixture
  order2 <- c(order,absent) # I add them to star the simulations with 30 species so that Ican compare it with the rest.
  
  Cmd_given_order(order2,length, yearstobejumped, timestep,site)
}


Cmd_given_order<-function(order2,length,yearstobejumped,timestep,site){ # distinct_tot is the output of function traits_dist (above)
 # order must be a vector with the Id of the species, and the ones we want to remove first  must be placed first in this vector
 # Writes Cmd file with species removed from the most to the least distinctive
  ord_decr<-order2
  
  sink(paste0("data/raw/cmd2_",site,"_abundance.txt"))
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

# Read the results of the simulations ####
site <- "GrandeDixence" # Bern ou Sion
# Au final, il me suffira de :
# ne pas toucher à la fonction productivity_specific, mais ajouter "abundance" comme élément du vecteur ORDER


#productivity_specific <- function(site){
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
# }


  PROD_tot <- aggregate(PROD_sp$mixture_t_ha,list(site = PROD_sp$site,order = PROD_sp$order,simul = PROD_sp$simul),FUN=sum)# sum of the productivities of the species
  
  prod_per_order <- spread(PROD_tot,order,x) # data frame with each order in column
  prod_per_order

  measure <- "productivity_tot"
  result <- read.table(paste0("data/processed/",measure,"_",site,"_with interval.txt"),header=T)
  abundance <- prod_per_order$abundance
  result <- cbind(result,abundance)
  
  ggplot(result,aes(x=simul-1,y=decreasing,color="Removing distinct species first")) +
    labs(x="Number of species removed",y=measure) +
    geom_line()+
    geom_ribbon(aes(ymin=int_min, ymax=int_max),fill="grey60", alpha=0.5,colour="black") +
    geom_line(aes(x=simul-1,y=increasing, color="Removing distinct species last")) +
    geom_line(aes(x=simul-1,y=abundance, color="Removing abundant species first")) +
    theme(legend.position = "bottom") +
    ggtitle(site) +
    scale_x_continuous(breaks = 2*c(1:15)) +
    ggsave(paste0("figures/productivity_abundance_",site,".png"))

  
  # petit check
  
BernAb <- read.table("data/raw/output-cmd2_Bern_abundance.txt/forceps.Bern.site_1_productivityScene.txt")
colnames(BernAb) <- colnames_prod
sum(subset(BernAb,date==3949)$adultProdBiomass)

for (i in c(1:30)){
  BernDecr <- read.table(paste0("data/raw/output-cmd2_Bern_random_",i,".txt/forceps.Bern.site_1_productivityScene.txt"))
  colnames(BernDecr) <- colnames_prod
  print(sum(subset(BernDecr,date==3949)$adultProdBiomass))
}


