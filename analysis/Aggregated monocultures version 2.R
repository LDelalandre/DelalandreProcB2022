source("R/Common variables.R")
source("R/Analysis_data_2.R")
# ATTENTION A NE PAS ECRASER LES FICHIERS PRECEDENTS!! C'EST TRES LONG A FAIRE TOURNER

# "Complete" artificial monocultures ####
# done : GrandeDixence(but random1 and not random_1, etc.), BERN, Bever
sit <- "Bern"

distinct_tot <- read.table("data/raw/distinctiveness of the species.txt",header = T)
distinct_tot$Id <- c(0:29)

ord_decr<-distinct_tot[order(distinct_tot$Di,decreasing=TRUE),]$Id
ord_incr<-distinct_tot[order(distinct_tot$Di,decreasing=FALSE),]$Id



# read the simul of the monocultures for each species ranked from i onwards 
# and compile them in the RESMONO data frame
dates <- 3950 - c(900,800,700,600,500,400,300,200,100,0) # keep one year every 100 years
# NB: I need it for temporal stability

RESMONO <- NULL
for (k in c(1:30)){ # k: every simulation in monoculture (so every one of the 30 species)
  resmono <- try(read.table(paste0("data/raw/output-cmd2_",sit,"_monoculture.txt/forceps.",sit,".site_",k,"_complete.txt"),header=T),silent=T)
  if (class(resmono) != "try-error"){
    colnames(resmono) <- colnames_res 
    resmono2 <- resmono %>% filter(date %in% dates)
        RESMONO <- rbind(RESMONO,resmono2)
  }
}
RESMONO <- select(RESMONO,date,id,speciesId,speciesShortName,age,dbh.cm.,height.m.,biomass.kg.)

ART.MONO <- data.frame(matrix(0,0,45) )
colnames(ART.MONO) <- c(colnames_res,c("site","order","simul"))
ART.MONO <- select(ART.MONO,date,id,speciesId,speciesShortName,age,dbh.cm.,height.m.,biomass.kg.,site,order,simul)

for (j in 1:32){ # go through every order (random 1 to 30, decreasing, and increasing)
  print("order") ; print(j)
  if (j < 31){
    set.seed(j)
    ord_random<-sample(c(0:29))
    ord <- ord_random
    orde <- paste0("random_",j)
  } else if (j == 31){
    ord <- ord_decr
    orde <- "decreasing"
  } else if (j==32){
    ord <- ord_incr
    orde <- "increasing"
  }

  for (i in c(1:30)){   # i is the number of the simulation
    print("simulation") ; print(i)
    species_kept <- ord[i:30]
    RESMONO2 <- filter(RESMONO,speciesId %in%species_kept)
    
    # read the results of the simul in which the species ranked before i have been removed
    res <- try(read.table(paste0("data/raw/output-cmd2_",sit,"_",order,".txt/forceps.",sit,".site_",i,"_complete.txt"),header=T),silent=T)
    if (class(res) != "try-error"){
      colnames(res) <- colnames_res
      res2 <- res %>% 
        filter(date %in% dates)
    }

    nb.trees <- dim(res2)[1] # number of trees in the mixture
    if(nb.trees<dim(RESMONO2)[1]){
      art.mono <- sample_n(RESMONO2,size=nb.trees,replace=F) # artificial monoculture with the same nb of individuals as in the mixture
    } else {
      art.mono <- RESMONO2 # CAREFUL: here, if there are less individuals in the pooled monocultures than
      # in the mixture, the aggregated monoculture that I will use for the rest of the analysis will have less individuals than the mixture
    }
    
    if (dim(art.mono)[1] > 0){
      art.mono$site <- sit
      art.mono$order <- orde
      art.mono$simul <- i
    }
    
    
    ART.MONO <- rbind(ART.MONO,art.mono)
  }
}

# write.table(ART.MONO,paste0("data/processed/Aggregated monocultures_method 2/artificial monocultures_",sit,".txt"),row.names=F)

biomass <- biomass_specific_aggreg_mono(sit)
write.table(biomass,paste0("data/processed/Aggregated monocultures_method 2/biomass_specific_",sit,"_aggreg_mono.txt"),row.names=F)

biomtot <- biomass_tot_aggreg_mono(sit)
write.table(biomtot,paste0("data/processed/Aggregated monocultures_method 2/biomass_tot_",sit,"_aggreg_mono.txt"),row.names=F)

# add confidence interval ####
# for (sit in SITE){
  table <- read.table(paste0("data/processed/Aggregated monocultures_method 2/biomass_tot_",sit,"_aggreg_mono.txt"),header=T)
  table3 <- median_conf_int(table)
  
  write.table(table3,paste0("data/processed/Aggregated monocultures_method 2/biomass_tot_",sit,"_aggreg_mono_with_interval_median.txt"))
# }

# Plot biomass erosion scenarios ####

result <- read.table(paste0("data/processed/Aggregated monocultures_method 2/biomass_tot_",sit,"_aggreg_mono_with_interval_median.txt"),header=T)
  
ggplot(result,aes(x=simul-1,y=decreasing,color="Removing distinct species first")) +
    labs(x="Number of species removed",y="total biomass") +
    # geom_line(size=1)+
    geom_ribbon(aes(ymin=int_min, ymax=int_max),fill="grey60", alpha=0.5,colour="black") +
    geom_line(aes(x=simul-1,y=increasing, color="#8766D"),size=2) +
    geom_line(aes(x=simul-1,y=decreasing, color="#00BFC4"),size=2) +
    theme(legend.position = "bottom") +
    ggtitle(sit) +
    scale_x_continuous(breaks = 2*c(1:15)) +
    theme(legend.title = element_blank()) +
    # theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none" ) # virer tous les titres
    ggsave(paste0("figures/Aggregated monocultures/Biomass_",sit,".png"))
  



    