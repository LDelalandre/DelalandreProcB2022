source("R/Common variables.R")
source("R/Analysis_data_2.R")

# Indications de ce qu'il y a à coder
# 1. Ne pas le faire (pour l'instant) sur les random, pour avoir un premier aperçu (mais il suffir juste de changer deux lignes de code)
# 2. Pouvoir utiliser les fonctions qui servent partout --> bien choisir les noms de mes fichiers de sortie !!
# --> notamment, c'est dans la variable "order" que doit se trouver l'info
# par exemple, j'avais "3degrees_increasing" pour l'augmentation de température
# Ici, je peux choisir order= "ag.mono_increasing", par exemple.
# Avec la fonction: productivity_specific(site,order)
# je lis le fichier : read.table(paste0("data/raw/output-cmd2_",site,"_",order,".txt/forceps.",site,".site_",number,"_productivityScene.txt")
# pas cool... sinon, je mets un if dans ma fonction, avec un argument supplémntaire, genre ag.mono = T ou F,
# et s'il est TRUE, je lis complètement autre chose. C'est le plus simple... mais il faut que je remodifie ma fonction et les codes où je l'appelle.

# ATTENTION A NE PAS ECRASER LES FICHIERS PRECEDENTS!! C'EST TRES LONG A FAIRE TOURNER

# "Complete" artificial monocultures ####
# done : GrandeDixence(but random1 and not random_1, etc.), BERN, Bever
sit <- "Bern"

prod_art_mono <- function(sit){
  # reads the monocultures in one site
  # compiles them into an artificial monoculture
  # write a data frame similar to productivity_specific
  distinct_tot <- read.table("data/raw/distinctiveness of the species.txt",header = T)
  distinct_tot$Id <- c(0:29)
  
  ord_decr<-distinct_tot[order(distinct_tot$Di,decreasing=TRUE),]$Id
  ord_incr<-distinct_tot[order(distinct_tot$Di,decreasing=FALSE),]$Id
  
  
  
  # read the simul of the monocultures for each species ranked from i onwards 
  # and compile them in the RESMONO data frame
  dates <- 3949 - c(900,800,700,600,500,400,300,200,100,0) # keep one year every 100 years
  # NB: I need it for temporal stability
  
  RESMONO <- NULL
  for (k in c(1:30)){ # k: every simulation in monoculture (so every one of the 30 species)
    resmono <- try(read.table(paste0("data/raw/output-cmd2_",sit,"_monoculture.txt/forceps.",sit,".site_",k,"_productivityScene.txt"),header=T),silent=T)
    if (class(resmono) != "try-error"){
      colnames(resmono) <- colnames_prod 
      resmono2 <- resmono %>% filter(date %in% dates)
      RESMONO <- rbind(RESMONO,resmono2)
    }
  }
  # RESMONO <- select(RESMONO,date,speciesId,speciesShortName,age,dbh.cm.,height.m.,biomass.kg.)
  RESMONO$totProdBiomass_t_ha <- RESMONO$adultProdBiomass_t_ha + RESMONO$saplingBiomass_t_ha
  
  
  ART.MONO <- data.frame(matrix(0,0,45) )
  colnames(ART.MONO) <- c(colnames_res,c("site","order","simul"))
  ART.MONO <- select(ART.MONO,date,id,speciesId,speciesShortName,age,dbh.cm.,height.m.,biomass.kg.,site,order,simul)
  
  for (j in 31:32){ # go through every order (random 1 to 30, decreasing, and increasing)
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
      res <- try(read.table(paste0("data/raw/output-cmd2_",sit,"_",orde,".txt/forceps.",sit,".site_",i,"_complete.txt"),header=T),silent=T)
      if (class(res) != "try-error"){
        colnames(res) <- colnames_res
        res2 <- res %>% 
          filter(date %in% dates)
      }
      
      nb.trees <- dim(res2)[2] # number of trees in the mixture
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
  ART.MONO
}

# Generate the artificial monocultures ####
for (sit in SITE){
  productivity_specific <- prod_art_mono(sit)
  write.table(productivity_specific,paste0("data/processed/Aggregated monocultures/productivity_specific_",sit,"_art mono.txt"),row.names=F,sep="\t")
}


# Data treatment to plot the artificial monocultures ####

# I have the artificial monoculture in one site
# (it's the usual "productivity_specific data frame)
# I sum the productivities of adults and saplings
# I compute average productivities per species (and per treatment : order, simul) on the 10 years every 100 years 
# Finally, I spread the data frame so that I can plot it

for (sit in SITE){
  ART.MONO <- read.table(paste0("data/processed/Aggregated monocultures/productivity_specific_",sit,"_art mono.txt"),header=T)
  ART.MONO$totProdBiomass_t_ha <- ART.MONO$adultProdBiomass_t_ha + ART.MONO$saplingBiomass_t_ha
  productivities <- ART.MONO %>% 
    group_by(site,order,simul,speciesShortName) %>% 
    summarize(mixture_t_ha=mean(totProdBiomass_t_ha))
  # I could compute relative productivities too
  PROD_sp <- productivities
  PROD_tot <- aggregate(PROD_sp$mixture_t_ha,list(site = PROD_sp$site,order = PROD_sp$order,simul = PROD_sp$simul),FUN=sum)# sum of the productivities of the species
  
  prod_per_order <- spread(PROD_tot,order,x) # data frame with each order in column
  write.table(prod_per_order,paste0("data/processed/Aggregated monocultures/productivity_tot_",sit,"_art mono.txt"),row.names=F,sep="\t")
}


# plot hte result ####

graph_site <- function(prod_per_order){
  ggplot(data=prod_per_order,aes(x=simul))+
    geom_line(aes(y=`decreasing`,colour="#00BFC4"),size=2) +
    geom_line(aes(y=`increasing`,colour="#8766D"),size=2) +
    ggtitle(paste(sit,"Aggredated monoculture")) +
    theme(  plot.title = element_text( size = 17, face = "bold")) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(face="bold", size=22),
          axis.text.y = element_text(face="bold", size=18))+
    xlab(NULL) + ylab(NULL) +
    ggsave(paste0("figures/Aggregated monocultures/Productivity_Aggreg mono_",sit,".png"))
}

for (sit in SITE){
  prod_per_order <- read.table(paste0("data/processed/Aggregated monocultures/productivity_tot_",sit,"_art mono.txt"),header=T)
  graph_site(prod_per_order)
}
  
