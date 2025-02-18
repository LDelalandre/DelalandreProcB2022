library(tidyverse)
source("R/comp_fct_dist.R")
source("R/Analysis_data.R")
traits <- read.table("data/raw/Traits of the species_complete.txt",header=T)
rownames(traits) <- NULL
ftraits <- traits %>% 
  select(SName,HMax,La,Ly,G) %>% 
  column_to_rownames("SName")
sname_id <- read.table("data/raw/correspondence_SName_Id.txt",header=T)

colnames_prod <- colnames(read.table("data/raw/colnames_productivityScene.txt",header=T))
colnames_prod2 <- colnames(read.table("data/raw/colnames_productivity.txt",header=T))

biomass_random <- read.csv2("data/processed/biomass_random_pour_cyrille.csv")


sit <- "Bever"
Nbpatches <- 50
length <- 2000
yearstobejumped <- 999
timestep <- 100

pers <- biomass_random %>% 
  filter(site == sit & simul == 1) %>% 
  group_by(order) %>% 
  filter(persists==T) %>% 
  summarize(n=n()) %>% 
  arrange(desc(n))
# 14 à 18 espèces persistent en mélange de 30 espèces à Bern.
# 11 espèces persistent en mélange de 30 sp à Bever


# Sélection des plus distinctes sur celles qui persistent en mélange à 30 espèces
persistent_sp <- biomass_random %>% 
  filter(site ==sit & simul == 1) %>% 
  group_by(order) %>% 
  filter(persists==T) %>% 
  # filter(order=="random_1") %>% # ça doit être pareil quel que soit l'ordre : je prends que la simul à 30 espèces initialement. Non : il y a de la variabilité !
  pull(speciesShortName) %>% 
  unique() # je prends les espèces qui persistent des fois dans les mélanges à 30 espèces. (on passe à 20 espèces)

# Generate Cmd files ####
df_di <- ftraits %>% 
  filter(rownames(.) %in% persistent_sp) %>% 
  comp_fct_dist() %>% # ATTENTION il faut aussi calculer la Di autrement qu'en faisant une ACP (direct sur les traits). Est-ce que ça change beaucoup ?
  merge(sname_id,by="SName") %>% 
  arrange(desc(Di))

ord <- df_di$Id
local_sp_pool <- df_di %>% arrange(Id) %>% pull(Id)

sink(paste0("analysis_subset_traits/data/ForCEEPS/cmd_files_temporary/cmd2_",sit,"_subset_traits_2sp.txt"))
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

# Retrait dans l'ordre décroissant de Di
# for (i in 1:length(ord)){
#   cat(paste0("forceps.",sit,".site	forceps.climate.",sit,".txt	",length,"	"),ord[i:length(ord)])
#   cat("\n")
# }

pairs_removed <- utils::combn(ord,2) # toutes les combinaisons de trois espèces parmi le pool local 

# Toutes les combinaison de toutes les espèces du pool local moins deux espèces
for (i in 1:dim(pairs_removed)[2] ){
  pair_removed <- pairs_removed[,i]
  cat(paste0("forceps.",sit,".site	forceps.climate.",sit,".txt	",length,"	"),
      local_sp_pool[-which(local_sp_pool%in% pair_removed)] )
  cat("\n")
}

sink()


# Process data ####
# Il faut que je sorte les infos pour les mélanges avec mes nouvelles simul
# Et que je voie qui persiste (en biomasse), même si c'est peut-être pas
# nécessaire pour la productivité ?


# Productivity - removal experiments ####

productivity_specific_adapted <- function(file,site,order,number){
  PROD <- NULL
  
  prod <- try(read.table(file ),silent=T)
  if (class(prod) != "try-error"){# sometimes, the files are empty, and it returns an error message
    colnames(prod)<-colnames_prod
    years_to_keep <- max(prod$date) - c(900,800,700,600,500,400,300,200,100,0)
    prod_to_keep <- subset(prod,date %in% years_to_keep)
    prod_to_keep$totProdBiomass_t_ha <- prod_to_keep$adultProdBiomass_t_ha + prod_to_keep$saplingBiomass_t_ha
    
    # specific productivities
    # productivities <- aggregate(prod_to_keep$totProdBiomass_t_ha, list(prod_to_keep$speciesShortName), mean)
    # colnames(productivities) <- c("species","mixture_t_ha")
    
    # idem with dplyr:
    productivities <- prod_to_keep %>%
      rename(SName = speciesShortName) %>% 
      group_by(SName) %>% 
      summarize(mixture_t_ha = mean(totProdBiomass_t_ha)) 
    
    productivities$site <- site
    productivities$order <- order
    productivities$simul <- number
    productivities$mixture_relative <- productivities$mixture_t_ha/sum(productivities$mixture_t_ha)
    PROD <- rbind(PROD,productivities)
    
  }
  PROD
}

file2sp <- paste0("analysis_subset_traits/data/raw/Output_ForCEEPS/",
                  sit,
                  "/output-cmd2_",site,"_subset_traits_",order,".txt/forceps.",
                  sit,".site_1_productivityScene.txt")
productivity_specific_adapted(file2sp,sit,order,number = 1)

file2sp <- paste0("analysis_subset_traits/data/raw/Output_ForCEEPS/",sit,
                  "/output-cmd2_",site,"_subset_traits_",order,".txt/forceps.",
                  sit,".site_105_productivityScene.txt")
productivity_specific_adapted(file2sp,sit,order,number = 1)


# pairs removed ####
pairs_removed <- utils::combn(ord,2) # toutes les combinaisons de deux espèces parmi le pool local 
pairs_removed_2 <- pairs_removed %>% 
  t() %>% 
  as.data.frame()
pairs_removed_2$number = c(1:dim(pairs_removed_2)[1])

nb_simul <- dim(pairs_removed_2)[1]
# attention, pour l'instant j'ai 150 simul et pas 190 car j'avais pris moins d'espèces
nb_simul <- 105 # temporaire

#_______________________________________________________________________________
# Read biomass data ####
ALL <- NULL
for (number in c(1:nb_simul)){
  file_biom2sp <- file2sp <- paste0("analysis_subset_traits/data/raw/Output_ForCEEPS/",sit,
                                    "/output-cmd2_",site,"_subset_traits_",order,".txt/forceps.",sit,".site_",numb,"_complete.txt")
  res<-try(read.table(file_biom2sp),silent=T) 
  if (class(res) != "try-error"){# sometimes, the files are empty, and it returns an error message
    colnames(res) <- colnames_res
    # gives mean biomass (in t/ha) of each species on 10 years sampled every 100 years
    # (done for one species pool (simul = 1 : all the species, simul = 30 : 1 sp remaining), for one
    # order of removal, and in one site)
    res2 <- res %>% 
      group_by(date,speciesShortName) %>% 
      select(date,speciesShortName,biomass.kg.) %>% 
      summarise(biom_tot=sum(biomass.kg.)) %>% # sum of the biomass per species per date
      mutate(biom_tot=biom_tot/(1000*0.08*Nbpatches)) %>% # so that the unit becomes t/ha
      mutate(site=site,order=order,simul=number) 
    
    # Apply a biomass threshold and keep species whose biomass in the last year is above a given fraction (the variable "threshold")
    # of the total biomass of the community.
    # Filtered contains the filtered species:
    filtered <- res2 %>% 
      group_by(date) %>% 
      filter(date==3950) %>% 
      mutate(sum=sum(biom_tot)) %>% # sum of total biomass used to compare with the biomass of one species
      filter(biom_tot>threshold*sum) %>% 
      pull(speciesShortName)
    
    
    # Average species biomass on 10 years every 100 years, and indicate if species are in the realized pool (in which case persists = TRUE)
    res3 <- res2 %>% 
      group_by(site,order,simul,speciesShortName) %>% 
      summarise(mean_biom_tot=mean(biom_tot)) %>%  # average of biomass on 10 years every 100 years
      mutate(
        persists = case_when(
          speciesShortName %in% filtered ~ TRUE,
          !(speciesShortName %in% filtered) ~ FALSE
        )
      ) # Add a column "persists" indicating if species is above biomass threshold
    
  }
  ALL <- rbind(ALL,res3)
}
# they always persist??


#_______________________________________________________________________________
# Read productivity data ####
sit <- "Bern"
order <- "2sp"
numb <- 1
file2sp <- paste0("analysis_subset_traits/data/raw/Output_ForCEEPS/",sit,
               "/output-cmd2_",site,"_subset_traits_",order,".txt/forceps.",sit,".site_",numb,"_productivityScene.txt")

PROD_2sp <- productivity_specific_adapted(file2sp,sit,order,number = numb) %>% 
  mutate(number = 1)

for (i in c(2:105)){
  file2sp <- paste0("analysis_subset_traits/data/raw/Output_ForCEEPS/",sit,
                     "/output-cmd2_",sit,"_subset_traits_2sp.txt",
                     "/forceps.",sit,".site_",i,"_productivityScene.txt")
  
  prod_2sp <- productivity_specific_adapted(file2sp,sit,order,number = i) %>% 
    mutate(number = i)
  PROD_2sp <- rbind(PROD_2sp,prod_2sp)
}



# Simulations parmi lesquelles il y a au moins une des 4 espèces les plus Di
sp_most_di <- df_di[1:4,] %>% 
  pull(Id)

simul_sp_most_di <- pairs_removed_2 %>% 
  filter(V1 %in% sp_most_di | V2 %in% sp_most_di) %>% 
  pull(number)

toplot <- PROD_2sp %>%
  group_by(simul) %>% 
  summarize(prod = sum(mixture_t_ha)) %>% 
  mutate(di_status = if_else(simul %in% simul_sp_most_di, "Di_absent","Di_present"))

ggplot(toplot,aes(x=di_status, y = prod))+
  geom_point()+
  geom_boxplot() 

# Loreau-Hector ####
species <- read.table("data/processed/correspondence_SName_Id.txt",header=T)

MONOCULTURES <- read.table("data/processed/productivity_monoculture_ALL sites.txt",header=T) %>% 
  mutate(persists_mono = persists) %>% 
  select(-c(persists))


MONO_site <- MONOCULTURES %>% 
  filter(site == sit)

MIXTURES <- PROD_2sp

merged <- merge(MIXTURES,MONO_site, by = "SName" ) %>% 
  group_by(order,simul) %>% 
  rename(site = site.x)

LH1 <- merged %>% 
  rename(YOi=mixture_t_ha,Mi=monoculture) %>% 
  filter(Mi>0) %>% 
  group_by(site,order,simul) %>%
  mutate(persists_mixt=TRUE) %>% # par defaut, mais à modifier!!
  # IL FAUT QUE JE FASSE UN SEUIL SUR LES BIOMASSES POUR LES ESPECES QUI PERSISTENT EN MELANGE
  mutate(nb_sp_realized_pool=sum(persists_mixt)) %>% # nb of species in a community (which is a 2000-year-simul defined by a given regional pool of species)
  mutate(nb_sp_regional_pool=31-simul) %>% # For simul 30, N=31-30=1 sp. For simul 1, N=31-1=30 sp, etc.
  mutate(nb_sp_local_pool = n()) # number of rows. cf. Chauvet 2017.

LH2 <- LH1 %>% 
  mutate(YO=sum(YOi)) %>%
  mutate(RYEi=1/nb_sp_regional_pool) %>% 
  # I divide by the regional pool, which represents the grains "seeded" (See Loreau & Hector, 2001, Nature): 
  # "RYEi = expected relative yield of species i in the mixture, which is simply its proportion seeded or planted"  
  mutate(RYOi=YOi/Mi) %>%
  mutate(YEi=RYEi*Mi) %>%
  mutate(YE=sum(YEi)) %>%
  mutate(DeltaY=YO-YE) %>%
  mutate(DeltaRYi=RYOi-RYEi)

LH3 <- LH2 %>% 
  mutate(Selection = nb_sp_regional_pool * ( mean(DeltaRYi * Mi) - mean(DeltaRYi) * mean(Mi) ),
         Complementarity = nb_sp_regional_pool * mean(DeltaRYi)*mean(Mi),
         DeltaY = nb_sp_regional_pool * mean(DeltaRYi * Mi)) %>%  # idem as YO - YE
  mutate(sum = Selection + Complementarity) %>% 
  # NB: N is the number elements on which we sum (cf. Loreau & Hector, 2001), i.e. "the number of species seeded or planted"
  
  group_by(site,order,simul)

LH4 <- LH3 %>% 
  select(SName,site,order,simul,persists_mono,Mi,persists_mixt,YOi,YE,DeltaY,Selection,Complementarity)

# per species -> aggregated
LH_aggregated <- LH4 %>% 
  group_by(site,order,simul) %>% 
  summarize(YO = mean(YOi),YE = mean(YE),DeltaY=mean(DeltaY),
            Complementarity=mean(Complementarity),Selection=mean(Selection),sum=mean(sum))



# Simulations dans lesquelles il y a au moins une des 4 Sp les plus distinctes ####
toplot_LH <- LH_aggregated %>%
  mutate(di_status = if_else(simul %in% simul_sp_most_di, "Di_absent","Di_present"))

ggplot(toplot_LH,aes(x=di_status, y = DeltaY ))+
  geom_boxplot()+
  geom_point()



