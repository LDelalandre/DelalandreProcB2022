# Script to compute functional distinctiveness and generate ForCEEPS command files #

source("final/0. Packages.R")
source("R/Before simulations.R")
source("R/Common variables.R")


# 1) Compute distinctiveness ####
# First, I make a PCA, then I compute the distinctiveness on the first four axis
traits<-read.table("data/raw/Traits of the species_complete.txt",header=T)
c1<-select_traits(traits) # data.frame with the traits of the species
ACP1<-PCA(c1,graph=F)

# position of the species on the first four axis
distACP <- ACP1$ind$coord %>% 
  as.data.frame() %>%
  select(Dim.1    ,   Dim.2     ,   Dim.3      ,  Dim.4) 

distinct_tot <-  distACP %>%
  traits_dist() %>% # here we compute the functional distinctiveness
  tibble::rownames_to_column("SName")

write.table(select(distinct_tot,Di,SName),"data/code_ForCEEPS_simulations/distinctiveness of the species.txt",row.names=F)

# 2) Write command files ####
distinct_tot <- read.table("data/code_ForCEEPS_simulations/distinctiveness of the species.txt",header = T)
for (site in SITE){
  Cmd_decr(distinct_tot,length, yearstobejumped, timestep,site)
  Cmd_incr(distinct_tot,length, yearstobejumped, timestep,site)
  # Cmd_rand(distinct_tot,length, yearstobejumped, timestep,site)
  Cmd_mono(distinct_tot,length, yearstobejumped, timestep,site)
}



# Order of removal and names of the species
data <- read.table("data/raw/distinctiveness of the species.txt",header=T)

decr <- data[order(data$Di,decreasing=T),]
colnames(decr) <- c("Di_decr","species_decr")

incr <- data[order(data$Di,decreasing=F),]
colnames(incr) <- c("Di_incr","species_incr")
orders <- cbind(incr,decr)

orders$lost_at_simul <- c(2:31)
orders$nb_sp_lost <- orders$lost_at_simul - 1
orders
write.table(orders,"data/raw/removal order and sp names.txt",row.names=F,sep="\t")