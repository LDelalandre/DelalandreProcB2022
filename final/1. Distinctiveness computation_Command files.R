source("R/Before simulations.R")
source("R/Common variables.R")

# Compute distinctiveness ####
# The distinctiveness calculation that I used
# First, I make a PCA, then I calculate the distinctiveness on the first axis

traits<-read.table("data/Traits of the species_complete.txt",header=T)
c1<-choice_traits_1(traits) # data.frame with the traits of the species
# c1-select(c1,-DDMin) 
# cprime <- select(c1,-c(A1max,A2))
ACP1<-PCA(c1)
# ACPprime <- PCA(cprime)
library(factoextra)
jpeg("figures/PCA on traits.png",width = 800, height = 800, units = "px")
fviz_pca_biplot(ACP1, repel = TRUE, # biplot
                col.var = "#2E9FDF", # Couleur des variables
                col.ind = "#696969"  # Couleur des individus 
)
dev.off()

# position of the species on the first four axis
distACP <- ACP1$ind$coord %>% 
  as.data.frame() %>%
  select(Dim.1    ,   Dim.2     ,   Dim.3      ,  Dim.4) 

distACP$SName<-rownames(distACP) 

distinct_tot <-  distACP %>%
  traits_dist() %>% # here we compute the functional distinctiveness
  tibble::rownames_to_column("SName")

write.table(select(distinct_tot,Di,SName),"data/raw/distinctiveness of the species.txt",row.names=F)

# Write command files ####
distinct_tot <- read.table("data/raw/distinctiveness of the species.txt",header = T)
for (site in SITE){
  Cmd_decr(distinct_tot,length, yearstobejumped, timestep,site)
  Cmd_incr(distinct_tot,length, yearstobejumped, timestep,site)
  Cmd_rand(distinct_tot,length, yearstobejumped, timestep,site)
  Cmd_mono(distinct_tot,length, yearstobejumped, timestep,site)
}



# File: order of removal and names of the species ####
data <- read.table("data/raw/distinctiveness of the species.txt",header=T)

decr <- data[order(data$Di,decreasing=T),]
colnames(decr) <- c("Di_decr","species_decr")
# write.table(decr,"data/decreasing order")

incr <- data[order(data$Di,decreasing=F),]
colnames(incr) <- c("Di_incr","species_incr")
orders <- cbind(incr,decr)

orders$lost_at_simul <- c(2:31)
orders$nb_sp_lost <- orders$lost_at_simul - 1
orders
write.table(orders,"data/removal order and sp names.txt",row.names=F,sep="\t")


# Exploration of other distinctiveness calculations ####
# 0: Gower's distance on all the traits
# 1: Euclidean distance on a subset of traits
# 2: Euclidean distance on a smaller subset of traits
c0 <- choice_traits_0(traits)
c0_di<-traits_dist_gower(c0)
c0_di$Id<-c(0:29)

c1<-choice_traits_1(traits) # data.frame with the traits of the species
c1_di<-traits_dist(c1) # same data.frame plus the distinctiveness
c1_di$Id <- c(0:29)

c2<-choice_traits_2(traits)
c2_di<-traits_dist(c2)
c2_di$Id <- c(0:29)

order0 <- c0_di[order(c0_di$Di),] # species ordered in increasing distinctiveness
order1 <- c1_di[order(c1_di$Di),] # species ordered in increasing distinctiveness
order2 <- c2_di[order(c2_di$Di),] # species ordered in increasing distinctiveness

# sink("data/raw/incr order dist_0.txt"); cat(rownames(c0_di[order(c0_di$Di),])) ; sink()
# sink("data/raw/incr order ditt_1.txt"); cat(rownames(c1_di[order(c1_di$Di),])) ; sink()
# sink("data/raw/incr order dist_2.txt"); cat(rownames(c2_di[order(c2_di$Di),])) ; sink()
