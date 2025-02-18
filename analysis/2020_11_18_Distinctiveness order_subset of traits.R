# Following Annette's remarks

source("R/Common variables.R")
source("R/Analysis_data.R")
source("R/Monocultures_functions.R")
source("R/Before simulations.R")

comp_fct_dist <- function(traits){
  c1<-choice_traits_1(traits) # data.frame with the traits of the species
  c1 <- select(c1,-c(A1max,A2))# remove traits that say explicitely that a species is a gymnosperm
  ACP1<-PCA(c1,graph=F)
  distACP <- ACP1$ind$coord %>% 
    as.data.frame() %>%
    select(Dim.1    ,   Dim.2     ,   Dim.3      ,  Dim.4) %>%
    traits_dist() # here we compute the functional distinctiveness
  distACP$SName<-rownames(distACP)
  
  distinct_tot <-  distACP
  # distinct_tot[order(distinct_tot$Di),]
  distinct_tot
}

traits <- read.table("data/traits of the species_complete.txt",header=T)


# no envt resp traits / Di ####
# Compute functional distinctiveness on traits not related to envt response (Annette's remark)
traits.simulations <- select(traits,Name,SName,S, HMax, AMax,   G, DDMin, WiTN, WiTX, DrTol, NTol, Brow,   Ly, La,A1max, A2)
dist.all.traits <- comp_fct_dist(traits.simulations)

traits.not.envt <- select(traits,Name,SName,S, HMax, AMax,   G, Brow,   Ly, La,A1max, A2)
dist.not.envt <- comp_fct_dist(traits.not.envt)

dist.all.traits[order(dist.all.traits$Di),]
dist.not.envt[order(dist.not.envt$Di),]


all <- dist.all.traits$Di
no.envt <- dist.not.envt$Di
cor(all,no.envt,method="spearman")
plot(all,no.envt)

# envt resp traits / Di
# compute functional distinctiveness  related to envt response (Annette's remark)
traits.envt <- select(traits.simulations,Name,SName,DDMin, WiTN, WiTX, DrTol, NTol)
dist.envt <- comp_fct_dist(traits.envt)

envt <- dist.envt$Di
cor(a,envt,method="spearman")
plot(all,envt)
