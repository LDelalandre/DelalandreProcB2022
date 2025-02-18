# Function to compute functional distinctiveness from a set of traits
source("R/Before simulations.R")

# species name must be provided as rownames
comp_fct_dist <- function(traits){
  # if ( grep("Name", colnames(traits) )!=0 ){
  #   traits <- traits %>% 
  #     mutate(Name=NULL) %>% 
  #     mutate(SName=NULL)
  # }
  ACP1<-FactoMineR::PCA(traits,graph=F)
  if(dim(ACP1$ind$coord)[2]>=4){
    distACP <- ACP1$ind$coord %>% 
      as.data.frame() %>%
      select(Dim.1    ,   Dim.2     ,   Dim.3      ,  Dim.4) %>%
      traits_dist() # here we compute the functional distinctiveness
  } else if (dim(ACP1$ind$coord)[2]==3) {
    distACP <- ACP1$ind$coord %>% 
      as.data.frame() %>%
      select(Dim.1    ,   Dim.2     ,   Dim.3) %>%
      traits_dist() # here we compute the functional distinctiveness
  } else if (dim(ACP1$ind$coord)[2]==2){
    distACP <- ACP1$ind$coord %>% 
      as.data.frame() %>%
      select(Dim.1    ,   Dim.2) %>%
      traits_dist() # here we compute the functional distinctiveness
  } else {
    distACP <- ACP1$ind$coord %>% 
      as.data.frame() %>%
      select(Dim.1) %>%
      traits_dist() # here we compute the functional distinctiveness
  }
  
  distACP$SName<-rownames(distACP)
  
  distinct_tot <-  distACP
  distinct_tot
}
