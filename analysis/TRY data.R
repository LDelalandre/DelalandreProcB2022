# Packages ####
library(dplyr) ;library(magrittr);library(purrr);library(ggplot2)

# Functions ####
source("R/comp_fct_dist.R")

remove.na <- function(table){ # removes rows containing at least one NA from a data frame
  contains.na <- apply(table,1,function(x)any(is.na(x)))
  sp.removed <- as.numeric(which(contains.na==T))
  table[-sp.removed,]
}

shorten_sp_name <- function(x){ # change species name from "Genera species" to "GSpe" in the vector of names
  space <- unlist(gregexpr(" ", x)) # where the space between genera and species name is
  genera <- substr(x, start = 1, stop = 1) # first letter of genera
  sp1 <- toupper(substr(x, start = space+1, stop = space+1)) # upper case first letter of species
  sp23 <- substr(x, start = space+2, stop = space+3) # lower case two following letters
  paste0(genera,sp1,sp23)
}

change_sp_name <- function(df2){ 
  # applies shorten_sp_name to the data frame (and not just to a vector)
  # NB species name should be in rownames
  df3 <- df2 %>% 
    mutate(SName = purrr::map_chr(rownames(df2),shorten_sp_name) )
  df3$SName[which(df3$SName=="PMug")] <- "PMon" # P. mugo = P. montana (2 names for the same species)
  rownames(df3) <- df3$SName
  select(df3,-SName)
}

bind_forceeps <- function (dist.ForCEEPS,dist.focus) { # groupe distinctiveness output of forceeps and another dataset
  dist.forceeps <- filter(dist.ForCEEPS,SName %in% dist.focus$SName)
  toplot <- data.frame(dist.forceeps$SName,dist.forceeps$Di,dist.focus$Di)
  colnames(toplot) <- c("SName","Di.forceeps","Di.focus")
  toplot
}

# Files ####
extract <- read.csv("data/TRY/leo_dataset.csv") # full dataset

#_____________________________________________________________
# FORCEEPS #### 
# Compute Distinctiveness on ForCEEPS parameters
traits <- read.table("data/traits of the species_complete.txt",header=T)
row.names(traits) <- traits$SName
traits.simulations <- select(traits,Name,SName, HMax, AMax,   G, DDMin, WiTN, WiTX, DrTol, NTol, Brow,   Ly, La,S,A1max, A2)
traits.simulations <- select(traits.simulations,-c(Name,SName))
dist.ForCEEPS <- comp_fct_dist(traits.simulations)


#______________________
# Tidy TRY dataset ####

# Filter the data
extract2 <- extract %>% 
  filter(TraitName!="") %>% # remove lines where no trait is measured (and even defined)
  filter(!is.na(StdValue))   # remove lines where there is no value measured for the trait
  # mutate(TraitName = droplevels(TraitName))

# Organize it such as species are in rows and mean values for their traits in columns
pivot.table <- with(extract2,tapply(StdValue,list(AccSpeciesName,TraitName),mean)) # I compute the mean of trait values for each species
pivot.table <- data.frame(pivot.table)

# Observations
dim(extract2)[1]/dim(extract)[1] # I removed 97% of the original dataset!! But I still have 44183 observations
class(extract2$StdValue) # I only have numeric traits
Filter(function(x)!any(is.na(x)), pivot.table) # remove traits for which at least one species has no value. No trait remains!

#_______________________________________
# WESTOBY: SLA + Height + Seed Mass ####

# Test the results with the 3 possible measures of SLA
PVAL <- c()
STAT <- c()
SLAmeasures <- colnames(pivot.table)[grep("Leaf.area.per.leaf.dry.mass",colnames(pivot.table))]

for (mesSLA in SLAmeasures){
  # Select the traits used in Westoby
  west <- pivot.table %>% 
    select(c("Seed.dry.mass", starts_with("Leaf.area.per.leaf.dry.mass"), "Plant.height.vegetative") ) %>% 
    select(c("Seed.dry.mass",
             all_of(mesSLA),
             "Plant.height.vegetative")) %>%
    rename(SLA = starts_with("Leaf.area.per.leaf.dry.mass..specific.leaf.area..SLA.or.1.LMA...undefined.if.petiole.is.in..or.excluded"))
  
  # remove species with NA
  west2 <- remove.na(west) # RQ: it removes subspecies such as Sorbus aucuparia subsp...: lack of measurements
  # Change species names to fit with ForCEEPS SName
  west3 <- change_sp_name(west2)
  # Compute functional distinctiveness
  dist.west0 <- comp_fct_dist(west3)
  # Order it in the same way as Di computed on ForCEEPS
  dist.west <- dist.west0[order(match(dist.west0$SName, dist.ForCEEPS$SName)),] 
  # Compare it with ForCEEPS order
  toplot <- bind_forceeps(dist.ForCEEPS,dist.west)
  
  # Correlation test
  cortest <- cor.test(toplot$Di.forceeps,toplot$Di.focus,method="spearman")
  PVAL <- c(PVAL,cortest$p.value)
  STAT <- c(STAT,cortest$estimate)
}
TESTScor <- data.frame(SLAmeasures,PVAL,STAT)

TESTScor <- TESTScor %>% 
  mutate(SLAmeasures =  purrr::map_chr(TESTScor$SLAmeasures,
                               function(x) stringr:::str_remove(x,"Leaf.area.per.leaf.dry.mass..specific.leaf.area..SLA.or.1.LMA...")))

ggplot(toplot,aes(x=Di.forceeps,y=Di.focus,label=SName)) +
  labs(x="Distinctiveness computed with ForCEEPS parameters",y="Distinctiveness computed with TRY, using Westoby's SLA") +
  geom_point() +
  geom_label()

#____________________________________________
# DIAZ: SM + LMA + H + SSD + LA  + Nmass ####

# Test the results with the 3 (resp. 7) possible measures of SLA (resp. LA) 
PVAL <- c() ; STAT <- c() ; SLA <- c() ; LA <- c()

SLAmeasures <- colnames(pivot.table)[grep("Leaf.area.per.leaf.dry.mass",colnames(pivot.table))]
LAmeasures <- colnames(pivot.table)[grep("Leaf.area..in.case.of.compound.leaves.",colnames(pivot.table))]

for (mesSLA in SLAmeasures){
  for(mesLA in LAmeasures){
    diaz <- pivot.table %>% 
      select(c("Seed.dry.mass",
               all_of(mesSLA),
               "Plant.height.vegetative",
               "Stem.specific.density..SSD..or.wood.density..stem.dry.mass.per.stem.fresh.volume.",
               all_of(mesLA),
               "Leaf.nitrogen..N..content.per.leaf.dry.mass")) 
    colnames(diaz) <- c("SM","SLA", "H", "SSD", "LA", "Nmass")
    diaz$LMA <- purrr::map_dbl(diaz$SLA,function(x)1/x)
    diaznoSLA <- select(diaz,-SLA)
    
    # remove species with NA
    diaz2 <- remove.na(diaznoSLA)
    # Change species names to fit with ForCEEPS SName
    diaz3 <- change_sp_name(diaz2)
    # Compute functional distinctiveness
    dist.diaz0 <- comp_fct_dist(diaz3)
    # Order it in the same way as Di computed on ForCEEPS
    dist.diaz <- dist.diaz0[order(match(dist.diaz0$SName, dist.ForCEEPS$SName)),] 
    # Compare it with ForCEEPS order
    toplot <- bind_forceeps(dist.ForCEEPS,dist.diaz)
    
    # Correlation test
    cortest <- cor.test(toplot$Di.forceeps,toplot$Di.focus,method="spearman")
    PVAL <- c(PVAL,cortest$p.value)
    STAT <- c(STAT,cortest$estimate)
    SLA <- c(SLA,mesSLA)
    LA <- c(LA,mesLA)
  }
}

TESTScor <- data.frame(SLA,LA,PVAL,STAT)
TESTScor <- TESTScor %>% 
  mutate(SLA =  purrr::map_chr(TESTScor$SLA,
                               function(x) stringr:::str_remove(x,"Leaf.area.per.leaf.dry.mass..specific.leaf.area..SLA.or.1.LMA..."))) %>% 
  mutate(LA = purrr::map_chr(TESTScor$LA,
                             function(x) stringr:::str_remove(x,"Leaf.area..in.case.of.compound.leaves."))
  )


ggplot(TESTScor,aes(PVAL))+
  geom_histogram(bins=40)+
  geom_vline(xintercept=0.05,colour="red")

ggplot(toplot,aes(x=Di.forceeps,y=Di.focus,label=SName)) +
  labs(x="Distinctiveness computed with ForCEEPS parameters",y="Distinctiveness computed with TRY, using Westoby's SLA") +
  geom_point()+
  geom_label()
  geom_smooth(method=lm)


# PRODUCTIVITE: SLA + LNC + H ####
# "Leaf.nitrogen..N..content.per.leaf.dry.mass"
# "Plant.height.vegetative" 
  
  PVAL <- c() ; STAT <- c() ; SLA <- c() ; LA <- c()
  
  SLAmeasures <- colnames(pivot.table)[grep("Leaf.area.per.leaf.dry.mass",colnames(pivot.table))]
  
  for (mesSLA in SLAmeasures){
      pro <- pivot.table %>% 
        select(c(all_of(mesSLA),
                 "Plant.height.vegetative",
                 "Leaf.nitrogen..N..content.per.leaf.dry.mass")) 
      colnames(pro) <- c("SLA", "H", "Nmass")
      
      # remove species with NA
      pro2 <- remove.na(pro)
      # Change species names to fit with ForCEEPS SName
      pro3 <- change_sp_name(pro2)
      # Compute functional distinctiveness
      dist.pro0 <- comp_fct_dist(pro3)
      # Order it in the same way as Di computed on ForCEEPS
      dist.pro<- dist.pro0[order(match(dist.pro0$SName, dist.ForCEEPS$SName)),] 
      # Compare it with ForCEEPS order
      toplot <- bind_forceeps(dist.ForCEEPS,dist.pro)
      
      # Correlation test
      cortest <- cor.test(toplot$Di.forceeps,toplot$Di.focus,method="spearman")
      PVAL <- c(PVAL,cortest$p.value)
      STAT <- c(STAT,cortest$estimate)
      SLA <- c(SLA,mesSLA)
      LA <- c(LA,mesLA)
    }
  
  TESTScor <- data.frame(SLA,PVAL,STAT)
  TESTScor <- TESTScor %>% 
    mutate(SLA =  purrr::map_chr(TESTScor$SLA,
                                 function(x) stringr:::str_remove(x,"Leaf.area.per.leaf.dry.mass..specific.leaf.area..SLA.or.1.LMA..."))) 


# RESPONSE TRAITS to be defined ####

pp <- pivot.table %>% 
    select(!starts_with("Leaf.area.per.leaf.dry.mass")) %>% 
    select(!starts_with("Leaf.area..in.case.of.compound.leaves"))
    
  
