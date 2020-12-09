# Packages ####
library(dplyr) ;library(magrittr);library(purrr)

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

# FORCEEPS #### Compute Distinctiveness on ForCEEPS parameters
traits <- read.table("data/traits of the species_complete.txt",header=T)
row.names(traits) <- traits$SName
traits.simulations <- select(traits,Name,SName, HMax, AMax,   G, DDMin, WiTN, WiTX, DrTol, NTol, Brow,   Ly, La,S,A1max, A2)
traits.simulations <- select(traits.simulations,-c(Name,SName))
dist.ForCEEPS <- comp_fct_dist(traits.simulations)

bind_forceeps <- function (dist.ForCEEPS,dist.focus) { # groupe distinctiveness output of forceeps and another dataset
  dist.forceeps <- filter(dist.ForCEEPS,SName %in% dist.focus$SName)
  toplot <- data.frame(dist.forceeps$SName,dist.forceeps$Di,dist.focus$Di)
  colnames(toplot) <- c("SName","Di.forceeps","Di.focus")
  toplot
}

# Files ####
extract <- read.csv("data/TRY/leo_dataset.csv") # full dataset
traits.available <- read.table("data/TRY/Traits available.txt",sep="\t",header=T)


# Explore the traits available from TRY ####
# there are 889 traits available for these species
# but in come cases, little data is available
traits.available <- select(traits.available,-X)

test <- traits.available %>%
  select(-TraitID,-Trait) %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum>100)
hist(test$sum,breaks=50)
which(test$sum==295664)
traits.available[109,]

traits.kept.1 <- subset(traits.available,.Alnus.viridis.ssp..sinuata.>0)

traits.available2 <- select(traits.available,-starts_with(".Alnus.viridis"),-starts_with(".Sorbus.aucuparia"))
       
i <- 1
rows.to.keep <- c()
for (i in 1:dim(traits.available2)[1]){
  if ( length(which(traits.available2[i,-1]==0)) ==0 ) {
    rows.to.keep <- c(rows.to.keep, i )
  }
}
rows.to.keep

traits.available[c(71,  98, 131),]


# Tidy TRY dataset ####

# Filter the data
extract2 <- extract %>% 
  filter(TraitName!="") %>% # remove lines where no trait is measured (and even defined)
  filter(!is.na(StdValue)) %>%  # remove lines where there is no value measured for the trait
  mutate(TraitName = droplevels(TraitName))

# Organize it such as species are in rows and mean values for their traits in columns
pivot.table <- with(extract2,tapply(StdValue,list(AccSpeciesName,TraitName),mean)) # I compute the mean of trait values for each species
pivot.table <- data.frame(pivot.table)



# Observations
dim(extract2)[1]/dim(extract)[1] # I removed 97% of the original dataset!! But I still have 44183 observations
class(extract2$StdValue) # I only have numeric traits

Filter(function(x)!any(is.na(x)), pivot.table) # remove traits for which at least one species has no value. No trait remains!
pivot2 <- pivot.table %>% # reduce data frame size...
  select(!starts_with("Leaf.area")& !starts_with("Species.environmental"))



# WESTOBY: SLA + Height + Seed Mass ####
west <- pivot.table %>% 
  select(c("Seed.dry.mass", starts_with("Leaf.area.per.leaf.dry.mass"), "Plant.height.vegetative") ) %>% 
  select(c("Seed.dry.mass",
           "Leaf.area.per.leaf.dry.mass..specific.leaf.area..SLA.or.1.LMA...undefined.if.petiole.is.in..or.excluded",
           "Plant.height.vegetative")) %>% 
  mutate(Species=rownames(pivot.table)) %>% 
  mutate(SLA = Leaf.area.per.leaf.dry.mass..specific.leaf.area..SLA.or.1.LMA...undefined.if.petiole.is.in..or.excluded) %>% 
  select(-Leaf.area.per.leaf.dry.mass..specific.leaf.area..SLA.or.1.LMA...undefined.if.petiole.is.in..or.excluded) %>% 
  tibble::column_to_rownames('Species')

# remove species with NA
west2 <- remove.na(west)
# Change species names to fit with ForCEEPS SName
west3 <- change_sp_name(west2)
# Compute functional distinctiveness
dist.west0 <- comp_fct_dist(west3)
# Order it in the same way as Di computed on ForCEEPS
dist.west <- dist.west0[order(match(dist.west0$SName, dist.ForCEEPS$SName)),] 
# Compare it with ForCEEPS order
toplot <- bind_forceeps(dist.ForCEEPS,dist.west)

cor(toplot$Di.forceeps,toplot$Di.focus,method="spearman")

ggplot(toplot,aes(x=Di.forceeps,y=Di.focus,label=SName)) +
  labs(x="Distinctiveness computed with ForCEEPS parameters",y="Distinctiveness computed with TRY, using Westoby's SLA") +
  geom_point() +
  geom_label()

# DIAZ: SM + LMA + H + SSD + LA  + Nmass ####
diaz <- pivot.table %>% 
  select(c("Seed.dry.mass",
           "Leaf.area.per.leaf.dry.mass..specific.leaf.area..SLA.or.1.LMA...undefined.if.petiole.is.in..or.excluded",
           "Plant.height.vegetative",
           "Stem.specific.density..SSD..or.wood.density..stem.dry.mass.per.stem.fresh.volume.",
           "Leaf.area..in.case.of.compound.leaves.undefined.if.leaf.or.leaflet..undefined.if.petiole.is.in..or.excluded.",
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

cor(toplot$Di.forceeps,toplot$Di.focus,method="spearman")

ggplot(toplot,aes(x=Di.forceeps,y=Di.focus,label=SName)) +
  labs(x="Distinctiveness computed with ForCEEPS parameters",y="Distinctiveness computed with TRY, using Westoby's SLA") +
  geom_point()+
  geom_label()
  geom_smooth(method=lm)


# PRODUCTIVITE: SLA + LNC + H ####
# "Leaf.nitrogen..N..content.per.leaf.dry.mass"
# "Plant.height.vegetative" 

# RESPONSE TRAITS to be defined ####


