# Packages ####
library(dplyr) ;library(magrittr);library(purrr);library(ggplot2)

# Functions ####
source("R/Common variables.R")
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

bind_forceeps <- function (dist.ForCEEPS,dist.focus) { # group distinctiveness output of forceeps and another dataset
  dist.forceeps <- filter(dist.ForCEEPS,SName %in% dist.focus$SName)
  toplot <- data.frame(dist.forceeps$SName,dist.forceeps$Di,dist.focus$Di)
  colnames(toplot) <- c("SName","Di.forceeps","Di.focus")
  toplot
}

compare_Di <- function(dist.ForCEEPS,traitdata){
  # dist.ForCEEPS is the output of comp_fct_dist, with species as row names,
  # and colnames =  Dim.1       Dim.2        Dim.3        Dim.4       Di SName.
  
  # traitsdata is the data frame of the traits I want to use to compute distinctiveness,
  # species as rownames, and traits as columns, ex: SM       SLA         H.
  
  # output: data framer with the distinctiveness value per species computed with the two datasets.
  
  # remove species with NA
  traitdata2 <- remove.na(traitdata)
  # Change species names to fit with ForCEEPS SName
  traitdata3 <- change_sp_name(traitdata2)
  # Compute functional distinctiveness
  dist.traitdata0 <- comp_fct_dist(traitdata3)
  # Order it in the same way as Di computed on ForCEEPS
  dist.traitdata<- dist.traitdata0[order(match(dist.traitdata0$SName, dist.ForCEEPS$SName)),] 
  # Compare it with ForCEEPS order
  toplot <- bind_forceeps(dist.ForCEEPS,dist.traitdata)
  toplot
}

plot_Di_Di <- function(toplot){
  # toplot is the output of the compare_Di function
  # the function plots the Di of the species computed from the two datasets used to generate toplot.
  ggplot(toplot,aes(x=Di.forceeps,y=Di.focus,label=SName)) +
    labs(x="Distinctiveness computed on ForCEEPS parameters",y="Distinctiveness computed on TRY traits") +
    geom_point()+
    geom_label()
    # geom_smooth(method=lm)
}



#______________________
# I) Tidy TRY dataset - DO IT ONCE ONLY ####
# Files ####
extract <- read.csv("data/TRY/leo_dataset.csv") # full dataset

# Filter the data
extract2 <- extract %>% 
  filter(TraitName!="") %>% # remove lines where no trait is measured (and even defined)
  filter(!is.na(StdValue))   # remove lines where there is no value measured for the trait
  # mutate(TraitName = droplevels(TraitName))

# Organize it so that species are in rows and mean values for their traits in columns
pivot.table <- with(extract2,tapply(StdValue,list(AccSpeciesName,TraitName),mean)) # I compute the mean of trait values for each species
pivot.table <- data.frame(pivot.table)

write.table(pivot.table,"data/processed/traits from TRY.txt")

#_____________________________________________________________
# II) ANALYSIS ####
# readRDS("data/TRY/try_8066_50997.Rds") # sorbus aucuparia

# FORCEEPS
# Compute Distinctiveness on ForCEEPS parameters
traits <- read.table("data/traits of the species_complete.txt",header=T)
row.names(traits) <- traits$SName
traits.simulations <- select(traits,Name,SName, HMax, AMax,   G, DDMin, WiTN, WiTX, DrTol, NTol, Brow,   Ly, La,S,A1max, A2)
traits.simulations <- select(traits.simulations,-c(Name,SName))
dist.ForCEEPS <- comp_fct_dist(traits.simulations)

# TRY
pivot.table <- read.table("data/processed/traits from TRY.txt")


#____________________________________________
# DIAZ: SM + LMA + H + SSD + LA  + Nmass ####
# I use SLA and LA measures including petiole
diaz <- pivot.table %>% 
  select(c("Seed.dry.mass",
           "Leaf.area.per.leaf.dry.mass..specific.leaf.area..SLA.or.1.LMA...petiole.included",
           "Plant.height.vegetative",
           "Stem.specific.density..SSD..or.wood.density..stem.dry.mass.per.stem.fresh.volume.",
           "Leaf.area..in.case.of.compound.leaves..leaf..petiole.included.",
           "Leaf.nitrogen..N..content.per.leaf.dry.mass")) 
colnames(diaz) <- c("SM","SLA", "H", "SSD", "LA", "Nmass")
diaz$LMA <- purrr::map_dbl(diaz$SLA,function(x)1/x)
diaznoSLA <- select(diaz,-SLA)

# Compute distinctiveness on these traits and compare it to ForCEEPS
toplot <- compare_Di(dist.ForCEEPS,diaznoSLA)
cortest <- cor.test(toplot$Di.forceeps,toplot$Di.focus,method="spearman")
diazplot <- 
  plot_Di_Di(toplot)+
  ggtitle("6 traits of the Global Spectrum (SM + LMA + H + SSD + LA  + Nmass)")+
  ggpubr::stat_cor(method="spearman",label.x.npc = 0.3)

ggsave(diazplot,filename = paste0("paper_2/TRY_Diaz.png"),
   width = 23,
   height = 18,
   units = "cm",
   dpi = 300)

# see dependency of rho value to one species (CSat)
toplot %>% 
  subset(SName!="CSat") %>% 
  {cor.test(.$Di.forceeps,.$Di.focus,method="spearman")}

#_________________________________
# PRODUCTIVITE: SLA + LNC + H ####
# "Leaf.nitrogen..N..content.per.leaf.dry.mass"
# "Plant.height.vegetative" 
  
pro <- pivot.table %>% 
  select(c("Leaf.area.per.leaf.dry.mass..specific.leaf.area..SLA.or.1.LMA...petiole.included",
            "Plant.height.vegetative",
            "Leaf.nitrogen..N..content.per.leaf.dry.mass")) 
colnames(pro) <- c("SLA", "H", "Nmass")

# Compute distinctiveness on these traits and compare it to ForCEEPS
toplot <- compare_Di(dist.ForCEEPS,pro)
cortest <- cor.test(toplot$Di.forceeps,toplot$Di.focus,method="spearman")  
dim(toplot)[1] # number of species kept
prodplot <- 
  plot_Di_Di(toplot)+
  ggtitle("Effect traits (SLA, H, Nmass)")+
  ggpubr::stat_cor(method="spearman",label.y.npc = 0.99)

ggsave(prodplot,filename = paste0("paper_2/TRY_effect traits.png"),
   width = 23,
   height = 18,
   units = "cm",
   dpi = 300)

#_________________________________
# Combine the two plots ####

combined <- cowplot::plot_grid(prodplot,diazplot, labels=c("A", "B"), ncol = 2, nrow = 1)
ggsave(combined,filename="paper_2/TRY combined.png",
       width = 35,
       height = 18,
       units = "cm",
       dpi = 300)

cowplot::save_plot("paper_2/TRY combined.png", combined,
          base_aspect_ratio = 1.3 #Laisser de la place pour la légende
)
