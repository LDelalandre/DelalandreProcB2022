library(tidyverse)
source("R/Common variables.R")

# Di_order <- read.table("data/raw/distinctiveness of the species.txt",header=T) %>% 
#   arrange(desc(Di))

SName_Id <- read.table("data/correspondence_SName_Id.txt",header=T)
order_removal <- read.table("data/Orders of removal of cmd files.txt",header=F) %>% 
  column_to_rownames("V1")

# Attribute distinct species to 3 groups
small <- c("BPen","AVir")
cold <- c("PCem","PMon","LDec","PSyl","PAbi")
drought <- c("TBac","AAlb","QRob")
common <- Di_order[11:30,]$SName

get_group <- function(SName){
  if (SName %in% small){
    "small"
  } else if (SName %in% cold){ 
    "cold"
  } else if (SName %in% drought){
      "drought"
  } else {
      "common"
  }
}

Id_group <- SName_Id %>% 
  mutate(group = map_chr(SName,get_group))



# Look a drop in ecosystem productivity
prod_tot <- read.table("data/processed/productivity_tot_with interval_median.txt",header=T)

get_productivity_drop <- function(site,order){
  removal_Id <- order_removal[order,] %>% 
    as.vector()
  
  species_group <- Id_group %>% 
    arrange(match(Id, removal_Id))
  
  prod_tot_site_order <- prod_tot %>% 
    filter(site==sit) %>% 
    select(site,simul,order)
  
  prod_drop <- c()
  for (i in c(1:29)){
    prod_drop <- c(prod_drop,prod_tot_site_order[i+1,3]-prod_tot_site_order[i,3])
  }
  prod_drop <- c(prod_drop,prod_tot_site_order[30,3]) # removing last species
  
  species_group$prod_drop <- prod_drop
  
  # plot ####
  ggplot(species_group,aes(x=group,y=prod_drop)) +
    geom_boxplot() +
    ggtitle(paste(sit,order))
}

sit <- "GrandeDixence"
order <- "random_10"
for (sit in SITE){
  for (order in ORDER){
    fig <- get_productivity_drop(sit,order)
    ggsave(paste0("figures/Drop_f_group/",sit,"_",order,".png"),plot=fig) 
  }
}

