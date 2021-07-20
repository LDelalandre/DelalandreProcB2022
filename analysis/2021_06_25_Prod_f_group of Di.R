library(tidyverse)
require(dunn.test)
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
  prod_init <- prod_tot_site_order[1,3]
  for (i in c(1:29)){
    prod_i2 <- prod_tot_site_order[i+1,3]
    prod_i1 <- prod_tot_site_order[i,3]
    # prod_drop <- c(prod_drop, (prod_i1-prod_i2)/prod_i1 ) # drop relative to i-1 prod
    # prod_drop <- c(prod_drop, (prod_i2-prod_i1) ) # absolute drop
    prod_drop <- c(prod_drop, (prod_i1-prod_i2)/prod_init ) # drop relative to 30 sp prod
  }
  prod_drop <- c(prod_drop,prod_tot_site_order[30,3]) # removing last species
  
  species_group %>% 
    mutate(prod_drop = prod_drop) %>% 
    mutate(site=sit,order=order)
  

}

# Relative drop
DROP <- NULL
for (sit in SITE){
  for (order in ORDER){
    species_group <- get_productivity_drop(sit,order)
    DROP <- rbind(DROP,species_group)
  }
}
DROP$prod_drop <- as.numeric(DROP$prod_drop)

# Plot ####

# Boxplot
ggplot(DROP,aes(x=group,y=prod_drop)) +
  facet_wrap(~site)+
  geom_boxplot() +
  ggsave("figures/Drop_f_group/1. All the sites.png",height = 30,width=30, units="cm")

ggplot(DROP %>% filter(!(group=="common")),aes(x=site,y=prod_drop)) +
  geom_boxplot()

# Histogram
ggplot(DROP,aes(x=prod_drop,color=group)) +
  facet_wrap(~site)+
  geom_density() +
  ggsave("figures/Drop_f_group/2. All the sites_density.png",height = 30,width=30, units="cm")


# Tests ####
TESTS <- NULL
for (sit in SITE){
  DROP2 <- DROP %>% 
    filter(site==sit)
  
  ktest <- kruskal.test(x=DROP2$prod_drop,g=DROP2$group)
  ktest$p.value
  
  dtest <- dunn.test::dunn.test(x=DROP2$prod_drop,g=DROP2$group)
  TESTi <- data.frame(site=sit,
             k.pval = rep(ktest$p.value,6),
             comparisons = dtest$comparisons,
             Z = dtest$Z,
             P.adjusted = dtest$P.adjusted)
  TESTS <- rbind(TESTS,TESTi)
}

TESTS %>% 
  filter(Z>0 & P.adjusted<0.05)

table_tests <- TESTS %>%
  mutate(k.pval = round(k.pval,3)) %>% 
  mutate(Z = round(Z,2)) %>% 
  mutate(P.adjusted = round(P.adjusted,3)) %>% 
  mutate(P.adjusted = paste0("(",P.adjusted,")")) %>% 
  unite("Dunn",Z:P.adjusted,sep=" ") %>% 
  spread(comparisons,Dunn)

library("kableExtra")

tabletests<- table_tests %>%
  # transmute(
  #   site=site,
  #   # KW.p.value = k.pval,
  #   p.value = ifelse(k.pval < 0.05,
  #                    cell_spec(k.pval, bold = T),
  #                    cell_spec(k.pval, bold=F)),
  #  #  "cold-common" = "cold - common",
  #  # "cold-drought" = "cold - drought",
  #  #  "cold - small" = "cold - small",
  #  #  "common - drought" ="common - drought",
  #  #  "common - small" ="common - small",
  #  #  "drought - small" = "drought - small"
  # ) %>%
  select(site, everything()) %>%
  kable( escape = F) %>%
  kable_styling("hover", full_width = F)
cat(tabletests, file = "paper_2/#drop_prod.doc")


# just in decreasing order ###

orde <- "decreasing"
DROP3 <- DROP %>% filter(order==orde)

ggplot(DROP3,aes(x=group,y=prod_drop)) +
  facet_wrap(~site)+
  geom_boxplot()
