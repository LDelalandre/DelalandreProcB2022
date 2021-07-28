source("R/Common variables.R")
library("dplyr")

# Functions ####
bind_mono <- function(MONO_ALL,MONO1,ord,sim){
  # rbind monoculture data frame to itself, but adding a column giving the order of removal and number of the simul
  # It gives me a global data frame to which I can add the productivities in mixture of all the species.
  Id_sp_to_keep <-EXTR_ORDERS[which(rownames(EXTR_ORDERS)==ord),][sim:30] # Id of species present in the simul
  MONO1 <- MONO1 %>% 
    mutate(order=ord,simul=sim) %>% 
    filter(Id %in% Id_sp_to_keep)
  bind_rows(MONO_ALL,MONO1)
}

#___________________________________________________
# Extract orders of removal from command files ####
EXTR_ORDERS <- NULL
sit <- "Sion" # I can take any site

cmd <- read.csv(paste0("data/code_ForCEEPS_simulations/Cmd files/cmd2_",sit,"_increasing.txt"))
extracted <- extract_order(cmd)
extr_split <- stringr::str_split(extracted, " ")
EXTR_ORDERS <- rbind(EXTR_ORDERS, extr_split[[1]] )

cmd <- read.csv(paste0("data/code_ForCEEPS_simulations/Cmd files/cmd2_",sit,"_decreasing.txt"))
extracted <- extract_order(cmd)
extr_split <- stringr::str_split(extracted, " ")
EXTR_ORDERS <- rbind(EXTR_ORDERS, extr_split[[1]] )

for(j in c(1:30)){
  cmd <- read.csv(paste0("data/code_ForCEEPS_simulations/Cmd files/cmd2_",sit,"_random_",j,".txt"))
  extracted <- extract_order(cmd)
  extr_split <- stringr::str_split(extracted, " ")
  EXTR_ORDERS <- rbind(EXTR_ORDERS, extr_split[[1]] )
}
rownames(EXTR_ORDERS) <- ORDER
write.table(EXTR_ORDERS,"data/Orders of removal of cmd files.txt",col.names=F)
EXTR_ORDERS <- read.table("data/Orders of removal of cmd files.txt")

#____________________________________________________________________
# Generate the Monoculture data frame ####
MONO <-read.table("data/processed/productivity_monoculture_ALL sites.txt",header=T)
MONO_ALL <- NULL
for(ord in ORDER){
  for(sim in c(1:30)){
    MONO_ALL <- bind_mono(MONO_ALL,MONO,ord,sim)
    }
}
MONO_ALL2 <- MONO_ALL %>% 
  select(Di,SName,monoculture,site,order,simul,persists,Id) %>% 
  rename(monoculture_t_ha=monoculture,persists_mono=persists,species=SName)

# Build a data frame with all the mixtures togethers ####
MIXT_ALL <- NULL
for (sit in SITE){
  MIXT <- read.table(paste0("data/processed/productivity_specific_",sit,".txt"),header=T)
  MIXT_ALL <- bind_rows(MIXT_ALL,MIXT)
}

# Pool the data frames together ####
MIXT_ALL2 <- MIXT_ALL %>%
  tibble(.) %>% 
  arrange(site,order,simul,species)
 
MONO_ALL3 <- MONO_ALL2 %>%
  tibble(.) %>% 
  arrange(site,order,simul,species)


# combine them
ALL_POOLED <- MONO_ALL3 %>% 
  mutate(persists_mixt=MIXT_ALL2$persists) %>% 
  mutate(mixture_t_ha=MIXT_ALL2$mixture_t_ha)

write.table(ALL_POOLED,"data/processed/LH_productivity_specific_every condition.txt",row.names=F,sep="\t")

