# charge all the "4. main figyre removal expermients.R" except the "CODE" section.
# Then :

# Environmental conditions
sites <- read.table("data/Site description.txt",header=T)

# Sites and coordinates used hereafter
COORD <- sites %>% 
  select(Site,Temp_moy,Annual_ppt) %>% 
  arrange(Temp_moy,Annual_ppt)

measure <- "productivity_tot"
datatoplot <- read.table(paste0("data/processed/",measure,"_with interval_median.txt"),header=T)

# generate graphs for each site
PLOT <- list()
i <- 0
for (sit in COORD$Site){
  i <- i+1
  result <- 
    datatoplot %>% 
    filter(site==sit)
  PLOT[[i]] <- ggplot(result,aes(x=simul-1,y=decreasing)) +  
    gg_removal_productivity()+
    ggtitle(sit) +
    geom_line(aes(x=simul-1,y=random_1))+
    geom_line(aes(x=simul-1,y=random_2))+
    geom_line(aes(x=simul-1,y=random_3))+
    geom_line(aes(x=simul-1,y=random_4))+
    geom_line(aes(x=simul-1,y=random_5))+
    geom_line(aes(x=simul-1,y=random_6))+
    geom_line(aes(x=simul-1,y=random_7))+
    geom_line(aes(x=simul-1,y=random_8))+
    geom_line(aes(x=simul-1,y=random_9))+
    geom_line(aes(x=simul-1,y=random_10))+
    geom_line(aes(x=simul-1,y=random_11))+
    geom_line(aes(x=simul-1,y=random_12))+
    geom_line(aes(x=simul-1,y=random_13))+
    geom_line(aes(x=simul-1,y=random_14))+
    geom_line(aes(x=simul-1,y=random_15))+
    geom_line(aes(x=simul-1,y=random_16))+
    geom_line(aes(x=simul-1,y=random_17))+
    geom_line(aes(x=simul-1,y=random_18))+
    geom_line(aes(x=simul-1,y=random_19))+
    geom_line(aes(x=simul-1,y=random_20))+
    geom_line(aes(x=simul-1,y=random_21))+
    geom_line(aes(x=simul-1,y=random_22))+
    geom_line(aes(x=simul-1,y=random_23))+
    geom_line(aes(x=simul-1,y=random_24))+
    geom_line(aes(x=simul-1,y=random_25))+
    geom_line(aes(x=simul-1,y=random_26))+
    geom_line(aes(x=simul-1,y=random_27))+
    geom_line(aes(x=simul-1,y=random_28))+
    geom_line(aes(x=simul-1,y=random_29))+
    geom_line(aes(x=simul-1,y=random_30))
  
}


