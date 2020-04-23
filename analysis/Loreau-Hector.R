#en faire une fonction
data <- read.table("data/processed/specific_biom_prod_complete.txt",header=T)
if (measure=="biomass_tot"){
  data <- select(data,species, mixture.t.ha.,monoculture.t.ha.,site, order, simul)
  names(data)[names(data) == "mixture.t.ha."] <- "YOi"
  names(data)[names(data) == "monoculture.t.ha."] <- "Mi"
}else{
  data <- select(data,species, prod_mixture,prod_monoculture,site, order, simul)
  names(data)[names(data) == "prod_mixture"] <- "YOi"
  names(data)[names(data) == "prod_monoculture"] <- "Mi"
}
data <- data %>% group_by(site,order,simul)

data2 <- data %>% 
  mutate(YO=mean(YOi)) %>%
  mutate(RYEi=1/(31-simul)) %>%# 1/N, N being the initial nb of species. For simul 30, N=31-30=1 sp. For simul 1, N=31-1=30 sp, etc.
  mutate(RYOi=YOi/Mi) %>%
  mutate(YEi=RYEi*Mi) %>%
  mutate(YE=sum(YEi)) %>%
  mutate(DeltaY=YO-YE) %>%
  mutate(DeltaRYi=RYOi-RYEi) %>%
  mutate(Mavg = mean(Mi)) %>%
  mutate(DeltaRYavg = mean(DeltaRYi)) %>%
  
  mutate(Cpltarity = (31-simul)*DeltaRYavg*Mavg) %>%
  mutate(Selection = DeltaY - Cpltarity)
  

data3 <- summarise(data2,DeltaY=mean(DeltaY),Cpltarity=mean(Cpltarity),Selection=mean(Selection)) 
# I use mean, because all the values are already the same for one group for DeltaRY etc.

sit <- SITE[1]
orde <- ORDER[1]
simu <- 1

toplot <- subset(data3,site==sit & order == orde)
plot(toplot$simul,toplot$Cpltarity,type="l")#,ylim=c(-10000,300))
# abline(toplot$simul,toplot$Selection)
