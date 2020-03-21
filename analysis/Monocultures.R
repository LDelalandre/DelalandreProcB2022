library("dplyr")
library("magrittr")
library("ggplot2")
colnames_mean<-colnames(read.table("data/colnames_mean.txt",header=T)) # idem
colnames_res<-colnames(read.table("data/colnames_res.txt", header=T))

biomass<-c()
sd_biom<-c()
for(i in c(1:30)){
  mean <- read.table(paste0("data/raw/output-cmd2_Bern_monocultures.txt/forceps.Bern.site_",i,"_mean.txt"))
  colnames(mean)<-colnames_mean
  # biomass averaged
  years_to_keep <- max(mean$date) - c(900,800,700,600,500,400,300,200,100,0)
  meanbiom <- mean(subset(mean,date %in% years_to_keep)$totalBiomass.t.ha.)
  biomass<-c(biomass,meanbiom)
  sd_biom<-c(sd_biom,sd(mean$totalBiomass.t.ha.))
}


# specific_values: a data frame with final biomass of each monoculture, etc.
specific_values <- read.table("data/raw/distinctiveness of the species.txt",header=T)
specific_values$biomass <- biomass
specific_values$sd_biom <- sd_biom
specific_values$Id <- c(0:29)
tot_biomass <- sum(specific_values$biomass)

# Specific biomass f(distinctiveness)
ggplot(specific_values,aes(x=Di,y=biomass,label=SName))+
  geom_point() +
  geom_label()



DIST_ORDER <- list(dist[order(specific_values$Di,decreasing=T),] , dist[order(specific_values$Di,decreasing=F),])
dist_order <- DIST_ORDER[1]
removal_biomass <- c()
for (i in 1:30){
  B <- specific_values %>%
    subset(Id %in% dist_ordered[c(i:30),]$Id)%>% # select species in
    select(biomass)
  B <- sum(B)/length(c(i:30))
  # B <- sum(B)/ tot_biomass
  # B <- sum(B)
  removal_biomass <- c(removal_biomass,B)
}
plot(c(1:30),removal_biomass)


