library("dplyr")
library("magrittr")
library("ggplot2")
source("R/Analysis_data.R")
colnames_mean<-colnames(read.table("data/colnames_mean.txt",header=T)) # idem
colnames_res<-colnames(read.table("data/colnames_res.txt", header=T))

biomass<-c()
sd_biom<-c()
for(i in c(1:30)){
  mean <- read.table(paste0("data/raw/output-cmd2_Bern_monocultures.txt/forceps.Bern.site_",i,"_mean.txt"))
  # NB it is in the good order because we made the monocultures in the order of the species Id number.
  colnames(mean)<-colnames_mean
  # biomass averaged
  years_to_keep <- max(mean$date) - c(900,800,700,600,500,400,300,200,100,0)
  meanbiom <- mean(subset(mean,date %in% years_to_keep)$totalBiomass.t.ha.)
  biomass<-c(biomass,meanbiom)
  sd_biom<-c(sd_biom,sd(mean$totalBiomass.t.ha.))
}


# specific_values: a data frame with final biomass of each monoculture, etc. ####
specific_values <- read.table("data/raw/distinctiveness of the species.txt",header=T)
specific_values$biomass_monoculture <- biomass
specific_values$sd_biom <- sd_biom
specific_values$Id <- c(0:29)
specific_values$relative_biomass <- specific_values$biomass / sum(specific_values$biomass)


# Compare specific biomass in monocultures and in the complete mixture ####
res <- read.table("data/raw/output-cmd2_Bern_decreasing.txt/forceps.Bern.site_1_complete.txt")
colnames(res)<-colnames_res
temp_plot <- temporal_plot(res)
# années sur lesquelles on veut moyenner les résultats.
# NB : ici, je ne les ai pas dans mon export. A refaire de manière à homogénéiser tout ça. En attendant :
dates <- as.numeric(max(temp_plot$date))#- c(900,800,700,600,500,400,300,200,100,0)
years_to_keep <- subset(temp_plot,date %in%dates) # we compute a mean value of biomass on those years
biomasses <- aggregate(years_to_keep$biomass, list(years_to_keep$species), mean)
colnames(biomasses) <- c("species","mixture")
biomasses$mixture <- biomasses$mixture/800 # so that the unit becomes t/ha
biomasses$mixture <- biomasses$mixture/sum(biomasses$mixture) #to have relative biomass

biom_mono <- c()
dist <- c()
for (sp in biomasses$species){
  biom_mono <- c(biom_mono,specific_values[which(specific_values$SName == sp ),]$biomass)
  dist <- c(dist,specific_values[which(specific_values$SName == sp ),]$Di)
}
biomasses$monoculture <- biom_mono/sum(biom_mono)
biomasses$Di <- dist

ggplot(biomasses,aes(x=monoculture,y=mixture,label=species))+
  geom_point() +
  geom_label()


# Specific monoculture biomass f(distinctiveness) ####
ggplot(specific_values,aes(x=Di,y=biomass,label=SName))+
  geom_point() +
  geom_label()

# Mixture and monoculture specific biomasses f(distinctiveness) ####
ggplot(biomasses,aes(y=Di,x=monoculture,label=species))+
  geom_point() +
  labs(x="Relative biomass",y="Distinctiveness") +
  ggtitle("Monocultures") +
  geom_label()
ggplot(biomasses,aes(y=Di,x=mixture,label=species))+
  geom_point() +
  labs(x="Relative biomass",y="Distinctiveness") +
  ggtitle("Mixture") +
  # xlim(0, 0.11) +
  geom_label()


# Species removal ####
DIST_ORDER <- list(specific_values[order(specific_values$Di,decreasing=T),]$Id, specific_values[order(specific_values$Di,decreasing=F),]$Id)
dist_order <- DIST_ORDER[[2]]
removal_biomass <- c()
for (i in 1:30){
  B <- specific_values %>%
    subset(Id %in% dist_order[c(i:30)])%>% # select species in
    select(biomass)
  B <- sum(B)/length(c(i:30))
  # B <- sum(B)/ tot_biomass
  # B <- sum(B)
  removal_biomass <- c(removal_biomass,B)
}
plot(c(1:30),removal_biomass)

# Search for and explanation to the shape of the curves ####
x=c(1:30)

y1=4500-150*x # linear decrease
plot(x,y1)

y2=4500-5*x^2 # concave-down decrease
plot(x,y2)

y3=4500-821.5*x^0.5 # Concave_up decrease (fait tq à x=0, on ait y=0, soit 4500/(30^0.5) = coef directeur)
plot(x,y3)

z1=rep(0,30)
z2 <- z1
z3 <- z1
for (i in 1:30){
  z1[i]=y1[i]/length(i:30)
  z2[i]=y2[i]/length(i:30)
  z3[i]=y3[i]/length(i:30)
}

plot(x,y1) ; points(x,y2) ; points(x,y3)
plot(x,z1) 
plot(x,z2) 
plot(x,z3)

