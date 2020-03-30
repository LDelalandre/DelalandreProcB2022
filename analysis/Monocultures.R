library("dplyr")
library("magrittr")
library("ggplot2")
library(gridExtra)
source("R/Analysis_data.R")
source("R/Monocultures_functions.R")


site <- SITE[1]
order <- ORDER[2]
number <- 1 # number of the simulation. 1: All the species are present. 30: Just one species is left.

# Specific values of biomass and sd(biomass) for the monocultures ####
specific_val <- specific_values(site)

# Superposition of mixture and monoculture ####
removal_comparison <- data.frame(matrix(data=c(0:29),ncol=3,nrow=30,dimnames=list(c(1:30),c("Nb_sp_removed","Monoculture_biomass","Mixture_biomass"))))
for (number in c(1:30)){
  biomass_comp <- biomasses(site=site,specific_val=specific_val,number=number,order=order)
  removal_comparison[number,2] <- sum(biomass_comp$`monoculture(t/ha)`)/dim(biomass_comp)[1]
  removal_comparison[number,3] <- sum(biomass_comp$`mixture(t/ha)`)
}
plot = ggplot(removal_comparison,aes(x = Nb_sp_removed))+
  geom_line(aes(y = Monoculture_biomass, colour = "Monoculture")) + 
  geom_line(aes(y = Mixture_biomass, colour = "Mixture"))+
  labs(x="Number of species removed",y="Biomass")+
  ggtitle(paste0(site,"_",order))
ggsave(filename = paste0("figures/removal_mixt_mono_",site,"_",order,".png"),
       plot=plot)


# Species removal ####
DIST_ORDER <- list(specific_val[order(specific_val$Di,decreasing=T),]$Id, specific_val[order(specific_val$Di,decreasing=F),]$Id)
dist_order <- DIST_ORDER[[2]]
removal_biomass <- c()
for (i in 1:30){
  B <- specific_val %>%
    subset(Id %in% dist_order[c(i:30)])%>% # select species in
    select(biomass_monoculture)
  B <- sum(B)/length(c(i:30))
  # B <- sum(B)/ tot_biomass
  # B <- sum(B)
  removal_biomass <- c(removal_biomass,B)
}
plot(c(1:30),removal_biomass)


  


# Mixture and monoculture specific biomasses f(distinctiveness) ####
for (number in c(1:30)){
  # Comparison of biomasses between monocultures and mixtures
  biomass_comp <- biomasses(site=site,specific_val=specific_val,number=number,order=order)
  
  g1 <- ggplot(biomass_comp,aes(y=Di,x=monoculture_relative,label=species))+
    geom_point() +
    labs(x="Relative biomass",y="Distinctiveness") +
    ggtitle("Monocultures") +
    ylim(min(biomass_comp$Di),max(biomass_comp$Di)) +
    geom_label()
  g2 <- ggplot(biomass_comp,aes(y=Di,x=mixture_relative,label=species))+
    geom_point() +
    labs(x="Relative biomass",y="Distinctiveness") +
    ggtitle("Mixture") +
    ylim(min(biomass_comp$Di),max(biomass_comp$Di)) +
    geom_label()
  
  ggsave(filename = paste0("figures/Monocultures_",site,"/specific_biomass_mixt_mono_",site,"_",order,"_",number,".png"),
         plot=grid.arrange(g1,g2))
}

for (number in c(1:30)){
  # Comparison of biomasses between monocultures and mixtures
  biomass_comp <- biomasses(site=site,specific_val=specific_val,number=number,order=order)
  # Relative decrease in biomass when the species is put in mixture
  biomass_comp$relat_decr <- (biomass_comp$`monoculture(t/ha)`-biomass_comp$`mixture(t/ha)` )/biomass_comp$`monoculture(t/ha)`
}

plot(biomass_comp$relat_decr~biomass_comp$monoculture_relative)
plot(biomass_comp$relat_decr~biomass_comp$mixture_relative)

hist(biomass_comp$relat_decr)
subset(biomass_comp,relat_decr<mean(relat_decr))


# Compare specific biomass in monocultures and in the complete mixture ####
ggplot(biomass_comp,aes(x=monoculture_relative,y=mixture_relative,label=species))+
  geom_point() +
  geom_label()


# Specific monoculture biomass f(distinctiveness) ####
ggplot(specific_val,aes(x=Di,y=biomass_monoculture,label=SName))+
  geom_point() +
  geom_label()



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

