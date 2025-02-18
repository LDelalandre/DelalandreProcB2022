library(gridExtra)
source("R/Common variables.R")
source("R/Monocultures_functions.R")

# site <- SITE[1]
# order <- ORDER[2]
# number <- 1 # number of the simulation. 1: All the species are present. 30: Just one species is left.

# Superposition of mixture and monoculture ####
for (site in SITE){
  for (order in ORDER[1:2]){
    removal_comparison <- data.frame(matrix(data=c(0:29),ncol=3,nrow=30,dimnames=list(c(1:30),c("Nb_sp_removed","Monoculture_biomass","Mixture_biomass"))))
    tot <- read.table(paste0("data/processed/biomass_specific_",site,"_with monocultures.txt"),header=T)
    for (number in c(1:30)){
      ord <- order
      biomass_comp <- subset(tot,order==ord & simul==number)
      removal_comparison[number,2] <- sum(biomass_comp$monoculture.t.ha)/dim(biomass_comp)[1]
      removal_comparison[number,3] <- sum(biomass_comp$mixture.t.ha)
    }
    
    write.table(removal_comparison,paste0("data/processed/removal_mixt_mono_",site,"_",order,".txt"),row.names=FALSE)

    plot = ggplot(removal_comparison,aes(x = Nb_sp_removed))+
      geom_line(aes(y = Monoculture_biomass, colour = "Monoculture")) +
      geom_line(aes(y = Mixture_biomass, colour = "Mixture"))+
      labs(x="Number of species removed",y="Biomass")+
      ggtitle(paste0(site,"_",order))
    ggsave(filename = paste0("figures/Monocultures_",site,"/removal_mixt_mono_",site,"_",order,".png"),
           plot=plot)
  }
}

# Productivity
for (site in SITE){
  for (order in ORDER[1:2]){
    removal_comparison <- data.frame(matrix(data=c(0:29),ncol=3,nrow=30,dimnames=list(c(1:30),c("Nb_sp_removed","Monoculture_biomass","Mixture_biomass"))))
    tot <- read.table(paste0("data/processed/productivity_specific_",site,"_with monocultures.txt"),header=T)
    for (number in c(1:30)){
      ord <- order
      biomass_comp <- subset(tot,order==ord & simul==number)
      removal_comparison[number,2] <- sum(biomass_comp$monoculture)/dim(biomass_comp)[1]
      removal_comparison[number,3] <- sum(biomass_comp$mixture_t_ha)
    }
    
    write.table(removal_comparison,paste0("data/processed/prod_removal_mixt_mono_",site,"_",order,".txt"),row.names=FALSE)
    
    
    plot = ggplot(removal_comparison,aes(x = Nb_sp_removed))+
      geom_line(aes(y = Monoculture_biomass, colour = "Monoculture")) + 
      geom_line(aes(y = Mixture_biomass, colour = "Mixture"))+
      labs(x="Number of species removed",y="Productivity")+
      ggtitle(paste0(site,"_",order)) +
      scale_x_continuous(breaks = 2*c(1:15)) +
      theme(legend.title = element_blank()) 
    if (order=="increasing"){
      plot <- plot + scale_color_manual(values=c("#00BFC4","green4"))
    } else {
      plot <- plot + scale_color_manual(values=c("#F8766D","green4"))
    }
    ggsave(filename = paste0("figures/Monocultures_",site,"/prod_removal_mixt_mono_",site,"_",order,".png"),
           plot=plot)
  }
}


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


