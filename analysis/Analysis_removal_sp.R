library(ggplot2)
library("funrar")
library(forcats) # for reordering factors
library("dplyr")
library("magrittr")
source("R/Analysis_data.R")
# Plot specific biomasses of the species at the end of the simulations.

colnames_mean<-colnames(read.table("data/colnames_mean.txt",header=T)) # idem
colnames_res<-colnames(read.table("data/colnames_res.txt", header=T))

SITE <- c("Bern","Bever","Cottbus","Huttwil")

ORDER <- c("increasing","decreasing","random_1","random_2") # I take two random simulations as examples.
# NB: I could average all biomasses amongst all random simulations.


# Plot: abundance of each species at the end of each simul ####
for (site in SITE){
  for (order in ORDER){
    res<-read.table(paste0("data/raw/output-cmd2_",site,"_",order,".txt/forceps.",site,".site_",1,"_complete.txt"))
    colnames(res)<-colnames_res
    temp_plot <- temporal_plot(res)
    FINAL <- filter(temp_plot,date==max(unique(res$date))) # final year
    FINAL$simul <- rep(1,dim(FINAL)[1])
    
    for(i in c(2:30)){
      res<-try(read.table(paste0("data/raw/output-cmd2_",site,"_",order,".txt/forceps.",site,".site_",i,"_complete.txt")), silent=T)
      if (class(res) != "try-error"){
        colnames(res)<-colnames_res
        temp_plot <- temporal_plot(res)
        final <- filter(temp_plot,date==max(unique(res$date))) # final year
        final$simul <- rep(i,dim(final)[1])
        FINAL <- rbind(FINAL,final)
      }
    }
    
    
    plot <- ggplot(FINAL,aes(x=simul,y=biomass,label=species))+
      labs(x="Number of species removed",y="Specific biomass") +
      geom_point() +
      geom_label() +
      aes(color=species) +
      if ( order == "increasing") {
        ggtitle(site,"Removing distinct species last")
      } else if ( order == "decreasing") {
        ggtitle(site,"Removing distinct species first")
      } else {
        ggtitle(site,"Random order of removal of species")
      } 
    plot + ggsave(paste0("figures/specific_biomass_",site,"_",order,".png"))
    
  }
}
