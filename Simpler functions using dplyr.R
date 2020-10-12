source("R/Common variables.R")
# simpler temporal_plot function using dplyr ####
site <- "Bever"
order <- "decreasing"
number <- 5
number2 <- 6
res<-try(read.table(paste0("data/raw/output-cmd2_",site,"_",order,".txt/forceps.",site,".site_",number,"_complete.txt")),silent=T) 
res2 <- try(read.table(paste0("data/raw/output-cmd2_",site,"_",order,".txt/forceps.",site,".site_",number+1,"_complete.txt")),silent=T) 

res3 <- try(read.table(paste0("data/raw/output-cmd2_",site,"_",order,"_2.txt/forceps.",site,".site_",number,"_complete.txt")),silent=T) 
res4 <- try(read.table(paste0("data/raw/output-cmd2_",site,"_",order,"_2.txt/forceps.",site,".site_",number+1,"_complete.txt")),silent=T) 


biomfunct <- function(res,number){
  # res is the raw data frame. Here, you choose if you want decreasing, increasing, random, with temperature increase, or not.
  # number is the number of the simulation (ex decreasing, simul 5)
  # NB NB NB if I launch all the simulations in one for the 30 random orders, 
  # they will be ranked from 1 to 900, i.e. 30 orders with are all ranked from 1 to 30.
  # I will have to rename them: from 1 to 30: order "random_1" and simul 1 to 30. Then from 31 to 60: order="random_2", etc.
  # gives as output: abundance  biomass  site  order  simul  mixture.t.ha
  colnames(res) <- colnames_res
  res2 <- res %>%
    group_by(date,speciesShortName) %>%
    mutate (biomass=sum(biomass.kg.)) %>%
    mutate(ab=1)%>%
    mutate(abundance=sum(ab))
  
  temp_plot <- summarise(res2, abundance=mean(abundance),biomass=mean(biomass))
  # compare it to:
  # old <- temporal_plot(res)
  
  # and put a threshold on it ####
  temp_plot_threshold <- temp_plot %>%
    mutate(biom_tot=sum(biomass))%>%
    subset(biomass>0.001*biom_tot)
  # compare to:
  # old2 <- temporal_plot_threshold(old)
  
  # Both give the same results.
  
  # Biomass ####
  dates <- as.numeric(max(temp_plot_threshold$date))- c(900,800,700,600,500,400,300,200,100,0) # years on which we average the biomass
  years_to_keep <- subset(temp_plot_threshold,date %in%dates)
  output <- years_to_keep %>%
    ungroup()%>%
    group_by(speciesShortName) %>%
    mutate(abundance = mean(abundance))%>%
    mutate(biomass=mean(biomass)) %>%
    summarize(abundance=mean(abundance),biomass=mean(biomass)) %>%
    mutate(site=site) %>%
    mutate(order=order) %>%
    mutate(simul=number) %>%
    mutate(mixture.t.ha = biomass/(1000*0.08*Nbpatches)) %>% 
    mutate(mixture_relative = mixture.t.ha/sum(mixture.t.ha) )
  
  # Here again, gives the same results as my code. 
  # --> I could have a more efficient code, and I could have earned a lot of time of coding...
  # --> At least, it gives me confidence in the robustness of my code
  output
}

 
r1 <- biomfunct(res,number)
r2 <- biomfunct(res2,number2)

r3 <- biomfunct(res3,number)
r4 <- biomfunct(res4,number2)

plot(r1$mixture.t.ha~r1$speciesShortName)
plot(r2$mixture.t.ha~r2$speciesShortName)

sum(r1$mixture.t.ha)
sum(r2$mixture.t.ha)

sum(r3$mixture.t.ha)
sum(r4$mixture.t.ha)


