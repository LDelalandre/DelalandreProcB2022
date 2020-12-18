# packages ####
library(ggplot2)

# scripts ####
source("R/Common variables.R")
source("R/Before simulations.R")

# functions ####
oneplot <- function(df,legend){
  # Plot of species removal expe in one site.
  # Columns of the data frame df:  site simul decreasing increasing random_1 random_10 etc.
  # Legend: a boolean. Display a legend or not.
    ggplot(df,aes(x=simul-1,y=decreasing)) +
    labs(x="Number of species removed",y=measure) +
    geom_line(aes(x=simul-1,y=increasing, color="#00BFC4"),size=2) +
    geom_line(aes(x=simul-1,y=decreasing, color="#F8766D"),size=2) +
    geom_ribbon(aes(ymin=int_min, ymax=int_max), alpha=0.5,fill="grey60") +
    geom_line(aes(x=simul-1,y=int_max, color="grey60"),size=0) +
    geom_line(aes(x=simul-1,y=int_min, color="grey60"),size=0) +
    ggtitle(site) +
    theme(plot.title = element_text(size=24)) +
    scale_x_continuous(breaks = 5*c(1:6)) +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none" ) + # virer tous les titres
    theme(axis.text=element_text(size=20)) +
    scale_color_identity(name = "Order of species loss",
                         breaks = c("#00BFC4","#F8766D","grey60"),
                         labels = c("Common species lost first", "Distinct species lost first", "Species lost randomly"),
                         guide = "legend")+
    {if(legend==T)theme(legend.position = "right")}
  # theme(legend.title = element_blank())
}

extract_legend <- function(result){
  # Extract the legend alone, from the data frame of species removal expe
  plot <- ggplot(result,aes(x=simul-1,y=decreasing)) + 
    gg_removal() + 
    theme(legend.position = "right")
  leg <- ggpubr::get_legend(plot)
  ggpubr::as_ggplot(leg)
}


# Extract the legend
legend <- read.table(paste0("data/processed/productivity_tot_Bever_with interval_median.txt"),header=T) %>% 
  extract_legend()

# generate plots ####
measure <- "productivity_tot"

PLOT <- list()
i <- 0
for (site in ord_plots){
  i <- i+1
  result <- read.table(paste0("data/processed/",measure,"_",site,"_with interval_median.txt"),header=T)
  PLOT[[i]] <- oneplot(result,legend=F)
}
PLOT

gg_removal <- function(){
  list(
    labs(x="Number of species removed",y=measure),
    geom_line(aes(x=simul-1,y=increasing, color="#00BFC4"),size=2),
    geom_line(aes(x=simul-1,y=decreasing, color="#F8766D"),size=2) ,
    geom_ribbon(aes(ymin=int_min, ymax=int_max), alpha=0.5,fill="grey60") ,
    geom_line(aes(x=simul-1,y=int_max, color="grey60"),size=0) ,
    geom_line(aes(x=simul-1,y=int_min, color="grey60"),size=0) ,
    ggtitle(site) ,
    theme(plot.title = element_text(size=24)) ,
    scale_x_continuous(breaks = 5*c(1:6)) ,
    theme(legend.position = "none" ) ,
    xlab("Number of species removed"),
    ylab("Productivity (t/ha)"),
    theme(axis.text=element_text(size=20)) ,
    scale_color_identity(name = "Order of species loss",
                           breaks = c("#00BFC4","#F8766D","grey60"),
                           labels = c("Common species lost first", "Distinct species lost first", "Species lost randomly"),
                           guide = "legend")
  ) 
}


ggplot(result,aes(x=simul-1,y=decreasing)) +
  gg_removal() +
  theme(axis.text=element_text(size=12),
         axis.title=element_text(size=14,face="bold")) 
  theme(legend.position = "right") 
  
  
