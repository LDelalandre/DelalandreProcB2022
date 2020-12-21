# packages ####
library(ggplot2)
library(cowplot)

# scripts ####
source("R/Common variables.R")
source("R/Before simulations.R")

# functions ####
gg_removal <- function(){
  # Plot of species removal expe in one site.
  # Columns of the data frame df:  site simul decreasing increasing random_1 random_10 etc.
  # Legend: a boolean. Display a legend or not.
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
    theme(axis.title = element_blank()),
    scale_color_identity(name = "Order of species loss",
                         breaks = c("#00BFC4","#F8766D","grey60"),
                         labels = c("Common species lost first", "Distinct species lost first", "Species lost randomly"),
                         guide = "legend")
  ) 
}

gg_removal_small <- function(){
  # Plot of species removal expe in one site.
  # Columns of the data frame df:  site simul decreasing increasing random_1 random_10 etc.
  # Legend: a boolean. Display a legend or not.
  list(
    labs(x="Number of species removed",y=measure),
    geom_line(aes(x=simul-1,y=increasing, color="#00BFC4"),size=1),
    geom_line(aes(x=simul-1,y=decreasing, color="#F8766D"),size=1) ,
    geom_ribbon(aes(ymin=int_min, ymax=int_max), alpha=0.5,fill="grey60") ,
    geom_line(aes(x=simul-1,y=int_max, color="grey60"),size=0) ,
    geom_line(aes(x=simul-1,y=int_min, color="grey60"),size=0) ,
    ggtitle(site) ,
    # theme(plot.title = element_text(size=10)) ,
    scale_x_continuous(breaks = 10*c(1:3)) ,
    ylim(0,3.5),
    theme(legend.position = "none" ) ,
    xlab("Number of species removed"),
    ylab("Productivity (t/ha)"),
    theme(axis.text=element_text(size=10)) ,
    theme(axis.title = element_blank()),
    scale_color_identity(name = "Order of species loss",
                         breaks = c("#00BFC4","#F8766D","grey60"),
                         labels = c("Common species lost first", "Distinct species lost first", "Species lost randomly"),
                         guide = "legend")
  ) 
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

# Environmental conditions plot ####
sites <- read.table("data/Site description.txt",header=T)

COORD <- sites %>% 
  # arrange(ord_plots) %>% 
  select(Site,Temp_moy,Annual_ppt) %>% 
  arrange(Temp_moy,Annual_ppt)


# generate plots ####
measure <- "productivity_tot"

PLOT <- list()
i <- 0
for (site in COORD$Site){
  i <- i+1
  result <- read.table(paste0("data/processed/",measure,"_",site,"_with interval_median.txt"),header=T)
  PLOT[[i]] <- ggplot(result,aes(x=simul-1,y=decreasing)) +  gg_removal_small()
}
PLOT



# Size of the graphs inserted in the environmental plot
width <- 2.3 # width of the subplots on the x axis
height <- 200

# Position of the graphs in the environmental plot
# Each site belongs to a group, defined as follows:
# group 1: graph at the top of the point
# group 2: graph at the bottom of the point
# group 3: graph on the left of the point
# group 4: graph on the right of the point
group <- c(1,2,4,3,1,3,1,3,2,4,4) # group of each site

# Environmental plot
environment2 <-
  ggplot(sites, aes(x=Temp_moy,y=Annual_ppt,label=Site))+
  geom_point()+
  # geom_label()+
  # geom_text(fontface = "bold",position=position_jitter(width=1,height=1)) +
  xlab("Average temperature (°C)")+
  ylab("Annual precipitations (mm)") +
  xlim(0,12) +
  ylim(400,1500)+
  theme(axis.text=element_text(size=15))+
  theme(axis.title = element_text(size=15)) +
  theme_article()

# Fill the environmental plot with the graphs of the sites
environment <- environment2
for (i in 1:11){
  coord <- COORD[i,-1]
  X <- as.numeric(coord[1])
  Y <-  as.numeric(coord[2])
  if (group[i]==1){
    environment <- environment +
      annotation_custom(ggplotGrob(PLOT[[i]]), 
                        xmin = X-width/2, xmax = X+width/2, ymin = Y+10, ymax = Y+height+10)  
  } else if(group[i]==2){
    environment <- environment +
      annotation_custom(ggplotGrob(PLOT[[i]]), 
                        xmin = X-width/2, xmax = X+width/2, ymin = Y-height-10, ymax = Y-10)
  } else if(group[i]==3){
    environment <- environment +
      annotation_custom(ggplotGrob(PLOT[[i]]), 
                        xmin = X-width-0.2, xmax = X-0.2, ymin = Y-height/2, ymax = Y+height/2)
  } else if(group[i]==4){
    environment <- environment +
      annotation_custom(ggplotGrob(PLOT[[i]]), 
                        xmin = X+0.2, xmax = X+width+0.2, ymin = Y-height/2, ymax = Y+height/2)
  } else if(group[i]==5){
    environment <- environment +
      annotation_custom(ggplotGrob(PLOT[[i]]), 
                        xmin = X+0.2, xmax = X+width+0.2, ymin = Y, ymax = Y+height)
  }
  
}
# Add a legend
environment + annotation_custom(ggplotGrob(legend), 
                                xmin = 0, xmax = 2.5, ymin = 400, ymax = 600) 




# REMINDER ####
# My plot is modular: I can display the legend or not, display the name of the axis or nor.
# Here is the code:
ggplot(result,aes(x=simul-1,y=decreasing)) +
  gg_removal() +
  theme(axis.title=element_text(size=20,face="bold")) 
  theme(legend.position = "right") 
  
# tuto à faire https://www.r-spatial.org/r/2018/10/25/ggplot2-sf-3.html
# celui-là pour commencer: https://aosmith.rbind.io/2019/04/22/embedding-subplots/
  
