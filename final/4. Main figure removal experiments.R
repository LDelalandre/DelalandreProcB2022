# packages ####
library(ggplot2)
library(cowplot)
library(egg) # for theme_article

# scripts ####
source("R/Common variables.R")
source("R/Before simulations.R")

# functions ####
gg_removal_productivity <- function(){
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
    theme(legend.position = "none" ) ,
    xlab("Number of species removed"),
    ylab("Productivity (t/ha)"),
    theme(axis.text=element_text(size=10)) ,
    theme(axis.title = element_blank()),
    scale_color_identity(name = "Order of species loss",
                         breaks = c("#00BFC4","#F8766D","grey60"),
                         labels = c("Common species lost first", "Distinct species lost first", "Species lost randomly"),
                         guide = "legend"),
    if(measure=="productivity_tot"){
      ylim(0,3.5)
    }
  ) 
}


extract_legend <- function(result){
  # Extract the legend alone, from the data frame of species removal expe
  plot <- ggplot(result,aes(x=simul-1,y=decreasing)) + 
    gg_removal_productivity() + 
    theme(legend.position = "right") +
    theme(legend.title = element_text(size=14), #change legend title font size
          legend.text = element_text(size=10) #change legend text font size)
    )
  
  leg <- ggpubr::get_legend(plot)
  ggpubr::as_ggplot(leg)
}

environmental_plot <- function(PLOT,COORD){
  # Generates the final figure.
  # PLOT is a list with the 11 removal experiment graphs to insert in the envt plot.
  # COORD is a data frame with site names and the coordinates of the sites (T° and precipitations)
  #_______________________________________________________________________________
  # Coordinates and position ####
  
  # Size of the graphs inserted in the environmental plot
  width <- 2.3 # width of the subplots on the x axis
  height <- 240
  
  # Position of the graphs in the environmental plot
  # Each site belongs to a group, defined as follows:
  # group 1: graph at the top of the point
  # group 2: graph at the bottom of the point
  # group 3: graph on the left of the point
  # group 4: graph on the right of the point
  # group 5: Basel, on the right and a bit on the top
  group <- c(1,2,4,3,1,3,1,3,2,5,4) # group of each site
  
  #_______________________________________________________________________________
  # Extract the legend
  legend <- 
    read.table(paste0("data/processed/productivity_tot_with interval_median.txt"),header=T) %>% 
    filter(site=="Bever") %>% 
    extract_legend()
  
  #_______________________________________________________________________________
  # Environmental plot ####
  environment2 <-
    ggplot(sites, aes(x=Temp_moy,y=Annual_ppt,label=Site))+
    geom_point()+
    xlab("Average temperature (°C)")+
    ylab("Annual precipitations (mm)") +
    xlim(0,12) +
    ylim(380,1500)+
    theme_article(base_size=18)
  
  #_______________________________________________________________________________
  # Insert the graphs in the environmental plot ####
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
                          xmin = X+0.2, xmax = X+width+0.2, ymin = Y-height/3, ymax = Y+2*height/3)
    }
    
  }
  # Add a legend
  environment_final <- environment + annotation_custom(ggplotGrob(legend), 
                                                       xmin = 0, xmax = 4, ymin = 350, ymax = 600) 
  # save the plot
  ggsave(filename = paste0("paper_2/",measure,"_removal.png"), 
         plot = environment_final,
         width = 30, 
         height = 26,
         units = "cm",
         dpi = 300)
  
  
}

#_______________________________________________________________________________
# CODE ####

# Environmental conditions
sites <- read.table("data/Site description.txt",header=T)

# Sites and coordinates used hereafter
COORD <- sites %>% 
  select(Site,Temp_moy,Annual_ppt) %>% 
  arrange(Temp_moy,Annual_ppt)

for (measure in c("productivity_tot","TS_productivity_tot_unfiltered","TS_productivity_tot_filtered")){
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
      ggtitle(sit)
  }
  
  # Add axes labels for one site (Bever chosen here)
  result <- 
    datatoplot %>% 
    filter(site=="Bever")
  PLOT[[2]] <- ggplot(result,aes(x=simul-1,y=decreasing)) + 
    gg_removal_productivity() +
    theme(axis.title=element_text()) +
    ggtitle("Bever") +
    if (measure != "productivity_tot"){
      ylab("Temporal stability \n of poductivity")
    }
  
  # generate the final graph
  environmental_plot(PLOT,COORD)
}


