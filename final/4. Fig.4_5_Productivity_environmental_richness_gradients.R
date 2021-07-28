# packages ####
source("Final/0. Packages.R")

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
    # ggtitle(site) ,
    # theme(plot.title = element_text(size=10)) ,
    scale_x_continuous(breaks = 10*c(1:3)) ,
    theme_light(),
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
    egg::theme_article(base_size=18)
  
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
measure <- "productivity_tot"
datatoplot <- read.table(paste0("data/processed/",measure,"_with interval_median.txt"),header=T)
# Environmental conditions
sites <- read.table("data/raw/Site description.txt",header=T)
sites$number <- c(5,6,4,7,8,9,3,2,10,1,11) # number the sites in the environmental graph
sites_order_plot <- sites %>% 
  arrange(number) %>% 
  mutate(group = c("cold","cold","cold","warm-wet","warm-wet","warm","warm","warm","warm-dry","warm-dry","warm-dry"))
sites_order_plot$group <- factor(sites_order_plot$group ,levels=c("cold","warm-wet","warm","warm-dry"))

# Fig. 4 ####
# all the plots
PLOT <- list()
i <- 0
for (sit in sites_order_plot$Site){
  nb <- sites_order_plot %>% 
    filter(Site==sit) %>% 
    pull(number)
  i <- i+1
  result <- 
    datatoplot %>% 
    filter(site==sit)
  PLOT[[i]] <- ggplot(result,aes(x=simul-1,y=decreasing)) +  
    gg_removal_productivity()+
    ggtitle(paste0(nb,". ",sit))
}

legend <- 
  read.table(paste0("data/processed/productivity_tot_with interval_median.txt"),header=T) %>% 
  filter(site=="Bever") %>% 
  extract_legend() 

environment_number <-
  ggplot(sites_order_plot, aes(x=Temp_moy,y=Annual_ppt,label=number,shape=group))+
  geom_point()+
  scale_shape_manual(values = c(4, 1, 2, 8), name="Environment") +
  xlab("Mean annual temperature (°C)")+
  ylab("Annual precipitation (mm)") +
  xlim(1,10) +
  ylim(380,1500)+
  ggrepel::geom_text_repel() +
  egg::theme_article(base_size=10)
  

complete_plot3 <- grid.arrange( PLOT[[1]],PLOT[[2]], PLOT[[3]],
     PLOT[[4]], PLOT[[5]],PLOT[[6]],
     PLOT[[7]],PLOT[[8]],PLOT[[9]],
     PLOT[[10]], PLOT[[11]],legend,
     ncol = 3, nrow = 4, 
     layout_matrix = rbind( c(0,1,2),c(3,4,5),c(6,7,8),
                           c(9,10,11)))

complete_plot4 <- annotate_figure(complete_plot3,
                # top = text_grob("Visualizing len", color = "red", face = "bold", size = 14),
                # bottom = text_grob("Data source: \n ToothGrowth data set", color = "blue",
                #                    hjust = 1, x = 1, face = "italic", size = 10),
                bottom = text_grob("Number of species lost", rot = 0,size=15),
                left = text_grob("Ecosystem productivity (t/ha)",rot=90,size=15),
                # fig.lab = "Figure 1", fig.lab.face = "bold"
)

width_A <- 1
height_A <- .3
p <- ggdraw() +
  draw_plot(environment_number, x = (1-width_A)/2, y = 1-height_A, width = width_A, height = height_A) +
  draw_plot(complete_plot4, x = 0, y = 0, width = 1, height = 1-height_A) +
  draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0, 0), y = c(1, 1-height_A))

ggsave(filename = paste0("figures_tables/Main figure_",measure,"_removal.png"), 
       plot = p,
       width = 20, 
       height = 30,
       units = "cm",
       dpi = 300)




# Fig. 5 ####
data_AUC1 <- datatoplot %>% 
  group_by(site) %>% 
  summarize(decreasing=sum(decreasing)/max(decreasing),increasing=sum(increasing)/max(increasing)) %>%
  arrange(match(site,sites_order_plot$Site)) %>%
  mutate(group=sites_order_plot$group) %>% 
  mutate(number=sites_order_plot$number) %>% 
  mutate(number=c(1:11))  %>% 
  mutate(Temp_moy = sites_order_plot$Temp_moy)

# sites2 <- sites %>%
#   column_to_rownames("Site") %>% 
#   select(Temp_moy,Annual_ppt,Na.kg.ha.an.,ProdMax.T.ha.an.)
# PCA <- FactoMineR::PCA(sites2)
# coord_sites <- as.data.frame(PCA$ind$coord) %>% 
#   rownames_to_column() %>% 
#   rename(site=rowname) %>% 
#   arrange(match(site,sites_order_plot$Site)) %>% 
#   mutate(group=data_AUC1$group)
#   
# 
# site_var <- as.data.frame(PCA$var$coord)%>% 
#   rownames_to_column() %>% 
#   rename(variable=rowname) 

order_loss <- function(orde){
  if (orde=="decreasing"){
    "A. Distinct species lost first"
  } else {
    "B. Common species lost first"
  }
}

data_AUC2 <- data_AUC1 %>% 
  gather(order,AUC,-site,-group,-number,-Temp_moy) %>% 
  mutate(order = map_chr(order,order_loss))
data_AUC2$group <- factor(data_AUC2$group,levels=c("cold","warm-wet","warm","warm-dry"))
data_AUC2$order <- factor(data_AUC2$order,levels=c("A. Distinct species lost first","B. Common species lost first"))


plot_AUC <- ggplot(data_AUC2,aes(x=Temp_moy,y=AUC,label=number))+
  facet_wrap(~order)+
  geom_point()+
  scale_shape_manual(values = c(4, 1, 2, 8),name="Environment") +
  labs(color="Order of species loss") +
  egg::theme_article(base_size=11) +
  ggrepel::geom_text_repel() +
  xlab("Mean annual temperature (°C)") +
  ylab("Relative area under the curve (AUC)") +
  theme(strip.text = element_text(size=11))

ggsave(filename = paste0("figures_tables/Main figure 2_AUC.png"), 
       plot = plot_AUC,
       width = 17, 
       height = 10,
       units = "cm",
       dpi = 300)