# Price analysis ####

table_threshold_price <- function(site,order){
  # Writes a table with the final biomass and abundance of each species (above a biomass threshold), the number of the simulation, and the order of removal
  # Simulation 1: all the species are present in the regional pool.
  # simul = 30: juste one species remains in the regional pool.
  # Order is an element of c("increasing", "decreasing", "random_1" ,  "random_2")
  res<-read.table(paste0("data/raw/output-cmd2_",site,"_",order,".txt/forceps.",site,".site_1_complete.txt"),header=F) 
  colnames(res) <- colnames(read.table(here::here("data","colnames_res.txt"),header=T))
  colnames(res)<-colnames_res
  temp_plot<-filter(temporal_plot(res),date==max(unique(res$date)))
  temp_plot$simul <- rep(1,dim(temp_plot)[1])
  data<-temporal_plot_threshold(temp_plot)
  for(i in c(2:30)){
    res<-try(read.table(paste0("data/raw/output-cmd2_",site,"_",order,".txt/forceps.",site,".site_",i,"_complete.txt"),header=F),silent=T) 
    if (class(res) != "try-error"){ # sometimes, the files are empty, and it returns an error message
      colnames(res)<-colnames_res
      temp_plot<-filter(temporal_plot(res),date==max(unique(res$date)))
      temp_plot$simul <- rep(i,dim(temp_plot)[1])
      data<-rbind(data,temporal_plot_threshold(temp_plot))
    }
  }
  data$order <- rep(order,dim(data)[1])
  write.table(data,paste0("data/processed/table_price_threshold_",site,"_",order,".txt"))
}


pairwise_analysis_price <- function(table_price){
  # Returns a table of CAFE values (SRE.L, etc.) for pairs of simulations. 
  # Simulation 1: all the species are present in the regional pool. Simulation 30: juste one species remains in the regional pool.
  # I follow the steps from the priceTools package examples. To install it from GitHub, see the appendix of Bannar Martin's 2017 paper.
  data<- table_price
  grouped.data<-group_by(data,simul,order)
  res1 <- grouped.data %>%
    pairwise.price(species="species",func="biomass")
  pp1<-group.columns(res1,gps=c('simul'),drop=F)
  dat<-filter(pp1,simul.x < simul.y) # so that I don't have the symmetric comparisons of simul 1 with 2, and of simul 2 with 1, etc.
  
  dat
}


CAFE_graphics <- function(dat1){
  # Returns the three plots from Bannar-Martin's paper.
  # dat1 is just one line from the pairwise.price analysis, i.e. a comparison between two communities.
  # e.g. dat1 <- read.table("data/processed/Pairwise_CAFE_values_Bern.txt",header=T)
  
  # Richness-Composition vectors
  s1<-leap.zig(dat1,type='bef',standardize=FALSE,
               # xlim=c(10,40),ylim=c(50,300),
               error.bars=F,
               vectors=T,raw.points = F,legend=TRUE) +
    scale_y_continuous("Ecosystem function \n(biomass (g))",
                       breaks=c(50,100,150,200,250,300))+
    annotate("text", x = mean(dat1$x.rich), y = mean(dat1$x.func), 
             label = "*",size=8)+
    annotate("segment", x = mean(dat1$y.rich)-1, xend = mean(dat1$y.rich)+1, 
             y = mean(dat1$y.func), yend = mean(dat1$y.func),colour = "black") +
    ggtitle("Richness-Composition")
  # s1
  
  # Community assembly vectors
  s2<-leap.zig(dat1,type='cafe',standardize=FALSE,
               # xlim=c(10,40),ylim=c(50,300),
               error.bars=F,
               vectors=T,raw.points = F,legend=TRUE)+
    scale_y_continuous("",breaks=c(50,100,150,200,250,300))+
    annotate("text", x = mean(dat1$x.rich), y = mean(dat1$x.func), 
             label = "*",size=8)+
    annotate("segment", x = mean(dat1$y.rich)-1, xend = mean(dat1$y.rich)+1, 
             y = mean(dat1$y.func), yend = mean(dat1$y.func),colour = "black")+
    ggtitle("Community Assembly")
  # s2
  
  # 5-part Price vectors
  s3<-leap.zig(dat1,type='price',standardize=FALSE,
               # xlim=c(10,40),ylim=c(50,300),
               error.bars=F,
               vectors=T,raw.points = F,legend=TRUE)+
    scale_y_continuous("",breaks=c(50,100,150,200,250,300))+
    annotate("text", x = mean(dat1$x.rich), y = mean(dat1$x.func), 
             label = "*",size=8)+
    annotate("segment", x = mean(dat1$y.rich)-1, xend = mean(dat1$y.rich)+1, 
             y = mean(dat1$y.func), yend = mean(dat1$y.func),colour = "black")+
    ggtitle("5-part Price")
  # s3
  
  grid.arrange(s1,s2,s3,nrow=2)
}