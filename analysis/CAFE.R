# Done using the code from "BaylesRdInvasionExampleRCode_101.R"
library("dplyr")
library("priceTools")
library("ggplot2")
library(gridExtra)

SITE <- c("Bern","Bever","Cottbus","Huttwil")
site <- SITE[1]
ORDER <- c("increasing", "decreasing", "random_1" ,  "random_2")
order <- ORDER[1]

source("R/Analysis_data.R")
for (order in ORDER){
  table_threshold_price(site,order)
}


# Price analysis ####
# I follow the steps from the priceTools package examples. To install it from GitHub, see the appendix of Bannar Martin's 2017 paper.
table_price <- read.table(paste0("data/processed/table_price_threshold_",ORDER[1],".txt"),header=T)
table_price <- rbind(table_price, read.table(paste0("data/processed/table_price_threshold_",ORDER[2],".txt"),header=T) )


pairwise_analysis_price <- function(table_price){
  # gives a table of CAFE values (SRE.L, etc.) for pairs of simulations. 
  # Simulation 1: all the species are present in the regional pool. Simulation 30: juste one species remains in the regional pool.
  data<- table_price
  grouped.data<-group_by(data,simul,order)
  res1 <- grouped.data %>%
    pairwise.price(species="species",func="biomass")
  pp1<-group.columns(res1,gps=c('simul'),drop=F)
  dat<-filter(pp1,simul.x < simul.y) # so that I don't have the symmetric comparisons of simul 1 with 2, and of simul 2 with 1, etc.
}

dat <- pairwise_analysis_price(table_price)


# Look at a specific comparison
simulx<-1
simuly<-10
dat1<-subset(dat,simul.x==simulx & simul.y==simuly)

# Look at the different effects ####

baseline_1 <- filter(dat,simul.x==1)
hist(baseline_1$SL)
hist(baseline_1$SG)
hist(baseline_1$CDE) # CDE is in average positive. It means that the remaining species increase thier biomass. This is not surprising.

plot(baseline_1$SL~baseline_1$simul.y)
plot(baseline_1$SRE.L~baseline_1$simul.y)
plot(baseline_1$y.rich~baseline_1$simul.y,type='l') # Globally, there is a decrease in richness (but with fluctuations) along the simulations. Consequently:
plot(baseline_1$y.func~baseline_1$simul.y,type='l')

plot(baseline_1$SRE.L~baseline_1$simul.y) # The richness effect increases in absolute value (because more and more sp are removed).
plot(baseline_1$SIE.L~baseline_1$simul.y)
# plot(baseline_1$SRE.G~baseline_1$simul.y)
# plot(baseline_1$SIE.G~baseline_1$simul.y)




# Graphics ####
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




