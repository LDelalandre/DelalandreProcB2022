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


# Write the tables ####
# NB: Takes a few minuts!
for (site in SITE){
  # Write tables of final biomasses to be used for the Price analysis
  for (order in ORDER){
    table_threshold_price(site,order)
  }
  
  # Compute all pairwise CAFE analysis inside a site and export them in a .txt
  table_price <- read.table(paste0("data/processed/table_price_threshold_",site,"_",ORDER[1],".txt"),header=T)
  for (order in ORDER){
    table_price <- rbind(table_price, read.table(paste0("data/processed/table_price_threshold_",site,"_",order,".txt"),header=T) )
  }
  dat <- pairwise_analysis_price(table_price)
  write.table(dat,paste0("data/processed/Pairwise_CAFE_values_",site,".txt"))
}



# Price analysis ####
dat <- read.table("data/processed/Pairwise_CAFE_values.txt")


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
plot(baseline_1$SRE.L~baseline_1$simul.y,main="Richness effect of species lost")
plot(baseline_1$SIE.L~baseline_1$simul.y,main="Identity effect of species lost")

plot(baseline_1$y.rich~baseline_1$simul.y,type='l') # Globally, there is a decrease in richness (but with fluctuations) along the simulations. Consequently:
plot(baseline_1$y.func~baseline_1$simul.y,type='l')

plot(baseline_1$SRE.L~baseline_1$simul.y) # The richness effect increases in absolute value (because more and more sp are removed).
plot(baseline_1$SIE.L~baseline_1$simul.y)
# plot(baseline_1$SRE.G~baseline_1$simul.y)
# plot(baseline_1$SIE.G~baseline_1$simul.y)




