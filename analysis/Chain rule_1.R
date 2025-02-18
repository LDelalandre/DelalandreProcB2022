source("R/Chain rule_functions.R")
change_temp("Bern",-1)

colnames_res<-colnames(read.table("data/colnames_res.txt", header=T))
SITE <- c("Bern","Bever","Cottbus","Huttwil")
site <- SITE[1]
temp_increase <- "1"

decreased <- read.table("data/raw/output-cmd2_Bern_4sp_Temp+-1.txt/forceps.Bern.site_1_complete.txt")
baseline <- read.table("data/raw/output-cmd2_Bern_4sp_Temp.txt/forceps.Bern.site_1_complete.txt")
increased <- read.table("data/raw/output-cmd2_Bern_4sp_Temp+1.txt/forceps.Bern.site_1_complete.txt")
                              
colnames(decreased) <- colnames_res
decreased_final <- subset(decreased,date==max(decreased$date))

#  on a zi = ni * wi. ni=abondance de l'espèce i ; wi = biomasse moyenne d'un individu de l'espèce i. zi = biomasse de toute l'espèce.
# on veut calculer le delta, différence entre deux simulations

SP <- c()
NI <- c()
WI <- c()

for (sp in unique(decreased_final$speciesShortName)){
  sub <- subset(decreased_final,speciesShortName==sp)
  SP <- c(SP,sp)
  WI <- c(WI,mean(sub$biomass.kg.) )
  NI <- c(NI,length(sub$biomass.kg.))
}
data.frame(sp=SP,ni=NI,wi=WI)
