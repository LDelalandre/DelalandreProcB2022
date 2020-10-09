##############################################################################################################################################################
##############################################################################################################################################################
#### Libraries 
library(effects)
library(car)
library(multcomp)
library(lme4) # mixed models
library(nlme) # mixed models
library(doBy)
library(Hmisc)
library(lsmeans)
library(smatr)
library(reshape)

digeco.data<-read.table("/Users/iprieto/Dropbox/2.CEFE/PRAISE/Base_donnees_Praise/digeco.data.csv", dec = ".", sep = ";", na="NA", header=TRUE)
names(digeco.data)

### Mean and SE functions ###
mean1 <- function(x) mean(x,na.rm=TRUE)
se <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
n<-function(x) length(na.omit(x))

## For the cumulated biomasss in both drought and water treatments (communities with 10 genotypes) 
# mscum2.6 and mscum(sp) are calculated withougt the second measurement (date 2) since we had no data for individual species ####
CE<-digeco.data[digeco.data$nbgeno.sp == 10 ,c("com","com2","sp.com","traitH2O","nbsp","nbgeno.sp","msdac6","msrga6","msfet6","msluz6","mstre6","mstot6","mscum6","mscum2.6","mscumdac6","mscumluz6","mscumfet6","mscumrga6","mscumtre6")]
CE[,c("msdac6","msrga6","msfet6","msluz6","mstre6","mstot6","mscum6","mscum2.6","mscumdac6","mscumluz6","mscumfet6","mscumrga6","mscumtre6")]=CE[,c("msdac6","msrga6","msfet6","msluz6","mstre6","mstot6","mscum6","mscum2.6","mscumdac6","mscumluz6","mscumfet6","mscumrga6","mscumtre6")]/576*10000

## Calculate M (mean yield of species i in monoculture) ######
Mdac<-mean(CE[CE$sp.com == "DAC",]$mscum2.6)
Mluz<-mean(CE[CE$sp.com == "LUZ" ,]$mscum2.6)
Mrga<-mean(CE[CE$sp.com == "RGA" ,]$mscum2.6)
Mfet<-mean(CE[CE$sp.com == "FET" ,]$mscum2.6)
Mtb<-mean(CE[CE$sp.com == "TB" ,]$mscum2.6)
M<-(Mdac+Mluz+Mrga+Mfet+Mtb)/5

## Calculate Ye #######
CE$RYej<-CE$mscum2.6/5
RYej<-0.2

Ye.dac<-RYej*Mdac
Ye.luz<-RYej*Mluz
Ye.fet<-RYej*Mfet
Ye.rga<-RYej*Mrga
Ye.tb<-RYej*Mtb
Ye<- Ye.dac+Ye.luz+Ye.fet+Ye.rga+Ye.tb

## Calculate DeltaRY ######
## Calculate RYoj (mean observed relative yield of species i in the mixture) #######
CE$RYoj.dac<-CE$mscumdac6/Mdac
CE$RYoj.luz<-CE$mscumluz6/Mluz
CE$RYoj.rga<-CE$mscumrga6/Mrga
CE$RYoj.fet<-CE$mscumfet6/Mfet
CE$RYoj.tb<-CE$mscumtre6/Mtb

## Calculate DeltaRY ######
CE$DeltaRY.dac<-CE$RYoj.dac-RYej
CE$DeltaRY.luz<-CE$RYoj.luz-RYej
CE$DeltaRY.rga<-CE$RYoj.rga-RYej
CE$DeltaRY.fet<-CE$RYoj.fet-RYej
CE$DeltaRY.tb<-CE$RYoj.tb-RYej

CE$DeltaRY <- (CE$DeltaRY.dac+CE$DeltaRY.luz+CE$DeltaRY.rga+CE$DeltaRY.fet+CE$DeltaRY.tb)/5

## Calculate Net effect (Yo - Ye) ###
Net<-CE[CE$sp.com == "toutes",]
Net$Neteffect<-Net$mscum2.6-Ye

Net.mean<-summaryBy(Neteffect~sp.com, Net, FUN=c(mean1,se,n)) # Mean net effect of all communities 

#### Complementary effect ############
Net$CE<-5*Net$DeltaRY*M

#### Selection effect (Net effect = CE + SE so it follows that SE = Net effect - CE ############
Net$SE<-Net$Neteffect-Net$CE

Net$logSE_CE<-log(abs(Net$SE)/abs(Net$CE))

Net.mean<-summaryBy(Neteffect~sp.com, Net, FUN=c(mean1,se,n)) # Mean net effect of all communities 
CE.mean<-summaryBy(CE~sp.com, Net, FUN=c(mean1,se,n)) # Mean complementarity effect of all communities 
SE.mean<-summaryBy(SE~sp.com, Net, FUN=c(mean1,se,n)) # Mean selection effect of all communities 
SE_CE.mean<-summaryBy(logSE_CE~sp.com, Net, FUN=c(mean1,se,n)) # Mean selection effect of all communities 

plot(Net$SE~Net$CE)

regCE<-lm(Net$Neteffect~Net$CE)
summary(regCE)

regSE<-lm(Net$Neteffect~Net$SE)
summary(regSE)

plot(Net$Neteffect~Net$CE, xlab="Complementarity effect, CE (g)", ylab="Net diversity effect (g)", xlim=c(18,95), ylim=c(-15,52), pch=19, las=1)
abline(regCE, lty=2)

plot(Net$Neteffect~Net$SE, xlab="Selection effect, SE (g)", xlim=c(-82,-20), ylim=c(-15,52), ylab="Net diversity effect (g)",  pch=19, las=1)
abline(regSE, lty=2)

plot(logSE_CE.mean1~1, xaxt="n", xlab="", data=SE_CE.mean, pch=5, ylim=c(-0.4,0.1))
arrows(1,SE_CE.mean$logSE_CE.mean1, x1, SE_CE.mean$logSE_CE.mean1 + SE_CE.mean$logSE_CE.se, length=0)
arrows(1,SE_CE.mean$logSE_CE.mean1, x1, SE_CE.mean$logSE_CE.mean1 - SE_CE.mean$logSE_CE.se, length=0)

## Stats for net effect ######

Net.test<-t.test(Net$Neteffect, mu=0, na.rm=TRUE, var.equal = TRUE, alternative = c("greater"))
Net.test.wat<-t.test(Net[Net$traitH2O == 2,]$Neteffect, mu=0, na.rm=TRUE, var.equal = TRUE, alternative = c("greater"))
Net.test.drought<-t.test(Net[Net$traitH2O == 1,]$Neteffect, mu=0, na.rm=TRUE, var.equal = TRUE, alternative = c("greater"))

CE.test<-t.test(Net$CE, mu=0, na.rm=TRUE, var.equal = TRUE, alternative = c("greater"))
CE.test.wat<-t.test(Net[Net$traitH2O == 2,]$CE, mu=0, na.rm=TRUE, var.equal = TRUE, alternative = c("greater"))
CE.test.drought<-t.test(Net[Net$traitH2O == 1,]$CE, mu=0, na.rm=TRUE, var.equal = TRUE, alternative = c("greater"))

SE.test<-t.test(Net$SE, mu=0, na.rm=TRUE, var.equal = TRUE, alternative = c("less"))
SE.test.wat<-t.test(Net[Net$traitH2O == 2,]$SE, mu=0, na.rm=TRUE, var.equal = TRUE, alternative = c("less"))
SE.test.drought<-t.test(Net[Net$traitH2O == 1,]$SE, mu=0, na.rm=TRUE, var.equal = TRUE, alternative = c("less"))

SE_CE.test<-t.test(Net$logSE_CE, mu=0, na.rm=TRUE, var.equal = TRUE, alternative = c("less"))
SE_CE.test.wat<-t.test(Net[Net$traitH2O == 2,]$logSE_CE, mu=0, na.rm=TRUE, var.equal = TRUE, alternative = c("less"))
SE_CE.test.drought<-t.test(Net[Net$traitH2O == 1,]$logSE_CE, mu=0, na.rm=TRUE, var.equal = TRUE, alternative = c("less"))

glm.Net<-lm(Neteffect~traitH2O, data=Net)
Anova(glm.Net)

glm.CE<-lm(CE~traitH2O, data=Net)
Anova(glm.CE)

glm.SE<-lm(SE~traitH2O, data=Net)
Anova(glm.SE)




### Figure histogram for CE, SE and net effects
FigNet<-melt(Net[,c(1,2,3,4,5,6,32,33,34)], id=c("com","com2","sp.com","traitH2O","nbsp","nbgeno.sp"))
plot(value~variable,data=FigNet)

# ## For the cumulated biomasss in both drought and water treatments ponderated by biomass ####
# CE<-digeco.data[digeco.data$nbgeno.sp == 10 ,c("com","com2","sp.com","traitH2O","nbsp","nbgeno.sp","msdac6","msrga6","msfet6","msluz6","mstre6","mstot6","mscum6")]
# CE$ab.dac<-CE$msdac6
# 
# CE<-melt(CE, id=c("com","com2","sp.com","traitH2O","nbsp","nbgeno.sp"))
# colnames(CE)=c("com","com2","sp.com","traitH2O","nbsp","nbgeno.sp","sp","ms")
# CE<-na.omit(CE)
# 
# CE$cwms<-CE$ms/5
# 
# cwms_mean<-summaryBy(cwms~sp.com, data=CE[CE$sp.com !="toutes" & CE$sp =="mscum6",], FUN=c(mean1,se,n))
# # colnames(cwms_mean)=c("sp.com", "ms","ms.se","ms.n")
# # colnames(cwms_mean.mixed)=c("sp.com", "ms","ms.se","ms.n")
# x<-sum(as.vector(cwms_mean$cwms.mean1))
# 
# Net<-CE[CE$sp.com == "toutes" & CE$sp == "mscum6",]
# Net$Neteffect<-Net$ms-x
# 
# Net.mean<-summaryBy(Neteffect~sp.com, Net, FUN=c(mean1,se,n))
# 
# Net.test<-t.test(Net$Neteffect, mu=0, na.rm=TRUE, var.equal = TRUE, alternative = c("greater"))
# 

# ## For the cumulated biomasss in both drought and water treatments by genotype ####
# CE<-digeco.data[digeco.data$sp.com == "toutes" ,c("com","com2","sp.com","traitH2O","nbsp","nbgeno.sp","msdac6","msrga6","msfet6","msluz6","mstre6","mstot6","mscum6")]
# 
# CE<-melt(CE, id=c("com","com2","sp.com","traitH2O","nbsp","nbgeno.sp"))
# colnames(CE)=c("com","com2","sp.com","traitH2O","nbsp","nbgeno.sp","sp","ms")
# # CE<-na.omit(CE)
# 
# CE$cwms<-CE$ms/5
# 
# cwms_mean<-summaryBy(cwms~sp, data=CE[CE$nbgeno.sp==1 & CE$sp !="mscum6",], FUN=c(mean1,se,n))
# colnames(cwms_mean)=c("sp", "ms","ms.se","ms.n")
# 
# x<-sum(as.vector(cwms_mean$ms))
# 
# Net<-CE[CE$sp == "mstot6" & CE$nbgeno.sp != 1,]
# Net$Neteffect<-Net$ms-x
# 
# Net.mean<-summaryBy(Neteffect~sp.com, Net, FUN=c(mean1,se,n))
# 
# Net.test<-t.test(Net$Neteffect, mu=0, na.rm=TRUE, var.equal = TRUE, alternative = c("greater"))

# CE<-digeco.data[digeco.data$traitH2O == 2 & digeco.data$nbgeno.sp == 10 ,c("com","com2","sp.com","traitH2O","nbsp","nbgeno.sp","msdac6","msrga6","msfet6","msluz6","mstre6","mstot6")]
# 
# CE<-melt(CE, id=c("com","com2","sp.com","traitH2O","nbsp","nbgeno.sp"))
# colnames(CE)=c("com","com2","sp.com","traitH2O","nbsp","nbgeno.sp","sp","ms")
# CE<-na.omit(CE)
# 
# CE$cwms<-CE$ms/5
# 
# cwms_mean<-summaryBy(cwms~sp.com, data=CE[CE$sp.com !="toutes" & CE$sp =="mstot6",], FUN=c(mean1,se,n))
# # colnames(cwms_mean)=c("sp.com", "ms","ms.se","ms.n")
# # colnames(cwms_mean.mixed)=c("sp.com", "ms","ms.se","ms.n")
# sum(as.vector(cwms_mean$cwms.mean1))
# 
# rbind(cwms_mean,cwms_mean.mixed)
# 

# ## For both drought and water treatments####
# CE<-digeco.data[digeco.data$nbgeno.sp == 10 ,c("com","com2","sp.com","traitH2O","nbsp","nbgeno.sp","msdac6","msrga6","msfet6","msluz6","mstre6","mstot6")]
# 
# CE<-melt(CE, id=c("com","com2","sp.com","traitH2O","nbsp","nbgeno.sp"))
# colnames(CE)=c("com","com2","sp.com","traitH2O","nbsp","nbgeno.sp","sp","ms")
# CE<-na.omit(CE)
# 
# CE$cwms<-CE$ms/5
# 
# cwms_mean<-summaryBy(cwms~sp.com, data=CE[CE$sp.com !="toutes" & CE$sp =="mstot6",], FUN=c(mean1,se,n))
# colnames(cwms_mean)=c("sp.com", "ms","ms.se","ms.n")
# y<-sum(as.vector(cwms_mean$ms))
# 
# Net_tot<-CE[CE$sp.com == "toutes" & CE$sp == "mstot6",]
# Net_tot$Neteffect<-Net_tot$ms-y
# 
# Net.test.tot<-t.test(Net_tot$Neteffect, mu=0, na.rm=TRUE, var.equal = TRUE)

plot.net <- summaryBy(value~variable, data=FigNet, FUN=c(mean1,se,n))

x1<-c(1,2,3)
plot(value.mean1~x1, type="n", xaxt="n", yaxt="n", ylab=expression("Biodiversity effects (g"~m^-2*")"), xlab="", ylim=c(-950,1300), xlim=c(0.5,3.5),  las=1, data=plot.net)
points(value.mean1~x1, pch=19, data=plot.net)
axis(2, labels=TRUE, tick = TRUE, lwd = 0, lwd.ticks = 1, tck=0.025, las=2)
axis(1,at=x1, labels=c("Net","CE","SE"), tick = TRUE, lwd = 0, lwd.ticks = 1, tck=0.025)
minor.tick(ny=5, nx=0, tick.ratio=-0.5)
arrows(x1,plot.net$value.mean1,x1, plot.net$value.mean1 + plot.net$value.se, length=0)
arrows(x1,plot.net$value.mean1,x1, plot.net$value.mean1 - plot.net$value.se, length=0)
text(x1,1250,labels=c("**","***","***"), cex=1.1)
abline(0,0, lty=1)

plot.net2 <- summaryBy(value~variable*traitH2O, data=FigNet, FUN=c(mean1,se,n))
x1<-c(0.9,1.1,1.9,2.1,2.9,3.1)
plot(value.mean1~x1, type="n", xaxt="n", yaxt="n", ylab=expression("Biodiversity effects (g"~m^-2*")"), xlab="", ylim=c(-950,1300), xlim=c(0.5,3.5),  las=1, data=plot.net2)
points(value.mean1~c(0.9,1.9,2.9), pch=1, data=plot.net2[plot.net2$traitH2O==1,])
points(value.mean1~c(1.1,2.1,3.1), pch=19, data=plot.net2[plot.net2$traitH2O==2,])
legend("topright", c("drought","irrigated"), pch = c(1,19), bty="n")
axis(2, labels=TRUE, tick = TRUE, lwd = 0, lwd.ticks = 1, tck=0.025, las=2)
axis(1,at=c(1,2,3), labels=c("Net","CE","SE"), tick = TRUE, lwd = 0, lwd.ticks = 1, tck=0.025)
minor.tick(ny=5, nx=0, tick.ratio=-0.5)
arrows(x1,plot.net2$value.mean1,x1, plot.net2$value.mean1 + plot.net2$value.se, length=0)
arrows(x1,plot.net2$value.mean1,x1, plot.net2$value.mean1 - plot.net2$value.se, length=0)
text(x1,1250,labels=c("**","***","***"), cex=1.1)
abline(0,0, lty=1)

