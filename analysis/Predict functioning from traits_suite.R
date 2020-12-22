library("lme4")
mono <- read.table("data/processed/data productivity~traits_monoculture.txt",header=T)
mixt <- read.table("data/processed/data productivity and biomass~traits_mixture.txt",header=T)

sit <- SITE[5]
test <- 
  mixt %>% # I work either on MONOCULTURES or MIXTURES here (choose mixt or mono)
  group_by(site) %>% 
  mutate(sum_biom=sum(biomass)) %>% 
  filter(biomass>0.001*sum_biom)

# centrer et réduire les valeurs de traits
test <- as.data.frame(test)
for(i in c(1:12)){
  test[,i+2] <- ( test[,i+2] - mean(test[,i+2]) )/sd(test[,i+2])
}

# Predict biomass ####
mod=lmer(log(biomass)~S +   HMax + AMax+  G  +  DDMin  +   WiTN  + DrTol    +        
           Brow + Ly +WiTX + NTol +(1|site) ,data=test) # passage de la biomasse au log améliore l'homoscédasticité et la normalité
mod=lmer(log(biomass)~S +   HMax + AMax+  G  +  DDMin  +   WiTN  + DrTol    +        
           Brow + Ly +WiTX + NTol +(1|site) ,data=test) # passage de la biomasse au log améliore l'homoscédasticité et la normalité


# 1. Vérification des hypothèses
summary(mod)
plot(mod) # homoscédasticité et linéarité
qqnorm(resid(mod)) # normalité des résidus
hist(resid(mod))

qqnorm(scale(resid(mod)))
abline(a=0,b=1, col=2) #permet de tracer droite

# 2. Test du modèle
mod0 <- lmer(log(biomass)~1 +(1|site) ,data=test)
anova(mod,mod0)

# 3. Test des effets fixes (les traits)
car::Anova(mod)

# 4. Test de l'effet aléatoire (site) " je ne sais pas comment faire
mod1 <- nlme::gls(log(biomass)~S +   HMax + AMax+  G  +  DDMin  +   WiTN  + DrTol    +        
                   Brow + Ly +WiTX + NTol  ,data=test)
anova(mod,mod1)




# Predict productivity ####
mod=lmer(log(prod)~S +   HMax + AMax+  G  +  DDMin  +   WiTN  + DrTol    +        
           Brow + Ly +WiTX + NTol +(1|site) ,data=test) # passage de la biomasse au log améliore l'homoscédasticité et la normalité
# je suis passé au log pour Ly et DrTol à cause du message : Warning message:
# Some predictor variables are on very different scales: consider rescaling 
summary(mod)

# 1. Vérification des hypothèses
plot(mod) # homoscédasticité et linéarité
qqnorm(resid(mod)) # normalité des résidus

qqnorm(scale(resid(mod)))
abline(a=0,b=1, col=2) #permet de tracer droite

# 2. Test du modèle
mod0 <- lmer(log(prod)~1 +(1|site) ,data=test)
anova(mod,mod0)

# 3. Test des effets fixes (les traits)
car::Anova(mod)


# interactions between traits that predict prod and sites ####
mod=lm(prod^0.5 ~ S*site +   HMax*site + Brow*site +WiTX*site + NTol*site ,data=test)
par(mfrow=c(2,2)) ; plot(mod)
summary(mod)
car::Anova(mod)

# Predict prod on three subset of sites ####

test2 <- filter(test,site %in% c("Bever","Davos","GrandeDixence")) # cold sites
# test2 <- filter(test,site %in% c("Adelboden","Huttwil","Bern")) # warm-wet sites
# test2 <- filter(test,!(site %in% c("Bever","Davos","GrandeDixence","Adelboden","Huttwil","Bern")))#warm-dry sites

mod=lmer(prod^0.5~S +   HMax + AMax+  G  +  DDMin  +   WiTN  + DrTol    +        
           Brow + Ly +WiTX + NTol +(1|site) ,data=test2) # passage de la biomasse au log améliore l'homoscédasticité et la normalité


# 1. Vérification des hypothèses
# summary(mod)
plot(mod) # homoscédasticité et linéarité
# qqnorm(resid(mod)) # normalité des résidus

qqnorm(scale(resid(mod))) ; abline(a=0,b=1, col=2) #permet de tracer droite

# 2. Test du modèle
# mod0 <- lmer(log(prod)~1 +(1|site) ,data=test2) ;anova(mod,mod0)
summary(mod)
# 3. Test des effets fixes (les traits)
car::Anova(mod)




step(mod, scope = ~ 1 ,
     direction="both", criterion = "AIC")

submod <- lm(formula = log(prod) ~ AMax + DDMin + DrTol + Brow + Ly + WiTX, 
             data = test2)
car::Anova(submod)

# Correlate biomass and productivity ####
for (sit in SITE){
  monosite <- filter(mono,site==sit)
  plot(monosite$prod~monosite$biomass,main=paste0("Monoculture ",sit))
}

for (sit in SITE){
  mixtsite <- filter(mixt,site==sit)
  plot(mixtsite$prod~mixtsite$biomass,main=paste0("Mixture ",sit))
  mod <- lm((mixtsite$prod)^0.5~mixtsite$biomass)
  par(mfrow=c(2,2)) ; plot(mod)
  mod0 <- lm(mixtsite$prod~1)
  an <- anova(mod,mod0)
  print(an$`Pr(>F)`)
  print(summary(mod)$r.squared)
}
# la biomasse et la productivité sont bien corrélés en mélange (modèles significatifs, et r² élevé)
# par contre, comme toujours, je galère avec le respect des hypothèses des modèles


# Predict persistence in mixture ####
mixt$persist <- mixt$biomass
mixt$persist[which(mixt$biomass!=0)] <- 1
mixt$persist <- as.factor(mixt$persist)

mixt2 <- filter(mixt,site %in% c("Bever","Davos","GrandeDixence")) # cold sites
mixt2 <- filter(mixt,site %in% c("Adelboden","Huttwil","Bern")) # warm-wet sites
mixt2 <- filter(mixt,!(site %in% c("Bever","Davos","GrandeDixence","Adelboden","Huttwil","Bern")))#warm-dry sites
# mixt2 <- filter(mixt,site %in% c("Sion","Cottbus")) 

modbin=glm(persist~ S +   HMax + AMax+  G  +  DDMin  +   WiTN  + DrTol    +        
             Brow + Ly +WiTX + NTol ,family=binomial,data=mixt2)
modbin3=glm(persist~  DDMin  ,family=binomial,data=mixt2)

modbin0=glm(persist~1,family=binomial,data=mixt2)
anova(modbin0,modbin,test="Chisq")
a <- car::Anova(modbin)


# write.table(a,"paper/anova_persist_cold.txt",sep="\t")
# write.table(a,"paper/anova_persist_warm-wet.txt",sep="\t")
# write.table(a,"paper/anova_persist_warm-dry.txt",sep="\t")
