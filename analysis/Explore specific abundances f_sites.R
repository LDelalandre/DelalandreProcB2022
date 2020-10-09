TOTAL <- read.table("data/processed/specific_biom_prod_complete.txt",header=T)

GD <- TOTAL %>% 
  subset(site=="GrandeDixence"&order=="increasing"&simul==1) %>% 
  select(site,species,mixture.t.ha.,prod_mixture)

Bern <- TOTAL %>% 
  subset(site=="Bern"&order=="increasing"&simul==1) %>% 
  select(site,species,mixture.t.ha.,prod_mixture)

Schwerin <- TOTAL %>% 
  subset(site=="Schwerin"&order=="increasing"&simul==1) %>% 
  select(site,species,mixture.t.ha.,prod_mixture)

plot(Bern$prod_mixture,GD$prod_mixture)
Bern[which(Bern$prod_mixture>mean(Bern$prod_mixture)),]

Bern[which(GD$prod_mixture>mean(GD$prod_mixture)),]
?order
?sort
Bern[order(Bern$prod_mixture,decreasing=T),]
GD[order(GD$prod_mixture,decreasing=T),]

Bern[order(Bern$mixture.t.ha.,decreasing=T),]
GD[order(GD$mixture.t.ha.,decreasing=T),]
Schwerin[order(Schwerin$mixture.t.ha.,decreasing=T),]


tr <- filter(traits,SName %in%c("PAbi","LDec","PCem"))
