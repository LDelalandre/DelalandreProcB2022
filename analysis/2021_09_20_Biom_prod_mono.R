library(dplyr)

# Correlate biomass and productivity in monoculture, in every site.
# (and log of biomass and productivity).

prod <- read.table("data/processed/productivity_monoculture_ALL sites.txt",header=T) # productivity
biom <- read.table("data/processed/biomass_mono_ALL sites.txt",header=T) %>% 
  rename(monoculture = monoculture.t.ha.) # biomass

sit <- "Cottbus"

BIOM_PROD <- NULL
for (sit in SITE){
  prodsite <- prod %>% 
    filter(site==sit) %>% 
    rename(prod = monoculture)  %>% 
    select(SName,prod)
  biomsite <-  biom %>% 
    filter(site==sit) %>% 
    rename(biom = monoculture)  %>% 
    select(SName,biom,Di,Id,site,persists)
  
  biom_prod <- merge(prodsite,biomsite,by="SName") %>% 
    mutate(site = sit)
  
  # ggplot(biom_prod,aes(x = biom, y = prod,label = SName))+
  #   geom_point()+
  #   ggrepel::geom_label_repel() 
  
  BIOM_PROD <- rbind(BIOM_PROD,biom_prod)
}
BIOM_PROD2 <- BIOM_PROD %>% 
  mutate(site = factor(site,levels=SITE)) %>% 
  mutate(log_biom = log(biom),
         log_prod = log(prod))


ggplot(BIOM_PROD2,aes(x = biom, y = prod,label = SName))+
  geom_point()+
  facet_wrap(~site)


ggplot(BIOM_PROD2,aes(x = log_biom, y = log_prod,label = SName))+
  geom_point()+
  facet_wrap(~site)



