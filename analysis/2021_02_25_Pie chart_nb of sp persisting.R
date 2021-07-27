library("tidyverse")


data <- read.table("data/processed/LH_productivity_specific_every condition.txt",header=T) %>% tibble()
data2 <- data %>% 
  filter(simul==1 & order=="decreasing") %>% # with 30 species
  filter(persists_mixt==T) %>% #▲ which species persist in mixture
  mutate(group_of_sp = map_chr(species,group_of_sp))


nb_persist <- data2 %>% group_by(site,group_of_sp) %>% 
  summarize(n=n()) %>% 
  spread(group_of_sp,n)

# Environmental plot ####
sites <- read.table("data/Site description.txt",header=T)

get_coordinates <- function(sit,metric){
  sites2 <- sites %>% 
    filter(Site==sit) %>% 
    select(Site,Temp_moy,Annual_ppt)
  if (metric == "T"){
    sites2 %>% pull(Temp_moy)
  } else if (metric =="P"){
    sites2 %>% pull(Annual_ppt)
  } else {
    "metric should be either T or P"
  }
}

nb_persist2 <- nb_persist %>% 
  mutate(Temp_moy = get_coordinates(site,"T")) %>% 
  mutate(Annual_ppt = get_coordinates(site,"P")) %>% 
  replace(is.na(.), 0) %>% 
  mutate(ntot=distinct10+common20)

environment2 <-  ggplot(sites, aes(x=Temp_moy,y=Annual_ppt,label=Site))+
  # geom_point()+
  xlab("Average temperature (°C)")+
  ylab("Annual precipitations (mm)") +
  xlim(0,12) +
  ylim(380,1500)+
  egg::theme_article(base_size=18)

environment2 + 
  geom_scatterpie(data = nb_persist2,mapping=aes(x=Temp_moy, y=Annual_ppt, group=site, r=ntot),
                                color=NA, alpha=.8) +
  geom_scatterpie_legend(nb_persist2$ntot, x=2.5, y=580)
  
  environment2 + 
    geom_scatterpie(data = nb_persist2,aes(x=Temp_moy, y=Annual_ppt, group=site, r=ntot),
                    cols=c("common20","distinct10"),color=NA,alpha=0.8) 
  
  nb_persist3 <- nb_persist2
  nb_persist3$common20 <- nb_persist2$common20+0.5
  nb_persist3$distinct10 <- nb_persist2$distinct10+0.5
  
environment2 + 
  geom_scatterpie(data = nb_persist3,aes(x=Temp_moy, y=Annual_ppt, group=site, r=ntot),
                  cols=c("common20","distinct10"),color=NA,alpha=0.8) 

#________________________________
library(scatterpie)

set.seed(123)
long <- rnorm(50, sd=100)
lat <- rnorm(50, sd=50)
d <- data.frame(long=long, lat=lat)
d <- with(d, d[abs(long) < 150 & abs(lat) < 70,])
n <- nrow(d)
d$region <- factor(1:n)
d$A <- abs(rnorm(n, sd=1))
d$B <- abs(rnorm(n, sd=2))
d$C <- abs(rnorm(n, sd=3))
d$D <- abs(rnorm(n, sd=4))
d$radius <- 6 * abs(rnorm(n))
head(d)

world <- map_data('world')
p <- ggplot(world, aes(long, lat)) +
  geom_map(map=world, aes(map_id=region), fill=NA, color="black") +
  coord_quickmap()
p + geom_scatterpie(aes(x=long, y=lat, group=region, r=radius),
                    data=d, cols=LETTERS[1:4], color=NA, alpha=.8) +
  geom_scatterpie_legend(d$radius, x=-160, y=-55)

