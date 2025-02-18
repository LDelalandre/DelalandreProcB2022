library("tidyverse")

source("R/Common variables.R")
source("R/Before simulations.R")
source("R/comp_fct_dist.R")


traits<-read.table("data/Traits of the species_complete.txt",header=T)

selected_traits1<-choice_traits_1(traits) # data.frame with the traits of the species
selected_traits2<-selected_traits1

decreasing_stepwise <- NULL
increasing_stepwise <- NULL
for(i in c(1:29)){
  # decreasing
  most_dist_sp <- 
    comp_fct_dist(selected_traits1) %>% 
    filter(Di==max(Di)) %>% 
    pull(SName)
  
  selected_traits1 <- 
    selected_traits1 %>% 
    filter( !(rownames(.) == most_dist_sp) )
  
  decreasing_stepwise <- c(decreasing_stepwise,most_dist_sp)
  
  # increasing
  least_dist_sp <- 
    comp_fct_dist(selected_traits2) %>% 
    filter(Di==min(Di)) %>% 
    pull(SName)
  
  selected_traits2 <- 
    selected_traits2 %>% 
    filter( !(rownames(.) ==least_dist_sp) )
  
  increasing_stepwise <- c(increasing_stepwise,least_dist_sp)
}

decreasing_stepwise
increasing_stepwise

stepwise_orders <- data.frame(dec_step = decreasing_stepwise, inc_step = increasing_stepwise)

# compare to the order we have
old_orders <- read.table("data/removal order and sp names.txt",header=T) %>% 
  select(species_incr,species_decr) %>% 
  rename(dec_old = species_decr,inc_old = species_incr)


orders <- cbind(stepwise_orders,old_orders)

species <- 
  traits %>% 
  arrange(SName) %>% 
  pull(SName)

# Have a data frame with each species and its ranking, in terms of distinctiveness, with both methods and in both orders of removal
ORDERS <- 
  data.frame(species) %>% 
  mutate(dec_step = map_int(species,function(sp) which(orders$dec_step==sp) )) %>%
  mutate(dec_old = map_int(species,function(sp) which(orders$dec_old==sp) )) %>%
  mutate(inc_step = map_int(species,function(sp) which(orders$inc_step==sp) )) %>% 
  mutate(inc_old = map_int(species,function(sp) which(orders$inc_old==sp) ))


corINC <- with(ORDERS,cor.test(dec_step,dec_old,method="spearman"))
corDEC <- with(ORDERS,cor.test(inc_step,inc_old,method="spearman"))

corINC$estimate

INC <- ggplot(ORDERS,aes(x=inc_old,y=inc_step,label=species))+
  geom_label(colour="#00BFC4") +
  annotate("text", x=25, y=5, label= paste("Spearman's rho =", round(corINC$estimate,digits=3), 
                                           "\n p.value =",round(corINC$p.value,digits=6) ))+
  ggtitle("Removing common species first") +
  xlab("Distinctiveness computed on the 30 species") +
  ylab("Distinctiveness computed each time a species is removed")

DEC <- ggplot(ORDERS)+
  geom_label(aes(x=dec_old,y=dec_step,label=species),colour="#F8766D") +
  annotate("text", x=25, y=5, label= paste("Spearman's rho =", round(corDEC$estimate,digits=3), 
                                           "\n p.value =",round(corDEC$p.value,digits=6) )) +
  ggtitle("Removing distinct species first") +
  xlab("Distinctiveness computed on the 30 species") +
  ylab("Distinctiveness computed each time a species is removed")

# combine the two plots
combined <- cowplot::plot_grid(DEC,INC, labels=c("A", "B"), ncol = 2, nrow = 1)
# If I want to add a title to it all : https://github.com/wilkelab/cowplot/blob/master/vignettes/plot_grid.Rmd
#  Section "joint plot title


ggsave(combined,filename="paper_2/Stepwise Di computation.png",
       width = 35,
       height = 18,
       units = "cm",
       dpi = 300)
