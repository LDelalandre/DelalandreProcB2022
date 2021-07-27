source("final/0. Packages.R")
source("R/Before simulations.R")
source("R/Common variables.R")
source("R/comp_fct_dist.R")

traits<-read.table("data/raw/Traits of the species_complete.txt",header=T)
selected_traits<-select_traits(traits) # data.frame with the traits of the species

#________________________________________________________________________________
# Table 2: Correlation prod-distinctiveness in monoculture ####
MIXTURES <- read.table("data/processed/LH_productivity_specific_every condition.txt",header=T) %>% 
  filter(persists_mixt==T)

# correlation tests
get_pval2 <- function(MIXT){
  cortest <- cor.test(MIXT$Di,MIXT$monoculture_t_ha,method="spearman",exact = FALSE)
  cortest$p.value
}

get_cor2 <- function(MIXT){
  cortest <- cor.test(MIXT$Di,MIXT$monoculture_t_ha,method="spearman",exact = FALSE)
  cortest$estimate
}


COR <- NULL
PVAL <- NULL
for (sit in SITE){
  M <- MIXTURES %>% 
    filter(site == sit) 
  COR <- c(COR,as.numeric(get_cor2(M)))
  PVAL <- c(PVAL,get_pval2(M))
  
}

correlations <- data.frame(site=SITE,cor=round(COR,digits=2),pval=PVAL) #pval=round(PVAL,digits=4))


correlations2 <- correlations %>% 
  mutate(p.value = map_dbl(pval,round,4)) %>% 
  arrange(desc(cor)) %>% 
  mutate(p.value = ifelse(p.value < 0.0001,
                          "<0.0001",
                          p.value),)



tablecor <- correlations2 %>%
  transmute(
    Site = site,
    Correlation = ifelse(p.value =="<0.0001",
                         cell_spec(cor, bold = T),
                         cell_spec(cor, bold=F)),
    p.value = ifelse(p.value =="<0.0001",
                     cell_spec(p.value, bold = T),
                     cell_spec(p.value, bold=F)),
    
  ) %>%
  select(Site, everything()) %>%
  kable( escape = F) %>%
  kable_styling("hover", full_width = F)
cat(tablecor, file = "figures_tables/Table_2_Corelation_prod_Di.doc")


#________________________________________________________________________________
# Table Sx: Correlations between traits ####
corr <- cor(selected_traits)
plot_cor <- ggcorrplot(corr,
                       hc.order = TRUE,
                       type = "lower",
                       lab = TRUE,digits=1)
ggsave ("figures_tables/Fig. SX_correlation between traits.png",plot=plot_cor,dpi="print",width=17,units="cm")

#________________________________________________________________________________
# Table Sx: Trait values and distinctiveness ####
snsn <- read.table("data/raw/correspondence_SName_Id.txt",header=T)
Species_name <- snsn$Name
Short_name <- snsn$SName

Distinctiveness <- distinct_tot$Di
table1 <- cbind(Species_name,Short_name,selected_traits,Di=round(Distinctiveness,digits=3))

tabletrait <- table1 %>%
  kable( escape = F) %>%
  kable_styling("hover", full_width = F)
cat(tabletrait, file = "figures_tables/Table x_Species trait values.doc")

