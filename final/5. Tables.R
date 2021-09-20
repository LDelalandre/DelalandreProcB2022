source("final/0. Packages.R")
source("R/Before simulations.R")
source("R/Common variables.R")
source("R/comp_fct_dist.R")

traits<-read.table("data/raw/Traits of the species_complete.txt",header=T)
selected_traits<-select_traits(traits) # data.frame with the traits of the species

#________________________________________________________________________________
# Table 2: Correlation prod-distinctiveness in monoculture ####
PROD_mono <- read.table("data/processed/productivity_monoculture_ALL sites.txt",header=T) # productivity
BIOM_mono <- read.table("data/processed/biomass_mono_ALL sites.txt",header=T) %>% 
  rename(monoculture = monoculture.t.ha.) # biomass

# correlation tests
get_pval2 <- function(MIXT){
  cortest <- cor.test(MIXT$Di,MIXT$monoculture,method="spearman",exact = FALSE)
  cortest$p.value
}

get_cor2 <- function(MIXT){
  cortest <- cor.test(MIXT$Di,MIXT$monoculture,method="spearman",exact = FALSE)
  cortest$estimate
}

# Biomass 
COR <- NULL
PVAL <- NULL
for (sit in SITE){
  M <- BIOM_mono %>% 
    filter(site == sit) 
  COR <- c(COR, as.numeric(get_cor2(M)) )
  PVAL <- c(PVAL,get_pval2(M))
  
}
correlations <- data.frame(site=SITE,cor=round(COR,digits=2),pval=PVAL) #pval=round(PVAL,digits=4))
correlations_biom <- correlations %>% 
  mutate(p.value = map_dbl(pval,round,2)) %>% 
  select(-pval) %>% 
  mutate(p.value = ifelse(p.value < 0.01,
                          "<0.01",
                          p.value),)

# Productivity 
COR <- NULL
PVAL <- NULL
for (sit in SITE){
  M <- PROD_mono %>% 
    filter(site == sit) 
  COR <- c(COR, as.numeric(get_cor2(M)) )
  PVAL <- c(PVAL,get_pval2(M))
  
}
correlations <- data.frame(site=SITE,cor=round(COR,digits=2),pval=PVAL) #pval=round(PVAL,digits=4))
correlations_prod <- correlations %>% 
  mutate(p.value = map_dbl(pval,round,2)) %>% 
  select(-pval) %>% 
  mutate(p.value = ifelse(p.value < 0.01,
                          "<0.01",
                          p.value),)

correlations2 <- merge(correlations_biom,correlations_prod,by="site") %>% 
  mutate(site = factor(site,levels=SITE)) %>%
  arrange(site)

tablecor <- correlations2 %>%
  transmute(
    Site = site,
    Correlation = ifelse(p.value.x =="<0.01"|p.value.x =="0.01",
                         cell_spec(cor.x, bold = T),
                         cell_spec(cor.x, bold=F)),
    p.value = ifelse(p.value.x =="<0.01"|p.value.x =="0.01",
                     cell_spec(p.value.x, bold = T),
                     cell_spec(p.value.x, bold=F)),
    Correlation2 = ifelse(p.value.y =="<0.01"|p.value.y =="0.01",
                          cell_spec(cor.y, bold = T),
                          cell_spec(cor.y, bold=F)),
    p.value2 = ifelse(p.value.y =="<0.01"|p.value.y =="0.01",
                      cell_spec(p.value.y, bold = T),
                      cell_spec(p.value.y, bold=F))
  ) %>%
  select(Site, everything()) %>%
  kable( escape = F,
         col.names = c("Site","Correlation","p.value","Correlation","p.value")) %>%
  kable_styling("hover", full_width = F) %>% 
  add_header_above(c(" " = 1,"Biomass" = 2,"Productivity" = 2))

cat(tablecor, file = "figures_tables/Table_2_Correlation_prod_Di.doc")


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
