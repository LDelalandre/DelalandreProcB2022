source("final/0. Packages.R")
source("R/Analysis_data.R")
source("R/Common variables.R")

# I) Generate data ####
df_maire <- function(LH_all_per_sp,sit,simul_to_keep){
  df_LH <- LH_all_per_sp %>% 
    filter(site == sit)
  
  dat <- data.frame(matrix(nrow=0,ncol=37))
  colnames(dat) <- c("site","order","simul",
                     as.character(read.table("data/raw/distinctiveness of the species.txt",header = T)$SName),
                     "productivity","DeltaY","Selection","Complementarity")
  j <- 1 # number of the assemblage = position (line) in the data.frame dat
  
  for (order in ORDER[3:32]){
    a <- order
    df_LH_order <- subset(df_LH,order==a)
    
    i = simul_to_keep # keep only simul 16, with 15 species
    sub <- subset(df_LH_order, simul==i) 
    # filter(persists_mixt == T) # NB SPECIES THAT ARE PRESENT IN THE SIMUL ARE ABOVE BIOMASS THRESHOLD
    # NB2: I did not filter on species persisting in mono and mixt to performe Loreau-Hector analysis.
    
    # dat[j,]$assemblage <- paste0(sit,"_",order,"_",i) # paste order and simul (an assemblage is a simul)
    dat[j,]$site <- sit
    dat[j,]$order <- order
    dat[j,]$simul <- i
    
    dat[j,]$productivity <- sum(sub$YOi)
    dat[j,]$DeltaY <- unique(sub$DeltaY)
    dat[j,]$Selection <- unique(sub$Selection)
    dat[j,]$Complementarity <- unique(sub$Complementarity)
    
    # all the species that we saw at the end of the simul are considered persistent (bad)
    # spnames <- as.character(sub$SName) 
    
    # only species above biomass threshold are considered persistent (better, I think)
    spnames <- sub %>% 
      pull(SName) %>% 
      as.character
    
    
    
    dat[j,which(!(colnames(dat) %in% c(spnames,"site","order","simul","productivity","DeltaY","Selection","Complementarity") ))] <- 0
    dat[j,which(colnames(dat) %in% sub$SName)] <- 1
    
    j <- j+1 # current position in the data.frame dat.
  }
  dat
} # end of df_maire

LH_all_per_sp <- read.csv("data/processed/Loreau-Hector coefficients_per species.csv") %>% 
  mutate(site = factor(site,levels=SITE)) %>% 
  filter(persists_mixt==T) #/!\ problem with Loreau-Hector computation 
# (but the problem might not be significant : test how much filtering on persist_mono and mixt
# changes LH coefficients).

for (site in SITE){
  for (nb_sp in c(10:20) ){
    simul_to_keep <- 30 - nb_sp +1
    dat <- df_maire(LH_all_per_sp,site,simul_to_keep)
    write.table(dat,paste0("data/processed/maire/occurrence_ppty",site,"_",30-simul_to_keep+1,"species.txt"),sep="\t",row.names=F)
  }
}


# II) Analyses ####
infer_key_species <- function(occurrence_ppty){
  # Candidate species ####
  max <- dim(occurrence_ppty)[2]-5 # /!\ SI JE REGARDE A RICHESSE CONSTANTE IL FAUT METTRE -4 !!!!!!
  Moccur <- occurrence_ppty[,c(4: max )]
  
  limit <- round(0.1*nrow(Moccur),0) # 10%
  occurrence <- apply(Moccur,2,sum)
  length(which(occurrence>=limit)) #we retained 381 fish species as candidates 
  required_threshold <- which(occurrence>=limit)
  
  # Keeping only candidates species and their presence/absence 
  presence_candidates <- Moccur[,required_threshold]
  candidates_species <- colnames(presence_candidates)
  
  
  # STEP 1: Computing initial fish biomass model (M0) ####
  data_M0 <- occurrence_ppty %>% 
    rename(property = !!property) 
    # mutate(richness = as.factor(richness)) # richness as a factor or a continuous predictor
  M0 <- lm(property~ richness,data=data_M0)
  
  # # Box-Cox
  # bc <- MASS::boxcox(M0,lambda=seq(-2,2,length=200))
  # lambda <- bc$x[which.max(bc$y)]
  # z <- (data_M0$property^lambda-1)/lambda
  # M0 <- lm(z ~ data_M0$richness)
  
  # # Test the model
  # par(mfrow=c(2,2)) ; plot(M0)
  # shapiro.test(residuals(M0)) # not normal
  # lmtest::bptest(M0) # Not homoscedastic
  # lmtest::dwtest(M0) # Independent
  
  
  AIC_M0 <- AIC(M0)
  
  #STEP 2: TESTING THE EFFECT OF EACH FISH SPECIES INDIVIDUALLY ####
  # For each candidate species, an alternative model (M1) is obtained by adding the presence/absence of this fish species to M0
  # AIC of M1 and the effect (positive or negative) of the candidate species are saved.
  
  #Results for all species will be stored in the Test_candidates matrix
  Test_candidates <- matrix(NA,length(candidates_species),3)
  colnames(Test_candidates) <- c("occurrence","AIC_M1","coeff_sp")
  rownames(Test_candidates) <- candidates_species
  
  
  for (k in 1:length(candidates_species)) {
    # Defining presence/absence of the species k
    candidate <- as.factor(presence_candidates[,k])
    Test_candidates[k,"occurrence"] <- length(which(candidate==1))
    
    # ppty_data <- occurrence_ppty %>% 
    #   pull(!!property)
    # data_with_candidate <-as.data.frame(cbind(ppty_data,as.character(candidate))) 
    
    #Adding presence/absence of species k to M0 to obtain M1
    cand_sp <- candidates_species[k]
    data_test <- occurrence_ppty %>% 
      select(!!property,!!cand_sp)
    colnames(data_test) <- c("property","candidate")
    
    M1 <- lm(property ~ candidate,data=data_test)
    
    Test_candidates[k,"AIC_M1"] <- AIC(M1)
    Test_candidates[k,"coeff_sp"] <- unique(coef(M1)[2])
    
  } # end of k
  
  Test_candidates <- as.data.frame(Test_candidates)
  
  
  # STEP 3: Determining species that improve the prediction ####
  # as ΔAIC (AIC_M0 - AIC_M1) > 4 with a positive effect
  Delta_AIC <- AIC_M0 - Test_candidates$AIC_M1
  better_AIC <- ifelse(Delta_AIC>=4,1,0)
  positive <- ifelse(Test_candidates$coeff_sp>0,1,0)
  
  key_species <- which(better_AIC==1 & positive ==1)
  
  #Extracting matrix that combines "key species" and their performances in the model
  Summary_key_species <- cbind(Test_candidates,Delta_AIC)[key_species,]
  nrow(Summary_key_species) 
  #26 fish species are significantly and positively related to fish biomass
  
  Summary_key_species
} # end of function infer_key_species

PROPERTIES <- c("productivity","DeltaY","Selection","Complementarity")
property <- "Complementarity"
site <- SITE[5]
nb_sp <- 15

# account for richness
SUMMARY_KEY_SPECIES <- NULL
for (site in SITE){
  for (property in PROPERTIES){
    OCCURRENCE_PPTY <- NULL
    for (nb_sp in c(10:20) ){
      occurrence_ppty <- read.table(paste0("data/processed/maire/occurrence_ppty",site,"_",nb_sp,"species.txt"),header=T)
      OCCURRENCE_PPTY <- rbind(OCCURRENCE_PPTY,occurrence_ppty)
    }
    OP <- OCCURRENCE_PPTY %>% 
      mutate(richness = 30 - simul + 1)
    
    Summary_key_species <- infer_key_species(OP) %>% 
      rownames_to_column("SName")
    if (dim(Summary_key_species)[1] > 0){
      Summary_key_species$site <- site
      Summary_key_species$property <- property
    }
    SUMMARY_KEY_SPECIES <- rbind(SUMMARY_KEY_SPECIES,Summary_key_species)
  }
}




# Previous (richness by richness) method
SUMMARY_KEY_SPECIES <- NULL
for (site in SITE){
  for (property in PROPERTIES){
    occurrence_ppty <- read.table(paste0("data/processed/maire/occurrence_ppty",site,"_",nb_sp,"species.txt"),header=T)
    Summary_key_species <- infer_key_species(occurrence_ppty) %>% 
      rownames_to_column("SName")
    if (dim(Summary_key_species)[1] > 0){
      Summary_key_species$site <- site
      Summary_key_species$property <- property
    }
    SUMMARY_KEY_SPECIES <- rbind(SUMMARY_KEY_SPECIES,Summary_key_species)
  }
}

# view(SUMMARY_KEY_SPECIES %>% filter(site == "Sion"))
SUMMARY_KEY_SPECIES %>% 
  filter(coeff_sp < 0)

#_______________________________________________________________________________
# Plot nb of key species ####
species <- read.table("data/raw/distinctiveness of the species.txt",header=T) %>% 
  pull(SName) 
distinct_sp <- species[1:10]
common_sp <- species[11:30]

SUMM2 <- SUMMARY_KEY_SPECIES %>%
  mutate(property = if_else(property=="productivity","Productivity",property)) %>% 
  mutate(status = if_else(SName %in% distinct_sp,"Distinct","Common")) %>% 
  arrange(factor(property, levels = c("Productivity","DeltaY","Selection","Complementarity"))) %>% 
  mutate(property = factor(property,levels=c("Productivity","DeltaY","Selection","Complementarity"))) %>% 
  arrange(factor(site, levels = SITE)) %>% 
  mutate(site = factor(site,levels=SITE)) %>% 
  filter(!(property == "DeltaY")) %>% 
  group_by(site,property,status) %>% 
  mutate(count = n()) %>% 
  mutate(status = if_else(status=="Common","Ordinary",status))

cols <- c("Distinct" = "#F8766D", "Ordinary" = "#00BFC4")
ggplot(SUMM2,aes(x=property,fill=status))+
  geom_histogram(stat="count") +
  facet_wrap(~site)+
  theme(axis.text.x = element_text(angle = 66)) +
  scale_color_manual(values = cols,
                     aesthetics = "fill") +
  theme(axis.text.x = element_text(hjust=1)) +
  ggsave(paste0("figures/2021_09_lm_maire/nb_key_sp.png"))

  
SUMM3 <- SUMM2 %>% 
  group_by(site,property,status) %>% 
  mutate(count = n()) %>% 
  summarize(coeff = mean(coeff_sp),count=mean(count))


ggplot(SUMM3,aes(x=property,y=coeff,fill=status))+
  # stat_summary(geom = "bar", fun = mean, position = "dodge") +
  # stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge",width = 0.5)+
  scale_color_manual(values = cols, aesthetics = "fill") +
  geom_histogram(stat="identity",position = "dodge") +
  theme(axis.text.x = element_text(angle = 66,hjust=1)) +
  ylim(c(0,2))+
  facet_wrap(~site) +
  geom_text(aes(label=count), vjust=0,position = position_dodge(width = 1)) +
  labs(fill = "Category of species") +
  xlab("Property") +
  ylab("Mean effect of key species") +
  ggsave(paste0("figures/2021_09_lm_maire/mean_effect_key_sp.png"),height=7,width=7)

  




# table
ppty <- "DeltaY"

SUMM2 %>% 
  select(site,property,SName,coeff_sp) %>% 
  filter(property == ppty) %>% 
  filter(SName %in% distinct_sp) %>% 
  spread(SName,coeff_sp) %>%  
  kable(escape = F) %>%
  kable_styling()

SUMM2 %>% 
  select(site,property,SName,coeff_sp) %>% 
  filter(property == ppty) %>% 
  filter(SName %in% common_sp) %>% 
  spread(SName,coeff_sp) %>%  
  kable(escape = F) %>%
  kable_styling()



#________________________________________
# PLots ####
# Represent the amount of each property in each site, 
# measured on the 30 random simulations containing all 15 species
PROPERTIES <- c("productivity","DeltaY","Selection","Complementarity")
PLOTS <- NULL
i <- 0
for (site in SITE){
    i <- i+1
    occurrence_ppty <- read.table(paste0("data/processed/maire/occurrence_ppty",site,"_15species.txt"),header=T)
    toplot <- occurrence_ppty %>% 
      select(site,order,productivity,DeltaY,Selection,Complementarity) %>% 
      gather(key = property,value = value,-c(site,order)) %>% 
      arrange(factor(property, levels = PROPERTIES)) %>% 
      mutate(property = factor(property,levels=PROPERTIES))
    
    plot <- ggplot(toplot,aes(x=property,y=value))+
      ggtitle(site) +
      geom_boxplot() +
      ylim(c(-4,6))
    PLOTS[[i]] <- plot
    # ggsave(plot = plot,filename = paste0("figures/2021_09_lm_maire/15species_",site,".png"))
}

PLOT <- PLOTS

complete_plot <- 
  grid.arrange( PLOT[[1]],PLOT[[2]], PLOT[[3]],
                                PLOT[[4]], PLOT[[5]],PLOT[[6]],
                                PLOT[[7]],PLOT[[8]],PLOT[[9]],
                                PLOT[[10]], PLOT[[11]],
                                ncol = 3, nrow = 4, 
                                layout_matrix = rbind( c(0,1,2),c(3,4,5),c(6,7,8),
                                                       c(9,10,11)))

ggsave(plot = complete_plot,
       filename = paste0("figures/2021_09_lm_maire/properties_15species.png"),
       width = 10, height = 10)

