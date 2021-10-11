# Infer key speciess by applying Maire's framework can be applied to identify key species for 
# productivity, selection, and complementarity.
# Adapted from Maire, E., Villéger, S., Graham, N., Hoey, A., Cinner, J.,    
# Ferse, S., Aliaume, C, Booth, D., Feary, D., Kulbicki, M., Sandin, S., Vigliola, L. & Mouillot, D. (2018).    
# Community-wide scan identifies fish species associated to coral reef services globally. Proceedings  B. 

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