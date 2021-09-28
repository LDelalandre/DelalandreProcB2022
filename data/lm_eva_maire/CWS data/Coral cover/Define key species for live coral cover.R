################################################################################################################
# R script to define species that are potential key contributor to a given ecosystem service.                  #
#                                                                                                              #
# Step 1: The goal is to model a given ecosystem service (Y) according to Environmental (E) and                #
# Socio-Economic conditions (SE) and species richness (R); to check its relevance according to its             #
# explanatory power and to save its Akaike Information Criterion (AIC_M0) as a reference for the next step     #
#                                                                                                              #
# Step2: The goal is to identify species key for the studied ecosystem service (Y) adding each candidate       #
# species (presence-absence) as an additional explanatory variable to M0 to compute model M1 and its           #
# associated AIC (AIC_M1). Finally, a species is declared as a key potential contributor to the ecosystem      #
# service if ΔAIC (AIC_M0-AIC_M1) > 4 and if its partial effect is positive (positive coefficient in the model)#
#                                                                                                              #
# Here we apply this framework to define the "key fish species" for live coral cover                           # 
# Code by Eva Maire (eva.maire@umontpellier.fr). Last update on 2018-06-26                                     #
# This code provides results of analysis used in Maire, E., Villéger, S., Graham, N., Hoey, A., Cinner, J.,    #
# Ferse, S., Aliaume, C, Booth, D., Feary, D., Kulbicki, M., Sandin, S., Vigliola, L. & Mouillot, D. (2018).   # 
# Community-wide scan identifies fish species associated to coral reef services globally. Proceedings  B.      #
#                                                                                                              #
################################################################################################################

#loading required libraries
require(MuMIn)
require(lme4)

# setting working directory
my_path<-"" #  <= "Coral cover" folder
setwd(my_path)

###################################################################################################################################
# LOADING AND PREPARING DATA
# (i) Loading complete dataset of environnemental and socioeconomic drivers, coral cover and species richness
# All the quantitative covariables are standardised. 

load("CoralCover_data.Rdata")
nrow(CoralCover_data) #741 reefs are used in the analysis (have coral cover data)

# (ii) Loading the matrix which describes presence/absence of fish species for each reef: Reefs x Species
load("presence_fish_species_CC.Rdata")

# We excluded fish species present on less than 1% of the reefs (i.e. 7 reefs for coral cover dataset)
limit <- round(0.01*nrow(CoralCover_data),0)
fish_sp <- colnames(presence_fish_species_CC)

occurrence <- apply(presence_fish_species_CC,2,sum)

length(which(occurrence>=limit)) #we retained 335 fish species as candidates 

required_threshold <- which(occurrence>=limit)

#Keeping only candidates species and their presence/absence 
candidates_species <- fish_sp[required_threshold]
presence_candidates <- presence_fish_species_CC[,required_threshold]

#Verifying that the 2 datasets have the same rows
length(is.element(CoralCover_data$Reef_ID,rownames(presence_candidates))) # 741 OK

###################################################################################################################################
#DEFINING KEY SPECIES

# STEP 1: Computing initial live coral cover (M0)

M0 <- lmer(logCoralCover~Management+Gravity_markets+Gravity_near_pop+Species_richness+
             Local_population_growth+HDI+Voice_accountability+Tourism+
             Population_size+Reef_fish_landings+Oceanic_productivity+Depth+Habitat+Census_method+
             (1|Nation)+(1|Site)+(1|Regional_location),data=CoralCover_data)

#Checking for the robustness of M0
AIC(M0)
r.squaredGLMM(M0)[2] #Radj = 0.63

#Storing AIC of M0 that will be the threshold to define which species improve the prediction of coral cover
AIC_M0 <- AIC(M0)

#STEP 2: TESTING THE EFFECT OF EACH FISH SPECIES INDIVIDUALLY
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
  
  data_with_candidate <-as.data.frame(cbind(CoralCover_data,candidate)) 
  
  #Adding presence/absence of species k to M0 to obtain M1
  M1 <- lmer(logCoralCover~candidate+Management+Gravity_markets+Gravity_near_pop+Species_richness+
               Local_population_growth+HDI+Voice_accountability+Tourism+
               Population_size+Reef_fish_landings+Oceanic_productivity+Depth+Habitat+Census_method+
               (1|Nation)+(1|Site)+(1|Regional_location),data=data_with_candidate)  
  
  Test_candidates[k,"AIC_M1"] <- AIC(M1)
  Test_candidates[k,"coeff_sp"] <- unique(coef(M1)$Nation$candidate)
  
} # end of k

Test_candidates <- as.data.frame(Test_candidates)

#Determining species that improve the prediction of coral cover as ΔAIC (AIC_M0 - AIC_M1) > 4 with a positive effect
Delta_AIC <- AIC_M0-Test_candidates$AIC_M1
better_AIC <- ifelse(Delta_AIC>=4,1,0)
positive <- ifelse(Test_candidates$coeff_sp>0,1,0)

key_species <- which(better_AIC==1 & positive ==1)

#Extracting matrix that combines "key species" and their performances in the model
Summary_key_species <- cbind(Test_candidates,Delta_AIC)[key_species,] 
nrow(Summary_key_species) 
#28 fish species are significantly and positively related to live coral cover

Summary_key_species

#Functional traits of key species
load("Functional_traits.Rdata")

Functional_traits_key_species <- Functional_traits[which(Functional_traits$Species_code%in%rownames(Summary_key_species)==T),]
Functional_traits_key_species #see Supplementary Material for more details on functional traits

###################################################################################################################################
# END OF SCRIPT
###################################################################################################################################
###################################################################################################################################



