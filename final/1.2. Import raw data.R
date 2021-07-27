# Script to download and extract raw data #

# NB: raw data weight ~ 120Go. They will be available as zipped files on Zenodo soon.
# Raw data curation (aggregate biomass and productivity from the individual level to species level) lasts ~2 to 3 hours on my computer.
# Intermediary data is available in the "processed" folder.
# If you want to avoid data curation step, jump directly to script 3. to 6. in the "final" folder. 

source("final/0. Packages.R")
source("R/Common variables.R")

# Check directory for ForCEEPS output is available, and create it if not ####
if(file.exists(file.path("data/raw/Output_ForCEEPS"))==F){
  dir.create(file.path("data/raw/Output_ForCEEPS"))  
}
   
# Download compressed raw data from Zenodo ####
download_zenodo("10.5281/zenodo.808257", path = "data/raw/Output_ForCEEPS") # change Zenodo doi when my data is uploaded

# Extract files ####
# Can take a few hours.
for (site in SITE){
  unzip(paste0("data/raw/Output_ForCEEPS",site,".zip"),
        overwrite=F,
        exdir= paste0("data/raw/Output_ForCEEPS",site,".zip")) # Takes some ours to run.
}
