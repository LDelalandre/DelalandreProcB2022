# Script to download and extract raw data #

# NB: raw data weighs ~ 120Go. They will be available as zipped files on Zenodo soon.
# Raw data curation (aggregate biomass and productivity from the individual level to species level) lasts ~2 to 3 hours on my computer.
# Intermediary data is available in the "processed" folder.
# If you want to avoid data curation step, jump directly to script 3. to 6. in the "final" folder. 

source("final/0. Packages.R")
source("R/Common variables.R")

ask = function( prompt ) {
  cat( paste( prompt, ':' ) )
  readLines( n=1 )
}

n = as.integer( ask( 'Please enter sample size n' ) )

ask = function( prompt ) {
  cat( paste0( prompt, ':' ) )
  readLines( n=1 )
}

choice = as.character( ask( 'Do you want to download raw data from Zenodo? (yes/no) \n
N.B. Raw data is heavy (~120 Go) and can be long to download and unzip. ' ) )

if (choice == "yes"){
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
} else if (choice=="no") {
  cat("An alternative option is to go on with the analysis from data aggregated at the species level. \n
To do so, please go to scripts 3. to 6. in the folder \"final\".")
}
