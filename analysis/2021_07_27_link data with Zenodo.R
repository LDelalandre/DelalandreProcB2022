# ZENODO ####
# Download data from zenodo
library("zen4R")
# download_zenodo("zenodo doi", path = "dossier_voulu")
path <- "test"
if(file.exists(file.path(path))==F){
  dir.create(file.path(path))  
}
download_zenodo("10.5281/zenodo.808257", path = path)


# Unzip data
download(url, dest="dataset.zip", mode="wb") 
unzip ("dataset.zip", exdir = "./")

unzip("TEST zip/Test2 to zip.zip",overwrite=F,exdir="TEST zip") # Works with files that are not heavy. Otherwise, very very very long.


# DRYAD #### bof bof...
library("rdryad")
list_dryad <- dryad_dataset(doi = "10.5061/dryad.f385721n")
