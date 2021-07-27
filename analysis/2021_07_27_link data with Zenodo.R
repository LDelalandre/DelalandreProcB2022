# Download data from zenodo
library("zen4R")
# download_zenodo("zenodo doi", path = "dossier_voulu")
download_zenodo("10.5281/zenodo.808257", path = "test")


# Unzip data
download(url, dest="dataset.zip", mode="wb") 
unzip ("dataset.zip", exdir = "./")

unzip("TEST zip/Test2 to zip.zip",overwrite=F,exdir="TEST zip") # Works with files that are not heavy. Otherwise, very very very long.
