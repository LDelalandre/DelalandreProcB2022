## Specify the packages of interest ####
packages = c("tidyverse", 
             "FactoMineR", "factoextra", # PCA computation
             "funrar", # Functional distinctiveness computation
             "cowplot", "gridExtra", "ggsignif", "egg", "ggpubr","ggcorrplot", # Packages used for plots
             "kableExtra", # Represent tables
             "zen4R", # download data from Zenodo repository
             "functClust" # functional clustering
)

## Load, or install and load, packages
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
    }
    library(x, character.only = TRUE)
  }
)


# Complete folder structure ####

# Create folder to export figures if it does not exist #
if(file.exists(file.path("figures_tables"))==F){
  dir.create(file.path("figures_tables"))  
}