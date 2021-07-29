### Functional distinctiveness and ecosystem productivity

This repository contains the code used for the article "Functionally distinct tree species support long-term productivity in extreme environments".

You can clone the repository with git by pasting in your terminal:

	git clone https://github.com/LDelalandre/Project-1_Distinct-sp_BEF.git
    
or 
just download the repository:
[LDelalandre/Project-1_Distinct-sp_BEF](https://github.com/LDelalandre/Project-1_Distinct-sp_BEF/archive/master.zip).

If you have [Rstudio](https://www.rstudio.com/) installed on your computer, you can then open `Project-1_Distinct-sp_BEF.Rproj` with Rstudio.

You can contact me at <leo.delalandre@cefe.cnrs.fr>


### Structure of the project repository

#### R

Contains functions that are used in the other scripts.

#### final

Scripts used for the article.

NB: If you don not want to download raw data (~120 Go), you can start from data aggregated at the species level by avoiding the script `2. Process data.R` and jump directly to the following scripts.

#### analysis

Exploratory scripts, that were not retained for the publication.

#### data 

##### > code_ForCEEPS_simulations 
Code to perform simulations using the ForCEEPS model, hosted on 
[Capsis](http://capsis.cirad.fr/capsis/help_en/forceeps) (contact them for any informations about the platform, how to install it on your machine, etc.).

##### > raw 
Contains raw ForCEEPS simulation outputs. These are not hosted on the present GitHub repository. 
Script to download them from [Zenodo](https://zenodo.org/record/5145755) is to be found in `final/2. Process data.R`.

##### > processed 

Processed data. Available from the present GitHub repository. NB: processed datasets are generated from: `final/2. Process data.R` 
only if raw data have been downloaded and de-zipped. Otherwise, you can start the analyses starting from processed data.




### Analyses performed

The present scripts were used to perform forest succession simulations with the ForCEEPS forest gap model, extracts biomass and productivity outputs and analyses them. Simulations were performed along a broad environmental gradient, either in monocultures or in mixtures whereby species richness gradients enabled us to test the consequences of the presence or absence of functionally distinct species on ecosystem productivity. See the article for more details. 
