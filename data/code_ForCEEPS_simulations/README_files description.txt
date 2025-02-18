* In a setup file, what changes between sites are two lines :

siteFileName = forceps.Bever.site
speciesFileName = forceps.species_PG
defaultClimateFileName = forceps.climate.Bever.txt

siteFileName = forceps.Bern.site
speciesFileName = forceps.species_PG
defaultClimateFileName = forceps.climate.Bern.txt


* In the Cmd files, I only call forceeps.setup, 
/!\ so in this file I have to change the two previous lines each time I change the site !