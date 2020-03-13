change_temp <- function(site,temp_increase){
  # A function that increases or decreases the mean temperature at a given site, in the climate file.
  # example: change_temp("Bern",-1)
  
  # site: a character string, the name of a site. For instance, "Bever".
  # temp_increase is a number, either positive or negative, which will be added to the mean temperature  in the climate file at every month of every year. 
  
  # Writes a climate file. If the temperature was increased by 1 degree, the file's name is of the form "forceps.climate.Bern1.txt"
  # If it was decreased by 1 degree, it is "forceps.climate.Bern-1.txt"
  
  climate <- read.table(paste0("data/raw/Files necessary to run the simulations/Climate files/forceps.climate.",site,".txt"),header=F)
  colnames(climate) <- c("Year","Month","Tmean","Tmin","Tmax","SumPrec","NbRainyDays")
  climate$Tmean <-  climate$Tmean + temp_increase
  write.table(climate,paste0("data/raw/Files necessary to run the simulations/Chain rule simulations/Climate files/forceps.climate.",site,"+",temp_increase,".txt")
              ,row.names = F,col.names = F,sep="\t")
}