source("R/Common variables.R")

# Biomass - species removal #####
biomass_specific_aggreg_mono <- function(site) {
  # 1. Takes the output of the simulations (complete files).
  # 2. Selects biomass above threshold
  # 3. Selects a subset of ten years and averages the biomass of each species on these years
  # 4. Returns data frame (e.g. specific_biomass_final_Bern.txt) whose columns are:
  # "species" "abundance"	"mixture(t/ha)"	"mixture_relative"	"site"	"order"	"simul"
  RES <- read.table(paste0("data/processed/Aggregated monocultures_method 2/artificial monocultures_",sit,".txt"),header=T)
  
  BIOMASSES <- as.data.frame(matrix(nrow=0,ncol=4,dimnames = list(NULL,c("species", "mixture(t/ha)", "mixture_relative", "simul"))))
  for (order in ORDER){
    orde <- order
    biomass_inc<-c()
    sd_biom_inc<-c()
    for(number in c(1:30)){
      res <- subset(RES,order==orde & simul==number)
      # res<-try(read.table(paste0("data/raw/output-cmd2_",site,"_",order,".txt/forceps.",site,".site_",number,"_complete.txt")),silent=T) 
      if (class(res) != "try-error"){# sometimes, the files are empty, and it returns an error message
        # colnames(res) <- colnames_res
        int <- temporal_plot(res)
        temp_plot <- temporal_plot_threshold(int) # NB: if error here, don't forget to charge the package "dplyr"
        
        dates <- as.numeric(max(temp_plot$date))- c(900,800,700,600,500,400,300,200,100,0) # years on which we average the biomass
        years_to_keep <- subset(temp_plot,date %in%dates) # we compute a mean value of biomass on those years
        
        if (dim(years_to_keep)[1]>0){
          biomasses <- aggregate(years_to_keep$biomass, list(years_to_keep$species), mean) # mean per species
          abundances <-  aggregate(years_to_keep$abundance, list(years_to_keep$species), mean)
          colnames(biomasses) <- c("species","mixture(t/ha)")
          biomasses$abundance <-  abundances$x
          biomasses$'mixture(t/ha)' <- biomasses$'mixture(t/ha)'/(1000*0.08*Nbpatches) # so that the unit becomes t/ha
          biomasses$mixture_relative <- biomasses$'mixture(t/ha)'/sum(biomasses$'mixture(t/ha)') #to have relative biomass
          biomasses$site <- site
          biomasses$order <- order
          biomasses$simul <- number
          BIOMASSES <- rbind(BIOMASSES,biomasses)
        }

      }
    }
  }
  BIOMASSES
}  

biomass_tot_aggreg_mono <- function(site){
  # 1. Reads the data frame printed by specific_biomass_final
  # 2. Sums the specific biomasses for each simul
  # 3. Returns a data.frame (e.g. total_biomass_final_Bern.txt)whose colnames are:
  # "site"	"simul"	"decreasing"	"increasing"	"random_1"	"random_10"	"random_2"	"random_3"	etc.
  # which gives the biomasses (t/ha) for each simul in each order in a given site.
  BIOMASSES_sp <- read.table(paste0("data/processed/Aggregated monocultures_method 2/biomass_specific_",site,"_aggreg_mono.txt"),header=T)
  BIOMASSES_tot <- aggregate(BIOMASSES_sp$mixture.t.ha.,list(site = BIOMASSES_sp$site,order = BIOMASSES_sp$order,simul = BIOMASSES_sp$simul),FUN=sum)# sum of the biomasses of the species
  
  biomass_per_order <- spread(BIOMASSES_tot,order,x) # data frame with each order in column
  biomass_per_order
}