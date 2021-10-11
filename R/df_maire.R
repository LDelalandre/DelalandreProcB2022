# Generate data on which Maire's framework can be applied to identify key species 
# for productivity, selection, and complementarity.


df_maire <- function(LH_all_per_sp,sit,simul_to_keep){
  # LH_all_per_sp provides Loreau-Hector's partitioning of biodiversity effects and its components at the species level.
  df_LH <- LH_all_per_sp %>% 
    filter(site == sit)
  
  dat <- data.frame(matrix(nrow=0,ncol=37))
  colnames(dat) <- c("site","order","simul",
                     as.character(read.table("data/raw/distinctiveness of the species.txt",header = T)$SName),
                     "productivity","DeltaY","Selection","Complementarity")
  j <- 1 # number of the assemblage = position (line) in the data.frame dat
  
  for (order in ORDER[3:32]){
    a <- order
    df_LH_order <- subset(df_LH,order==a)
    
    i = simul_to_keep # keep only simul 16, with 15 species
    sub <- subset(df_LH_order, simul==i) 
    # filter(persists_mixt == T) # NB SPECIES THAT ARE PRESENT IN THE SIMUL ARE ABOVE BIOMASS THRESHOLD
    # NB2: I did not filter on species persisting in mono and mixt to performe Loreau-Hector analysis.
    
    # dat[j,]$assemblage <- paste0(sit,"_",order,"_",i) # paste order and simul (an assemblage is a simul)
    dat[j,]$site <- sit
    dat[j,]$order <- order
    dat[j,]$simul <- i
    
    dat[j,]$productivity <- sum(sub$YOi)
    dat[j,]$DeltaY <- unique(sub$DeltaY)
    dat[j,]$Selection <- unique(sub$Selection)
    dat[j,]$Complementarity <- unique(sub$Complementarity)
    
    # all the species that we saw at the end of the simul are considered persistent (bad)
    # spnames <- as.character(sub$SName) 
    
    # only species above biomass threshold are considered persistent (better, I think)
    spnames <- sub %>% 
      pull(SName) %>% 
      as.character
    
    
    
    dat[j,which(!(colnames(dat) %in% c(spnames,"site","order","simul","productivity","DeltaY","Selection","Complementarity") ))] <- 0
    dat[j,which(colnames(dat) %in% sub$SName)] <- 1
    
    j <- j+1 # current position in the data.frame dat.
  }
  dat
} # end of df_maire