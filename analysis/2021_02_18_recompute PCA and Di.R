source("R/Common variables.R")
source("R/Before simulations.R")
source("R/comp_fct_dist.R")

library("cowplot")
library("gridExtra")
library("ggsignif")
library("tidyverse")

# COMMON VAR ####
traits<-read.table("data/Traits of the species_complete.txt",header=T)
selected_traits<-choice_traits_1(traits) # data.frame with the traits of the species used in the main analysis
distinct_tot <- read.table("data/raw/distinctiveness of the species.txt",header = T)

# FUNCTIONs ####
make_PCA <- function(selected,distinct_tot){ 
  # selected is the set of traits (in columns) and species (in rows) from which I make the PCA
  # I can subset the traits, or the species, or both
  # NB I can make a PCA on all the dataset, but recompute Di from a subset of traits and add it on the plot !
  # (I have to be sure to specify it so that I don't mess up)
  
  # make the PCA
  c1<-selected # data.frame with the traits of the species
  
  ACP1<-PCA(selected)
  
  c1$Dim.1 <- ACP1$ind$coord[, 1] 
  c1$Dim.2 <- ACP1$ind$coord[, 2] 
  c1$SName <- rownames(c1)
  c1$Name <- traits %>% 
    filter(SName %in% c1$SName) %>% 
    pull(Name)
  
  c1$Distinctiveness <- distinct_tot$Di
  c1$Distinctiveness <- c1$Distinctiveness/max(c1$Distinctiveness) # to have a relative distinctiveeness (good idea?)
  
  axis <- ACP1$var$coord[,c(1,2)] %>% 
    data.frame() %>% 
    rownames_to_column(var="varnames")
  
  axis2 <- transform(axis,
                     Dim.1 = 4.7 * Dim.1,
                     Dim.2 = 4.7* Dim.2)
  
  var.explain.dim1 <- round(ACP1$eig[1,2],digits=1)
  var.explain.dim2 <- round(ACP1$eig[2,2],digits=1)
  
  plot_pca <-
    ggplot(data=c1,aes(x=Dim.1,y=Dim.2)) + 
    geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
    geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
    coord_equal() + 
    geom_text(data=axis2, aes(x=Dim.1, Dim.2, label=varnames), size = 5, vjust=1, color="black")+
    geom_segment(data=axis2, aes(x=0, y=0, xend=Dim.1-0.2, yend=Dim.2-0.2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
    ggrepel::geom_label_repel(aes(label = Name),size = 5) + # or SName
    aes(color=Distinctiveness) +
    scale_colour_gradient(low="#00BFC4",high="#F8766D") +
    # scale_colour_gradient(low = "#132B43",high = "#56B1F7" ) +
    labs(x=paste0("Dim 1 (",var.explain.dim1,"%)"),
         y=paste0("Dim 2 (",var.explain.dim2,"%)") ) +
    theme_classic() +
    theme(axis.title=element_text(size=15),axis.text=element_text(size=15)) 
  
  plot_pca
}


PCA_and_correlation <- function(selected){ 
  # chose the traits and species of interest in selected
  # be careful not to have changed row (=species) order!
  distinct_new <- comp_fct_dist(selected)
  
  PCA <- make_PCA(selected,distinct_new)
  
  # correlate the two di orders
  A <- distinct_tot %>% 
    filter(SName %in% rownames(distinct_new)) %>% 
    pull(Di)
  B <- distinct_new %>% 
    pull(Di)
  cortest <- cor.test(A,B,method="spearman")
  
  list(PCA,distinct_new,cortest)
}

# CODE ####

# i) without conifers ####
conifers <- c("AAlb","TBac","PCem","PMon","LDec","PSyl","PAbi")

selected <- selected_traits %>%
  filter(!(rownames(.) %in% conifers))

PCA_and_correlation(selected)
# Very correlated : the same part of the functional space is still distinct!

# ii) without A1max ####
selected <- selected_traits %>%
  select(-c(A1max,DDMin,NTol,WiTX))

PCA_and_correlation(selected)


# iii) Annette's suggestion ####
selected <- selected_traits %>%
  select(c(HMax,AMax,G,Brow,Ly,La))

out <- PCA_and_correlation(selected)
newdi <- out[[2]] %>% 
  arrange(desc(Di))
# we keep the order rather well, but we lose PCem which gets common.
