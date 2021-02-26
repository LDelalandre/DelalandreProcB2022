source("R/Common variables.R")
source("R/Before simulations.R")
source("R/comp_fct_dist.R")

library("cowplot")
library("gridExtra")
library("ggsignif")
library("tidyverse")

distinct_tot <- read.table("data/raw/distinctiveness of the species.txt",header = T)


#________________________________________________________________________________
# Fig 2: PCA and bootstrap ####
# i) PCA ####
traits<-read.table("data/Traits of the species_complete.txt",header=T)
c1<-choice_traits_1(traits) # data.frame with the traits of the species
# cprime <- select(c1,-c(A1max,A2))
ACP1<-PCA(c1)
c1$Dim.1 <- ACP1$ind$coord[, 1] 
c1$Dim.2 <- ACP1$ind$coord[, 2] 
c1$SName <- rownames(c1)
c1$Name <- traits$Name
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
  

ggsave("paper_2/PCA.png",plot_pca,height=26,width=33,units="cm",dpi="print")

percent_var <- factoextra::fviz_eig(ACP1, addlabels = TRUE, ylim = c(0, 30))
ggsave(filename = "paper_2/#figur_variance_axis.png",plot = percent_var)

# try to draw circles around groups of species:
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

# And a demonstration of it's use:

dat <- circleFun(c(1,-1),2.3,npoints = 100)
#geom_path will do open circles, geom_polygon will do filled circles
ggplot(dat,aes(x,y)) + geom_path()

dat <- circleFun(c(1,-1),2.3,npoints = 100)
plot_pca +
  geom_path(data=dat,aes(x=x,y=y))
  


# ii) SENSITIVITY ANALYSIS: bootstrap ####
traits<-read.table("data/Traits of the species_complete.txt",header=T)
selected_traits<-choice_traits_1(traits) # data.frame with the traits of the species

perform_boostrap <- function(perform){
  if (perform == T){
    init_dist <- 
      comp_fct_dist(selected_traits) %>% 
      # arrange(desc(Di)) %>% 
      pull(Di)
    
    corstat <- c()
    for (i in c(1:10000)){
      col <- sample(c(1:14),size=14,replace=T) # sample with replacement the traits
      boot_traits <- selected_traits[,col]
      boot_dist <- comp_fct_dist(boot_traits) %>% 
        pull(Di)
      boot_cor <- cor(init_dist,boot_dist,method="spearman")
      corstat <- c(corstat,boot_cor)
    }
    write.table(corstat,"data/processed/bootstrap on Di correlations.txt",row.names=F,col.names=F)
    
  }
}
perform_boostrap(perform = F) # I did it once, and now I use the data frame saved


# corstat <- read.table("data/processed/bootstrap on Di correlations.txt") %>% 
#   pull(V1)
# summary(corstat)
# sd(corstat)
# 
# 
# jpeg("paper_2/histogram bootstrap.jpg", width = 350, height = 350)
# hist(corstat,main="Bootstrap distribution of rho", xlab="Rho",xlim=c(0,1))
# dev.off()

hist_bootstrap <- ggplot(data.frame(corstat),aes(x=corstat)) +
  geom_histogram(col="black",fill="grey",binwidth = 0.05) +
  labs(x="Rho",y="frequency") +
  theme_classic() +
  theme(axis.title=element_text(size=15),axis.text=element_text(size=15))


# Select a subset of the traits
# I remove the abiotic response traits to cold, and the correlation is still high
selected_traits2 <- select(selected_traits,-c(DDMin,WiTX,NTol,Ly,La))
newDi <- comp_fct_dist(selected_traits2) %>% 
  pull(Di)
cor(init_dist,newDi,method="spearman")
# Just give the info in the text 


# Combine the two plots ####
fig2 <- plot_pca +
  annotation_custom(ggplotGrob(hist_bootstrap), 
                    xmin = 2.6, xmax = 5.4, ymin = 1.5, ymax = 3.8) + 
  geom_text(  label="A", x=-4,y=3.9,size = 13,color = "black") +
  geom_text(  label="B", x=3.9,y=3.9,size = 13,color = "black")

ggsave("paper_2/#figur_pca.png",fig2,height=26,width=29,units="cm",dpi="print")



#________________________________________________________________________________
# table traits and species names ####
traits_desc <- read.table("paper_2/Parameters description.txt",sep="\t")
colnames(traits_desc) <- c("Trait","Description")

table4 <- mtcars[1:5, 1:4] %>%
  mutate(
    car = row.names(.),
    mpg = mpg,
    cyl = cyl,
    disp = ifelse(disp > 200,
                  cell_spec(disp, bold = T),
                  cell_spec(disp, bold=F)),
    hp = hp
  ) %>%
  select(car, everything()) %>%
  kable( escape = F) %>%
  kable_styling("hover", full_width = F) %>%
  column_spec(5, width = "3cm")
cat(table4, file = "Test.doc")
  
