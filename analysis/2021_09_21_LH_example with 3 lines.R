Mi <- c(1.002,0.708,0.015)
YOi <- c(0.155,0.013,0)
nb_sp_regional_pool <- c(3,3,3)
nb_sp_realized_pool <- c(2,2,2)
RYEi <- c(1/3,1/3,1/3)

LH <- data.frame(Mi,YOi,nb_sp_regional_pool,nb_sp_realized_pool,RYEi)

LH2 <- LH %>% 
  mutate(YO=sum(YOi)) %>%
  mutate(RYEi=1/nb_sp_regional_pool) %>% 
  
  mutate(RYOi=YOi/Mi) %>%
  mutate(YEi=RYEi*Mi) %>%
  mutate(YE=sum(YEi)) %>%
  mutate(DeltaY=YO-YE) %>%
  
  mutate(DeltaRYi=RYOi-RYEi)
  

#_________________________________________________
LH2 %>% 
  summarize(Selection = nb_sp_regional_pool * ( mean(DeltaRYi * Mi) - mean(DeltaRYi) * mean(Mi) ),
            Complementarity = nb_sp_regional_pool * mean(DeltaRYi)*mean(Mi),
            DeltaY = nb_sp_regional_pool * mean(DeltaRYi * Mi)) %>%  # idem as YO - YE
  mutate(sum = Selection + Complementarity)

#___________________________________________________
    # compute terms of decompo
LH2 %>% mutate(Cpltarity = 3*DeltaRYavg*Mavg) %>%
  mutate(Selection_diff = DeltaY - Cpltarity)%>%
  mutate(Selection_cov=3*cov(DeltaRYi,Mi)) %>% 
  mutate(DeltaRYM = DeltaRYi*Mi) %>% 
  mutate(DeltaY2 = sum(DeltaRYM)) %>% 
  mutate(covar = 1/3*sum(DeltaRYi*Mi) - 1/3*sum(DeltaRYi)*1/3*sum(Mi)) %>%  # compute covariance
  mutate(selection_cov_manual = 3*covar)
  
LH2 %>% 
  summarize(DeltaY=mean(DeltaY),
            Cpltarity=mean(Cpltarity),
            Selection=mean(Selection),
            Selection2=mean(Selection2),
            DeltaY2=mean(DeltaY2)
            ) %>% 
  mutate(sum = Cpltarity + Selection2)

  
cov(c(-0.178,-0.315,-1/3),c(1.002,0.708,0.015))

0.034*3 - 0.159*3

# Je veux deltay= -0.502