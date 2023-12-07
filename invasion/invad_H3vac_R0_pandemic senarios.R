rm(list=ls())
library(deSolve)
load('H3Vac.invad.RData')

sirmod <- function(t,y,parms){
  with(as.list(c(y, parms)), {
    beta0_1 <- R0_1*(gamma1+mu)
    beta0_2 <- R0_2*(gamma2+mu)
    
    beta1 <- beta0_1*(1+a1*cos(2*pi*t)) #seasonality forcing
    beta2 <- beta0_2*(1+a2*cos(2*pi*t)) #seasonality forcing
    
    dS <- mu-(rho+mu)*S-beta1*S*(I1+I21+IV1)-beta2*S*(I2+I12+IV12)+
      sigma1*R1+sigma2*R2+sigmaV*V #susceptible
    
    dI1 <- beta1*S*(I1+I21+IV1)-(gamma1+mu)*I1 #primary Strain1 infection
    dRT1 <- gamma1*I1+sigmaV2*RVT1-(mu+sigmaT+rho2)*RT1 #strain-transcending immunized + infection-induced immunized against strain1
    dR1 <- sigma2*R12 +sigmaT*RT1 +sigmaV2*RV1 -
      beta2*(1-theta1)*R1*(I2+I12+IV12) -(mu+sigma1+rho2)*R1 #infection-induced immunized against strain1
    dI12 <- beta2*(1-theta1)*R1*(I2+I12+IV12)-(mu+gamma2)*I12 #recovered from primary strain1 infection and infected by strain2 
    
    dI2 <- beta2*S*(I2+I12+IV12)-(gamma2+mu)*I2 #primary Strain2 infection
    # dRT2 <- gamma2*I2 + sigmaV1*RVT2 - (mu+sigmaT+rho1)*RT2 #strain-transcending immunized + infection-induced immunized against strain2
    # dR2 <- sigma1*R12 +sigmaT*RT2 +sigmaV1*RV2 - 
    #   beta1*(1-theta2)*R2*(I1+I21+IV1+IV21) -(mu+sigma2+rho1)*R2 #infection-induced immunized against strain2
    
    dRT2 <- gamma2*I2 - (mu+sigmaT)*RT2 #strain-transcending immunized + infection-induced immunized against strain2
    dR2 <- sigma1*R12 +sigmaT*RT2 - beta1*(1-theta2)*R2*(I1+I21+IV1) -
      (mu+sigma2)*R2 #infection-induced immunized against strain2
    
    dI21 <- beta1*(1-theta2)*R2*(I1+I21+IV1)-(mu+gamma1)*I21 #recovered from primary strain2 infection and infected by strain1
    
    # dV <- rho*S +sigma1*RV1+sigma2*RV2 -(mu+sigmaV)*V-
    #   beta2*(1-tau)*V*(I2+I12+IV2+IV12)-beta1*(1-tau2)*V*(I1+I21+IV1+IV21) #vaccinated
    
    dV <- rho*S +sigma1*RV1 -(mu+sigmaV)*V-beta1*(1-tau2)*V*(I1+I21+IV1) #vaccinated
    
    
    dIV1 <- beta1*(1-tau2)*V*(I1+I21+IV1) - (gamma1+mu)*IV1
    dRVT1 <- gamma1*IV1 + rho2*RT1 - (sigmaV2+sigmaT+mu)*RVT1
    dRV1 <- sigma2*RV12 +sigmaT*RVT1 + rho2*R1-
      beta2*(1-theta1)*(1-tau)*RV1*(I2+I12+IV12)-(sigma1+sigmaV2+mu)*RV1
    dIV12 <- beta2*(1-theta1)*(1-tau)*RV1*(I2+I12+IV12)-
      gamma2*IV12 -mu*IV12
    
    # dIV2 <- beta2*(1-tau)*V*(I2+I12+IV2+IV12)-(gamma2+mu)*IV2 #vaccinated individuals infected by strain2
    # dRVT2 <- gamma2*IV2 + rho1*RT2 - (sigmaV1+sigmaT+mu)*RVT2 #strain-transcending immunized - recovered form strain2 infection
    # dRV2 <- sigma1*RV12 +sigmaT*RVT2 +rho1*R2 -(sigma2+sigmaV1+mu)*RV2 -
    #   beta1*(1-theta2)*(1-tau2)*RV2*(I1+I21+IV1+IV21) #vaccinated individuals recovered from strain2 infection 
    
    # dIV21 <- beta1*(1-theta2)*(1-tau2)*RV2*(I1+I21+IV1+IV21) -(gamma1+mu)*IV21
    # dRV12 <- gamma1*IV21 + gamma2*IV12-(sigma1+sigma2+mu)*RV12
    dRV12 <- gamma2*IV12-(sigma2+mu)*RV12
    
    dR12 <- gamma2*I12 +gamma1*I21 - (mu+sigma1+sigma2)*R12
    
    
    res <- c(dS,dI2,dV,dR2,dRT2,dI1,dR12,dR1,dI12,dI21,dRT1,dRV12,dRV1,dIV1,
             dRVT1,dIV12)
    
    return(list(res,inci1=beta1*S*(I1+I21+IV1)+beta1*(1-theta2)*R2*(I1+I21+IV1)+
                  beta1*(1-tau2)*V*(I1+I21+IV1),
                inci2=beta2*S*(I2+I12+IV12)+beta2*(1-theta1)*R1*(I2+I12+IV12)))
  })
}
times.invad <- seq(0,5/52,by=1/(52*7))
invad.list.R0 <- lapply(c(1,1.3,1.6,3.2),function (r){
  invad.list.tau2 <- lapply(seq(0,1,by=0.1),function(tau2){
    parms.invad <- parms0
    parms.invad['tau2'] <- tau2
    # natural immunity aginst seasonal variant
    # parms.invad['theta2']=0.5
    # pandemic strain
    parms.invad['theta2']<-0
    parms.invad['theta1']<-0
    parms.invad['R0_1'] <- r
    parms.invad['gamma1'] <- parms.invad['gamma2']
    parms.invad['a1'] <- parms.invad['a2']
    parms.invad['beta0_1']<-r*(parms.invad['gamma2']+parms.invad['mu'])
    parms.invad['sigma1'] <- parms.invad['sigma2']
    invad.list <- lapply(simdf.list, function(simdf) {
      start.invad <- unlist(simdf[nrow(simdf),2:(ncol(simdf)-5)])
      start.invad <- start.invad[-5]
      start.invad['I1'] <- 1e-6
      start.invad['S'] <- start.invad['S']-1e-6
      start.invad['R12'] <- 0
      start.invad['R1'] <- 0
      start.invad['I12'] <-0
      start.invad['I21'] <-0
      start.invad['RT1'] <-0
      start.invad['RV12'] <-0
      start.invad['RV1'] <-0
      start.invad['IV1'] <-0
      start.invad['RVT1'] <-0
      # start.invad['IV21'] <-0
      start.invad['IV12'] <-0
      parms.invad['rho'] <- simdf[1,'H3vac.cov']
      parms.invad['sigmaV'] <- simdf[1,'H3vac.sigmaV']
      parms.invad['sigmaV1'] <- 0
      parms.invad['sigmaV2'] <- parms.invad['sigmaV']
      parms.invad['rho1'] <-0
      parms.invad['rho2'] <- parms.invad['rho']
      
      # if(start.invad['I2'] < 1e-100){
      #   start.invad['I2'] <- 0
      #   start.invad['R2'] <- 0
      #   start.invad['RT2'] <- 0
      # }
      invad.output <- as.data.frame(rk(y=start.invad, times=times.invad,
                                       func=sirmod,parms=parms.invad,hmax=1/365))

      invad.output['H3vac.cov']<-parms.invad['rho']
      invad.output['H3vac.tau2']<-parms.invad['tau2']
      invad.output['H3vac.sigmaV']<-parms.invad['sigmaV']
      invad.output['R0_1'] <- parms.invad['R0_1']
      return(invad.output)
    })
    return(invad.list)
  })
  return(invad.list.tau2)
})
save.image('H3endemicVac.pandStrInvad.NoNatImm.RData')

## 5 years trajectory
times.invad.5yr <- seq(0,5,by=1/(52*7))
invad.list.R0 <- lapply(c(1,1.3,1.6,3.2),function (r){
  invad.list.tau2 <- lapply(seq(0,1,by=0.1),function(tau2){
    parms.invad <- parms0
    parms.invad['tau2'] <- tau2
    # natural immunity aginst seasonal variant
    # parms.invad['theta2']=0.5
    # pandemic strain
    parms.invad['theta2']<-0
    parms.invad['theta1']<-0
    parms.invad['R0_1'] <- r
    parms.invad['gamma1'] <- parms.invad['gamma2']
    parms.invad['a1'] <- parms.invad['a2']
    parms.invad['beta0_1']<-r*(parms.invad['gamma2']+parms.invad['mu'])
    parms.invad['sigma1'] <- parms.invad['sigma2']
    invad.list <- lapply(simdf.list, function(simdf) {
      start.invad <- unlist(simdf[nrow(simdf),2:(ncol(simdf)-5)])
      start.invad <- start.invad[-5]
      start.invad['I1'] <- 1e-6
      start.invad['S'] <- start.invad['S']-1e-6
      start.invad['R12'] <- 0
      start.invad['R1'] <- 0
      start.invad['I12'] <-0
      start.invad['I21'] <-0
      start.invad['RT1'] <-0
      start.invad['RV12'] <-0
      start.invad['RV1'] <-0
      start.invad['IV1'] <-0
      start.invad['RVT1'] <-0
      # start.invad['IV21'] <-0
      start.invad['IV12'] <-0
      parms.invad['rho'] <- simdf[1,'H3vac.cov']
      parms.invad['sigmaV'] <- simdf[1,'H3vac.sigmaV']
      parms.invad['sigmaV1'] <- 0
      parms.invad['sigmaV2'] <- parms.invad['sigmaV']
      parms.invad['rho1'] <-0
      parms.invad['rho2'] <- parms.invad['rho']
      
      # if(start.invad['I2'] < 1e-100){
      #   start.invad['I2'] <- 0
      #   start.invad['R2'] <- 0
      #   start.invad['RT2'] <- 0
      # }
      invad.output <- as.data.frame(rk(y=start.invad, times=times.invad.5yr,
                                       func=sirmod,parms=parms.invad,hmax=1/365))
      
      invad.output['H3vac.cov']<-parms.invad['rho']
      invad.output['H3vac.tau2']<-parms.invad['tau2']
      invad.output['H3vac.sigmaV']<-parms.invad['sigmaV']
      invad.output['R0_1'] <- parms.invad['R0_1']
      return(invad.output)
    })
    return(invad.list)
  })
  return(invad.list.tau2)
})
save.image('H3endemicVac.pandStrInvad5yrs.NoNatImm.RData')



# pandSize.R0 <- lapply(invad.list.R0, function(x){
#   invad.list <- x[[1]]
#   trough <- ldply(invad.list,rbind) %>%
#     filter(time==1)
#   trough$R0 <- x[[2]]
#   return (trough)
# })

# pandSize.R0 <- ldply(invad.list.R0, function(x){
#   invad.list <- lapply(x, function(df){
#     a <- ldply(df,rbind)
#     return(a)
#   })
#   trough <- ldply(invad.list,rbind) %>%
#     filter(time==5)
#   return(trough)
# })

library(plyr)
library(dplyr)
trough.R0 <- ldply(invad.list.R0, function(x){
  invad.list <- lapply(x, function(df){
    a <- ldply(df,rbind)
    return(a)
  })
  
  # traj <-ldply(invad.list, function(x){
  #   simdf <- ldply(x,rbind)
  # })
  traj <-ldply(invad.list, rbind)
  # traj <- ldply(invad.list,rbind)
  
  trough <- traj %>%
    group_by(H3vac.cov,H3vac.tau2,H3vac.sigmaV,R0_1) %>%
    summarise(peak.time1=which.max(inci1),trough1=min(inci1),
              peak.time2=which.max(inci2),trough2=min(inci2),
              Itotal1=sum(inci1),
              Itotal2=sum(inci2))
  
  trough[trough$peak.time1==1,'trough1'] <- -1
  trough[trough$peak.time2==1,'trough2'] <- -1
  
  # trough$R0 <- x[[2]]
  
  return (trough)
})
troughR0 <- trough.R0 %>%
  mutate(immDur=1/H3vac.sigmaV,vac.rate=H3vac.cov/52*100)
save(troughR0,pandSize.R0,file='invasion/plot.data.RData')
library(ggplot2)
pl.troughR0 <- ggplot(troughR0,aes(x=vac.rate,y=H3vac.tau2,z=trough1)) +
  facet_grid(immDur~R0_1,scales='free',
             labeller = label_bquote(rows=sigma[V]:.(1/immDur), 
                                     cols=R[0]^1:.(R0_1)))+
  theme_classic()+
  theme(panel.spacing = unit(0.6, "lines"))+
  geom_contour_filled(bins=4,breaks=c(-2,0,10e-8,10e-6,10e-4)) +
  geom_contour(bins=4,breaks=c(-2,0,10e-8,10e-6,10e-4),
               aes(colour = factor(..level..==10e-6,levels = c(T, F),
              labels = c("Pandemic persistence","")))) +
  scale_colour_manual(values = c("black", "#00000000")) +
  labs(fill='Trough depth',color='')+
  xlab('Vaccination rate (%) per week')+
  scale_x_continuous(expand = c(0,0),labels=c('0','0.25','0.5','0.75','1'))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values=c('grey','#FFF4BD',
                             '#feb24c','#bd0026'),
                    labels=list("No invasion",expression(0-10^{-8}), 
                                expression(10^{-8}-10^{-6}),
                                expression(10^{-6}-10^{-4})))+
  labs(y=expression('susceptibility reduction by vaccines'~tau))
pl.troughR0
ggsave('H3endemicVac.pandStrR0varyRhoInvad.NoNatImm.trough.png',pl.troughR0, width=7,height=6)

# trough.R0big <- troughR0 %>%
#   filter(R0==3.2)
# 
# trough.R0big.rho <- trough.R0big %>%
#   mutate(subpl="vaccine immunity duration (yr)", 
#          pl.y=immDur, grp=vac.rate) %>%
#   filter(grp==1.730)
# 
# trough.R0big.sigmaV <- trough.R0big %>%
#   mutate(subpl="vaccination rate (%) per week", 
#          pl.y=vac.rate, grp=immDur) %>%
#   filter(grp==16)
# 
# trough.R0big.pl <- trough.R0big.sigmaV %>%
#   rbind(trough.R0big.rho)
# 
# pl <- ggplot(trough.R0big.pl,aes(x=vac.tau2,y=pl.y,z=trough)) +
#   facet_wrap(subpl~.,scales='free',nrow=2,strip.position = "left")+
#   theme_classic()+theme(strip.placement="outside",strip.background = element_blank(),
#                         strip.text.y = element_text(size = 10, angle = -90,colour = 'black'))+
#   geom_contour_filled(bins=10)+
#   labs(fill='Trough depth',color='')+
#   ylab(NULL)+
#   scale_x_continuous(expand = c(0,0),labels=c('0','0.25','0.5','0.75','1'))+
#   scale_y_continuous(expand = c(0,0))+
#   # scale_fill_manual(values=rev(heat.colors(12)))+
#   scale_fill_manual(values=c('#fed976','#feb24c','#e31a1c','#bd0026','#800026',
#   '#CCBAD7FF', '#BA9FC7FF','#A783B6FF','#8F5EA1FF','black'))+
#   # '#773A8BFF','#572872FF','#361959FF','black'))+
#   # '#e31a1c'
#   # ,'#bd0026','#800026'))+
#   # scale_color_manual(values=c('#1a9850','#ffffcc','#ffeda0','#fed976','#feb24c',
#   #                             'black','#fc4e2a','#e31a1c','#bd0026','#800026'))+
#   labs(x=expression('susceptibility reduction by vaccines'~tau))
# pl
# ggsave('H3endemicVac.pandStrR0_3.2_Invad_pandSize.NoNatImm.png',pl, width=6,height=6)
# 
# trough.R0big.sigmaV %>%
#   filter(vac.rate==1.730)
# 
traj.R0 <- ldply(invad.list.R0, function(x){
  invad.list <- lapply(x, function(df){
    a <- ldply(df,rbind)
    return(a)
  })
  # traj <-ldply(invad.list, function(x){
  #   simdf <- ldply(x,rbind)
  # })
  # traj <-ldply(invad.list, rbind)
  trough <- ldply(invad.list,rbind) %>%
    mutate(immDur=1/H3vac.sigmaV,vac.rate=signif(H3vac.cov/52*100,3)) %>%
    filter(immDur==2 & vac.rate==1.730)

  trough2 <- ldply(invad.list,rbind) %>%
    mutate(immDur=1/H3vac.sigmaV,vac.rate=signif(H3vac.cov/52*100,3)) %>%
    filter(immDur==2 & vac.rate==0 & H3vac.tau2==0.2)

  trough3 <- trough %>%
    rbind(trough2)
  # trough3$R0 <- x[[2]]
  return (trough3)
})
# 
traj <- traj.R0 %>% filter(R0_1==1.6)
# 
pl.traj <- ggplot(data=traj %>% filter(H3vac.tau2==0.2|H3vac.tau2==0.4|H3vac.tau2==0.8),mapping = aes(x=time))+
  geom_line(aes(y=inci1,lty=factor(H3vac.tau2)))+
  facet_wrap(~vac.rate,labeller = label_bquote(cols=rho:.(vac.rate))) +
  theme_classic() +
  scale_y_log10()+
  labs(lty='susceptibility reduction\nby vaccines',
       y = 'Incidence of the pandemic strain')
ggsave('invasion/traj.png',pl.traj,width=6,height=3)
# 
# 
# traj.R0.2 <- lapply(invad.list.R0.5yr, function(x){
#   invad.list <- x[[1]]
#   trough <- ldply(invad.list,rbind) %>%
#     mutate(immDur=1/vac.sigmaV,vac.rate=signif(vac.cov/52*100,3)) %>%
#     filter(vac.tau2==0.2 & vac.rate==0)
#   trough$R0 <- x[[2]]
#   return (trough)
# })
# 
# traj <- traj.R0.2[[4]]
# 
# pl.traj <- ggplot(data=traj,mapping = aes(x=time))+
#   geom_line(aes(y=inci1,lty=factor(immDur)))+
#   theme_classic() +
#   scale_y_log10()+
#   labs(lty='susceptibility reduction\nby vaccines',
#        y = 'Incidence of invading strain')
# 
# library(ggplot2)
# pandSize.R0 <- pandSize.R0 %>%
#   mutate(immDur=1/H3vac.sigmaV,vac.rate=H3vac.cov/52*100,
#          pandSize=I1+I21+IV1)
# 
# pl.pandSize.R0 <- ggplot(pandSize.R0,aes(x=vac.rate,y=H3vac.tau2,z=pandSize*1000)) +
#   
#   
#   ggplot(trough.R0,aes(x=vac.rate,y=H3vac.tau2,z=pandSize*1000)) +
#   facet_grid(immDur~R0_1,scales='free',labeller = label_bquote(rows=sigma[V]:.(1/immDur), 
#                                                                cols=R[0]^1:.(R0_1)))+
#   theme_classic()+
#   theme(panel.spacing = unit(0.6, "lines"))+
#   geom_contour_filled(bins=5)+
#   labs(fill='Pandemic attack \nrate (x10e-3)',color='')+
#   xlab('Vaccination rate (%) per week')+
#   scale_x_continuous(expand = c(0,0),labels=c('0','0.25','0.5','0.75','1'))+
#   scale_y_continuous(expand = c(0,0))+
#   # scale_fill_manual(values=rev(heat.colors(12)))+
#   scale_fill_manual(values=c('#fed976','#feb24c','#e31a1c','#bd0026','#800026'))+
#   labs(y=expression('susceptibility reduction by vaccines'~tau))
# pl.pandSize.R0
# ggsave('invasion/H3endemicVac.pandStrR0varyRhoInvad.NoNatImm.PandSize.png',pl.pandSize.R0, width=7,height=6)

# pandSize.R0big <- pandSize %>%
#   filter(R0==3.2)
# 
# pandSize.R0big.rho <- pandSize.R0big %>%
#   # filter(vac.cov==unique(simdf.t3$vac.cov)[2]) %>%
#   mutate(subpl="vaccine immunity duration (yr)", 
#          pl.y=immDur, grp=vac.rate) %>%
#   # filter(R0==1|R0==1.4|R0==1.8|R0==2.2|R0==2.6|R0==3) %>%
#   # filter(grp==0 | grp==0.577 | grp==1.15 | grp==1.73)
#   filter(grp==1.730)
# 
# pandSize.R0big.sigmaV <- pandSize.R0big %>%
#   # filter(vac.sigmaV==unique(simdf.t3$vac.sigmaV)[3]) %>%
#   mutate(subpl="vaccination rate (%) per week", 
#          pl.y=vac.rate, grp=immDur) %>%
#   # filter(grp==1 | grp==4 | grp==16 | grp==64)
#   filter(grp==16)
# pandSize.R0big.pl <- pandSize.R0big.sigmaV %>%
#   rbind(pandSize.R0big.rho)
# 
# pl <- ggplot(pandSize.R0big.pl,aes(x=vac.tau2,y=pl.y,z=Itotal1)) +
#   facet_wrap(subpl~.,scales='free',nrow=2,strip.position = "left")+
#   theme_classic()+theme(strip.placement="outside",strip.background = element_blank(),
#                         strip.text.y = element_text(size = 10, angle = -90,colour = 'black'))+
#   geom_contour_filled(bins=5)+
#   labs(fill='Pandemic attack rate',color='')+
#   ylab(NULL)+
#   scale_x_continuous(expand = c(0,0),labels=c('0','0.25','0.5','0.75','1'))+
#   scale_y_continuous(expand = c(0,0))+
#   # scale_fill_manual(values=rev(heat.colors(12)))+
#   scale_fill_manual(values=c('#fed976','#feb24c','#e31a1c','#bd0026','#800026'))+
#   # '#CCBAD7FF', '#BA9FC7FF','#A783B6FF','#8F5EA1FF',
#   # '#773A8BFF','#572872FF','#361959FF','black'))+
#   # '#e31a1c'
#   # ,'#bd0026','#800026'))+
#   # scale_color_manual(values=c('#1a9850','#ffffcc','#ffeda0','#fed976','#feb24c',
#   #                             'black','#fc4e2a','#e31a1c','#bd0026','#800026'))+
#   labs(x=expression('susceptibility reduction by vaccines'~tau))
# 
# ggsave('H3endemicVac.pandStrR0_3.2_Invad_pandSize.NoNatImm.png',pl, width=6,height=6)
