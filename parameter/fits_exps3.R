rm(list=ls())
library(deSolve)
library(dplyr)
library(gridExtra)
library(lubridate)
library(ggplot2)
theme_set(theme_bw())

sirmod <- function(t,y,parms){
  with(as.list(c(y, parms)), {
    beta0_1 <- R0_1*(gamma1+mu)
    beta0_2 <- R0_2*(gamma2+mu)
    
    beta1 <- beta0_1*(1+a1*cos(2*pi*t)) #seasonality forcing
    beta2 <- beta0_2*(1+a2*cos(2*pi*t)) #seasonality forcing
    
    dS <- mu-(rho+mu)*S-beta1*S*(I1+I21+IV1+IV21)-beta2*S*(I2+I12+IV2+IV12)+
      sigma1*R1+sigma2*R2+sigmaV*V #susceptible
    
    dI1 <- beta1*S*(I1+I21+IV1+IV21)-(gamma1+mu)*I1 #primary Strain1 infection
    dRT1 <- gamma1*I1+sigmaV2*RVT1-(mu+sigmaT+rho2)*RT1 #strain-transcending immunized + infection-induced immunized against strain1
    dR1 <- sigma2*R12 +sigmaT*RT1 +sigmaV2*RV1 -
      beta2*(1-theta1)*R1*(I2+I12+IV2+IV12) -(mu+sigma1+rho2)*R1 #infection-induced immunized against strain1
    dI12 <- beta2*(1-theta1)*R1*(I2+I12+IV2+IV12)-(mu+gamma2/(1-nu1))*I12 #recovered from primary strain1 infection and infected by strain2 
    
    dI2 <- beta2*S*(I2+I12+IV2+IV12)-(gamma2+mu)*I2 #primary Strain2 infection
    dRT2 <- gamma2*I2 + sigmaV1*RVT2 - (mu+sigmaT+rho1)*RT2 #strain-transcending immunized + infection-induced immunized against strain2
    dR2 <- sigma1*R12 +sigmaT*RT2 +sigmaV1*RV2 - 
      beta1*(1-theta2)*R2*(I1+I21+IV1+IV21) -(mu+sigma2+rho1)*R2 #infection-induced immunized against strain2
    dI21 <- beta1*(1-theta2)*R2*(I1+I21+IV1+IV21)-(mu+gamma1/(1-nu2))*I21 #recovered from primary strain2 infection and infected by strain1
    
    dV <- rho*S +sigma1*RV1+sigma2*RV2 -(mu+sigmaV)*V-
      beta2*(1-tau)*V*(I2+I12+IV2+IV12)-beta1*(1-tau2)*V*(I1+I21+IV1+IV21) #vaccinated
    
    dIV1 <- beta1*(1-tau2)*V*(I1+I21+IV1+IV21) - (gamma1/(1-eta2)+mu)*IV1
    dRVT1 <- gamma1/(1-eta2)*IV1 + rho2*RT1 - (sigmaV2+sigmaT+mu)*RVT1
    dRV1 <- sigma2*RV12 +sigmaT*RVT1 + rho2*R1-
      beta2*(1-theta1)*(1-tau)*RV1*(I2+I12+IV2+IV12)-(sigma1+sigmaV2+mu)*RV1
    dIV12 <- beta2*(1-theta1)*(1-tau)*RV1*(I2+I12+IV2+IV12) -
      gamma2/(1-nu1)/(1-eta)*IV12 -mu*IV12
    
    dIV2 <- beta2*(1-tau)*V*(I2+I12+IV2+IV12)-(gamma2/(1-eta)+mu)*IV2 #vaccinated individuals infected by strain2
    dRVT2 <- gamma2/(1-eta)*IV2 + rho1*RT2 - (sigmaV1+sigmaT+mu)*RVT2 #strain-transcending immunized - recovered form strain2 infection
    dRV2 <- sigma1*RV12 +sigmaT*RVT2 +rho1*R2 -(sigma2+sigmaV1+mu)*RV2 -
      beta1*(1-theta2)*(1-tau2)*RV2*(I1+I21+IV1+IV21) #vaccinated individuals recovered from strain2 infection 
    
    dIV21 <- beta1*(1-theta2)*(1-tau2)*RV2*(I1+I21+IV1+IV21) -
      (gamma1/(1-eta2)/(1-nu2)+mu)*IV21
    dRV12 <- gamma1/(1-eta2)/(1-nu2)*IV21 + gamma2/(1-nu1)/(1-eta)*IV12-
      (sigma1+sigma2+mu)*RV12
    dR12 <- gamma2/(1-nu1)*I12 +gamma1/(1-nu2)*I21 - (mu+sigma1+sigma2)*R12
    

    res <- c(dS,dI1,dI2,dV,dR12,dR1,dR2,dI12,dI21,dIV2,dIV12,dRT1,dRT2,
             dRV12,dRV1,dRV2,dRVT2,dIV1,dRVT1,dIV21)
    
    return(list(res,inci1=beta1*S*(I1+I21+IV1+IV21)+beta1*(1-theta2)*R2*(I1+I21+IV1+IV21)+
                  beta1*(1-tau2)*V*(I1+I21+IV1+IV21)+beta1*(1-theta2)*(1-tau2)*RV2*(I1+I21+IV1+IV21),
                inci2=beta2*S*(I2+I12+IV2+IV12)+beta2*(1-theta1)*R1*(I2+I12+IV2+IV12)+
                  beta2*(1-tau)*V*(I2+I12+IV2+IV12)+beta2*(1-theta1)*(1-tau)*RV1*(I2+I12+IV2+IV12)))
  })
}
sim <- function(rel.vac.rho,rel.parm,odefunc,start.cond,times.vec,parms.init){
  
  # strength/duration against both target and non-target clades are scaled
  parms.init['tau']<-parms.init['tau']*rel.parm['rel.tau']
  parms.init['tau2']<-parms.init['tau2']*rel.parm['rel.tau']
  parms.init['sigmaV']<-parms.init['sigmaV']*rel.parm['rel.sigmaV']
  parms.init['sigmaV1']<-parms.init['sigmaV1']*rel.parm['rel.sigmaV']
  parms.init['sigmaV2']<-parms.init['sigmaV2']*rel.parm['rel.sigmaV']
  
  if(parms.init['tau']>1){
    parms.init['tau']<-1
    parms.init['rho1'] <- 0
    parms.init['sigmaV1'] <- 0
  }
  
  if(parms.init['tau2']>1){
    parms.init['tau2']<-1
    parms.init['rho2'] <- 0
    parms.init['sigmaV2'] <- 0
  }
  
  parms.init['rho']<-rel.vac.rho*parms.init['rho']
  parms.init['rho1'] <- rel.vac.rho*parms.init['rho1']
  parms.init['rho2'] <- rel.vac.rho*parms.init['rho2']
  
  SIR.output <- as.data.frame(rk4(y=start.cond, times=times.vec,func=odefunc,
                                  parms=parms.init))
  
  sim.output <- SIR.output %>%
    mutate(time=time+1922,vac.cov=parms.init['rho'],vac1.cov=parms.init['rho1'],
           vac2.cov=parms.init['rho2'],vac.sigma=parms.init['sigmaV'],
           vac1.sigma=parms.init['sigmaV1'],vac2.sigma=parms.init['sigmaV2'],
           vac.tau1=parms.init['tau'],vac.tau2=parms.init['tau2'],
           rel.tau=rel.parm['rel.tau'],rel.sigmaV=rel.parm['rel.sigmaV'])
  
  return(sim.output)
}
sumry.cycl <- function(simdf){
    stats <- data.frame(strain=c('H3','H1'),avg.inci=0,epi.size=0,no.epi=0,
                        avg.epi.size=0,avg.max=0)
    
    stats[stats$strain=='H1','avg.inci']=mean(simdf$inci1)
    stats[stats$strain=='H3','avg.inci']=mean(simdf$inci2)
    stats[stats$strain=='H1','global.min']=min(simdf$inci1)
    stats[stats$strain=='H1','global.max']=max(simdf$inci1)
    stats[stats$strain=='H3','global.min']=min(simdf$inci2)
    stats[stats$strain=='H3','global.max']=max(simdf$inci2)
    
    simdf <- simdf %>%
      mutate(local.max1 = if_else(lag(inci1) < inci1 & lead(inci1) < inci1, TRUE, FALSE),
             local.max2 = if_else(lag(inci2) < inci2 & lead(inci2) < inci2, TRUE, FALSE))
    
    inci1.maxs <- simdf[which(simdf$local.max1),'inci1']
    inci2.maxs <- simdf[which(simdf$local.max2),'inci2']
    inci1.maxs <- inci1.maxs[which(inci1.maxs>mean(simdf$inci1))]
    inci2.maxs <- inci2.maxs[which(inci2.maxs>mean(simdf$inci2))]
    
    stats[stats$strain=='H1','epi.size']<-sum(inci1.maxs)
    stats[stats$strain=='H1','no.epi']<-length(inci1.maxs)
    stats[stats$strain=='H1','avg.max']<- if_else(length(inci1.maxs)==0,
                                                  0,sum(inci1.maxs)/length(inci1.maxs))
    stats[stats$strain=='H1','avg.epi.size']<- if_else(length(inci1.maxs)==0,
                                                       stats[stats$strain=='H1','avg.inci'],
                                                       sum(inci1.maxs)/length(inci1.maxs))
    
    stats[stats$strain=='H3','epi.size']<-sum(inci2.maxs)
    stats[stats$strain=='H3','no.epi']<-length(inci2.maxs)
    stats[stats$strain=='H3','avg.max']<- if_else(length(inci2.maxs)==0,
                                                  0,sum(inci2.maxs)/length(inci2.maxs))
    stats[stats$strain=='H3','avg.epi.size']<- if_else(length(inci2.maxs)==0,
                                                       stats[stats$strain=='H3','avg.inci'],
                                                       sum(inci2.maxs)/length(inci2.maxs))
    
    stats['vac.cov']=simdf[1,'vac.cov']
    stats['vac.sigma']=simdf[1,'vac.sigma']
    stats['vac.tau1']=simdf[1,'vac.tau1']
    stats['vac.tau2']=simdf[1,'vac.tau2']
    stats['rel.tau']=simdf[1,'rel.tau']
    stats['rel.sigmaV']=simdf[1,'rel.sigmaV']
    
    return(stats)
  }


bestfitsparms <- data.frame()
bestfitsparms[1,'theta1'] <- 0.8
bestfitsparms[1,'theta2'] <- 0.5
bestfitsparms[1,'a'] <- 0.04
bestfitsparms[2,'theta1'] <- 0.55
bestfitsparms[2,'theta2'] <- 0.45
bestfitsparms[2,'a'] <- 0.04
for (idx in 1:nrow(bestfitsparms)){
  idx = 2
  # initial parameter (w/ vaccine)
  parms0 <- c(R0_1=1.44,R0_2=1.60,a1=bestfitsparms[idx,'a'],
              a2=bestfitsparms[idx,'a'],sigma1=1/3.12,sigma2=1/2.28,
              sigmaV1=1/3.12,sigmaV2=1/2.28,sigmaV=1/2.7,
              gamma1=365/2.64,gamma2=365/3.03,mu=1/75,nu1=0,nu2=0,
              rho=1,rho1=1,rho2=1,eta=0,eta2=0,tau=0,tau2=0,sigmaT=365/(30*2),
              theta1=bestfitsparms[idx,'theta1'],
              theta2=bestfitsparms[idx,'theta2'])
  parms0['beta0_1'] <- parms0['R0_1']*(parms0['gamma1']+parms0['mu'])
  parms0['beta0_2'] <- parms0['R0_2']*(parms0['gamma2']+parms0['mu'])
  
  # bivalent vaccine
  parms1 <- parms0
  parms1['tau'] <- 1
  parms1['tau2'] <- 1
  
  # H1 vaccine
  parms2 <- parms0
  parms2['tau2'] <- 1
  parms2['tau'] <- parms0['theta1']
  parms2['rho2'] <- 0
  parms2['sigmaV2'] <- 0
  
  # H3 vaccine
  parms3 <- parms0
  parms3['tau'] <- 1
  parms3['tau2'] <- parms0['theta2']
  parms3['rho1'] <- 0
  parms3['sigmaV1'] <- 0
  
  start <- c(S=0.998,I1=0.001,I2=0.001,V=0,R12=0,R1=0,R2=0,I12=0,I21=0,IV2=0,
             IV12=0,RT1=0,RT2=0,RV12=0,RV1=0,RV2=0,RVT2=0,IV1=0,RVT1=0,IV21=0)
  times <- seq(0,100,by=1/104)
  
  # vaccine breadth
  rel.parm <- c()
  rel.parm['rel.sigmaV']<-1
  stats.sum <- data.frame()
  stats.sum2 <- data.frame()
  simdf.meta.list<-list()
  rel.rho<-seq(0,1,by=0.1)
  for (scaling in c(0.5,1,2)){
    rel.parm['rel.tau']<- scaling
    
    simdf.list1 <- lapply(rel.rho,rel.parm=rel.parm,parms.init=parms1,odefunc=sirmod,
                          start.cond=start,times.vec=times,sim)
    sumry.stats1 <- lapply(simdf.list1,sumry.cycl) %>% bind_rows()
    sumry.stats1 <- sumry.stats1 %>%
      mutate(scaling=scaling,vaccineImmunity='Bivalent vaccine')

    # H1 vaccine
    simdf.list2 <- lapply(rel.rho,rel.parm=rel.parm,parms.init=parms2,odefunc=sirmod,
                          start.cond=start,times.vec=times,sim)
    sumry.stats2 <- lapply(simdf.list2,sumry.cycl) %>% bind_rows()
    sumry.stats2 <- sumry.stats2 %>%
      mutate(scaling=scaling,vaccineImmunity='H1 vaccine')
    # H3 vaccine
    simdf.list3 <- lapply(rel.rho,rel.parm=rel.parm,parms.init=parms3,odefunc=sirmod,
                          start.cond=start,times.vec=times,sim)
    sumry.stats3 <- lapply(simdf.list3,sumry.cycl) %>% bind_rows()
    sumry.stats3 <- sumry.stats3 %>%
      mutate(scaling=scaling,vaccineImmunity='H3 vaccine')
    
    stats.sum <- bind_rows(stats.sum,sumry.stats1,sumry.stats2,sumry.stats3)
    tmp.list <- list(simdf.list1,scaling=scaling)
    
    simdf.meta.list <- append(simdf.meta.list,tmp.list)
  }
  stats.sum$scaling2 <- paste0('b==',stats.sum$scaling)
  
  # H3 vaccine: vaccine-immunity duration vs natural immunity duration
  rel.parm['rel.tau']<- 1
  stats.sum.vac.dur <- data.frame()
  simdf.meta.list.vac.dur<-list()
  rel.rho<-seq(0,1,by=0.1)
  for (scaling in c(0.5,1,2)){
    rel.parm['rel.sigmaV']<-scaling
    simdf.list1 <- lapply(rel.rho,rel.parm=rel.parm,parms.init=parms1,odefunc=sirmod,
                          start.cond=start,times.vec=times,sim)
    sumry.stats1 <- lapply(simdf.list1,sumry.cycl) %>% bind_rows()
    sumry.stats1 <- sumry.stats1 %>%
      mutate(scaling=scaling,
             vaccineImmunity='Bivalent vaccine')
    # H1 vaccine
    simdf.list2 <- lapply(rel.rho,rel.parm=rel.parm,parms.init=parms2,odefunc=sirmod,
                          start.cond=start,times.vec=times,sim)
    sumry.stats2 <- lapply(simdf.list2,sumry.cycl) %>% bind_rows()
    sumry.stats2 <- sumry.stats2 %>%
      mutate(scaling=scaling,vaccineImmunity='H1 vaccine')
    # H3 vaccine
    simdf.list3 <- lapply(rel.rho,rel.parm=rel.parm,parms.init=parms3,odefunc=sirmod,
                          start.cond=start,times.vec=times,sim)
    sumry.stats3 <- lapply(simdf.list3,sumry.cycl) %>% bind_rows()
    sumry.stats3 <- sumry.stats3 %>%
      mutate(scaling=scaling,vaccineImmunity='H3 vaccine')
    
    stats.sum.vac.dur <- bind_rows(stats.sum.vac.dur,sumry.stats1,sumry.stats2,
                                   sumry.stats3)
    
    tmp.list <- list(simdf.list1,scaling=scaling)
    
    simdf.meta.list.vac.dur <- append(simdf.meta.list.vac.dur,tmp.list)
  }
  
  stats.sum.vac.dur <- stats.sum.vac.dur %>%
    mutate(scaling2=factor(paste0('b==',scaling)))
  
  save(simdf.meta.list,stats.sum,simdf.meta.list.vac.dur,
       stats.sum.vac.dur,file=paste0('sim_scale2Clades_',idx,'.RData'))
}

