rm(list=ls())
# setwd('invasion/')
library(deSolve)
## 1. Simulate the situation where H3 is endemic
## in the presence of the vaccine (full immunity against H3) and H1 is absent.

## immunization
rho<-seq(from=0,to=1,by=0.1) # vaccine coverage
# cross protection by vaccination
# tau2<-seq(0,1,by=0.1)
sigmaV<-1/c(0.5,1,2,4,8,16)

# Updated Aug 2023
parms0 <- c(R0_2=1.60,a2=0.04,sigma2=1/2.7,sigmaV=1/2.7,gamma2=365/3.03,
            mu=1/75,rho=1,tau=1,sigmaT=365/(30*2))
parms0['beta0_2'] <- parms0['R0_2']*(parms0['gamma2']+parms0['mu'])
# H3 vaccine

start <- c(S=0.999,I2=0.001,V=0,R2=0,IV2=0,RT2=0,RV2=0,RVT2=0)
times <- seq(0,600,by=1/104)

# Update Aug 2023
sirmod <- function(t,y,parms){
  with(as.list(c(y, parms)), {
    beta0_2 <- R0_2*(gamma2+mu)
    beta2 <- beta0_2*(1+a2*cos(2*pi*t)) 
    dS <- mu-(rho+mu)*S-beta2*S*(I2+IV2)+sigma2*R2+sigmaV*V 
    dI2 <- beta2*S*(I2+IV2)-(gamma2+mu)*I2
    dRT2 <- gamma2*I2 - (mu+sigmaT)*RT2
    dR2 <- sigmaT*RT2 -(mu+sigma2)*R2
    dV <- rho*S+sigma2*RV2 -(mu+sigmaV)*V-beta2*(1-tau)*V*(I2+IV2)
    dIV2 <- beta2*(1-tau)*V*(I2+IV2)-(gamma2+mu)*IV2
    dRVT2 <- gamma2*IV2 - (sigmaT+mu)*RVT2
    dRV2 <- sigmaT*RVT2 -(sigma2+mu)*RV2
    res <- c(dS,dI2,dV,dR2,dIV2,dRT2,dRV2,dRVT2)
    
    return(list(res,inci2=beta2*S*(I2+IV2)+beta2*(1-tau)*V*(I2+IV2)))
  })
}
sim2 <- function(vac.parms,odefunc=sirmod,start.cond=start,
                 times.vec=times,parms.init=parms0){
  
  parms.init['rho'] <- vac.parms[1]
  parms.init['sigmaV'] <- vac.parms[2]
  
  SIR.output <- rk(y=start, times=times,func=odefunc,parms=parms.init,
                   hmax=1/365)
  
  sir.df <- as.data.frame(SIR.output)
  sir.df['H3vac.cov'] <- vac.parms[1]
  sir.df['H3vac.sigmaV'] <- vac.parms[2]
  
  res <- sir.df[sir.df["time"]>max(sir.df['time'])-50,]
  return(res)
}

# tau.m <- array(t(array(tau2,c(11,6))),c(11,11,6))
# rho.m <- array(unlist(lapply(rho,function(ary) array(ary,c(11,6)))),c(11,11,6))
# sigmaV.m <- array(matrix(array(sigmaV,c(6,6))),c(11,11,6))
# vac.parms.list <- mapply(function(x,y,z) list(c(x,y,z)),tau.m,sigmaV.m,rho.m)
rho.m <- t(array(rho,c(11,6)))
sigmaV.m <- array(sigmaV,c(6,11))
vac.parms.list <- mapply(function(x,y) list(c(x,y)),rho.m,sigmaV.m)

if(length(unique(vac.parms.list))!=length(vac.parms.list)) {
  print('Number of parameter settings are not correct. Please reset!')
} else{
  simdf.list <- lapply(vac.parms.list,sim2)
}
save.image('H3Vac.invad.RData')