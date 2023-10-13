library(deSolve)
# load('H3Vac.invad.short.RData')
load('H3Vac.invad.RData')
times.invad <- seq(0,5/52,by=1/(52*7))
invad.list.R0 <- lapply(as.array(seq(1,4,by=0.2)),function (r){
  parms.invad <- parms3
  # natural immunity
  parms.invad['theta2']=0.5
  # parms.invad['theta2']=0
  parms.invad['R0_1'] <- r
  parms.invad['beta0_1']<-r*(parms.invad['gamma1']+
                               parms.invad['mu'])
  invad.list <- lapply(simdf.list, function(simdf) {
    start.invad <- unlist(simdf[nrow(simdf),2:23])
    start.invad['I1'] <- 1e-6
    start.invad['S'] <- start.invad['S']-1e-6
    parms.invad['rho']=simdf[1,'vac.cov']
    parms.invad['tau2']=simdf[1,'vac.tau2']
    parms.invad['sigmaV']=simdf[1,'vac.sigmaV']
    invad.output <- as.data.frame(rk(y=start.invad, times=times.invad,
                       func=sirmod,parms=parms.invad,hmax=1/365))
    invad.output['vac.cov']<-simdf[1,'vac.cov']
    invad.output['vac.tau2']<-simdf[1,'vac.tau2']
    invad.output['vac.sigmaV']<-simdf[1,'vac.sigmaV']
    return(invad.output)
  })
  return(list(invad.list,r))
  })

# save.image('H3vac.H1invad.short.R0.results.RData')
save.image('H3vac.H1invad.R0vary_strongNatImm.RData')
# save.image('H3vac.H1invad.R0vary_noNatImm.RData')