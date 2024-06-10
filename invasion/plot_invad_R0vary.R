rm(list=ls())
load('H3endemicVac.pandStrInvad.NoNatImm.RData')
library(ggplot2)
library(dplyr)
library(plyr)

simdf.list <- lapply(invad.list.R0, function(x){
  invad.list <- lapply(x, function(df){
    a <- ldply(df,rbind)
    return(a)
  })
  return(invad.list)
})

simdf_all <-ldply(simdf.list, function(x){
  simdf <- ldply(x,rbind)
})

simdf_all <- simdf_all %>%
  mutate(H1log=log(I1+I21+IV1))

simdf.t0 <- simdf_all %>%
  filter(time==0)
simdf.t1 <- simdf_all %>%
  filter(time==1/(52*7))
simdf.t2 <- simdf_all %>%
  filter(time==2/(52*7))
simdf.t3 <- simdf_all %>%
  filter(time==3/(52*7))
simdf.t4 <- simdf_all %>%
  filter(time==4/(52*7))

simdf.t3$growth.rate.t3 <- (simdf.t3$H1log-simdf.t2$H1log)*52*7
simdf.t4$growth.rate.t4 <- (simdf.t4$H1log-simdf.t3$H1log)*52*7
simdf.t3.sigmaV <- simdf.t3 %>%
  mutate(subpl="vaccination rate", 
         pl.y=H3vac.tau2, grp=1/H3vac.sigmaV)

pl.sigmaV.R0 <- ggplot(simdf.t3.sigmaV,aes(x=H3vac.cov/52*100,y=pl.y,z=growth.rate.t3/52/7)) +
  facet_grid(grp~R0_1,scales='free',labeller = label_bquote(rows=sigma[V]:.(1/grp), 
                                                            cols=R[0]^1:.(R0_1)))+
  theme_classic()+
  theme(panel.spacing = unit(0.8, "lines"))+
  geom_contour_filled(bins=6,breaks = c(-150,-75,0,75,150,225,300)/52/7)+
  geom_contour(bins=6,breaks = c(-150,-75,0,75,150,225,300)/52/7,
               aes(colour = factor(..level..==0,levels = c(T, F), 
                                   labels = c("Pandemic emergence boundary","")))) + 
  scale_colour_manual(values = c("darkred", "#00000000")) + 
  labs(fill='Invading strain initial\ngrowth rate (per day)',color='')+
  xlab('Vaccination rate (%) per week')+
  labs(y=expression('susceptibility reduction by vaccines'~tau))+
  # scale_x_continuous(expand = c(0,0),labels=c('0','0.5','1','1.5','2'))+
  scale_y_continuous(expand = c(0,0),labels=c('0','0.25','0.5','0.75','1'))+
  # scale_fill_manual(values=rev(heat.colors(12)))+
  scale_fill_manual(values=c('#287A22FF','#8DBC80FF',
                             '#fed976','#feb24c','#fc4e2a','#bd0026'))
ggsave('figures/H3endemicVac.pandStrR0varyInvad.NoNatImm.png',pl.sigmaV.R0, width=9,height=7)

# simdf.t3.rho <- simdf.t3 %>%
#   # filter(vac.cov==unique(simdf.t3$vac.cov)[2]) %>%
#   mutate(subpl="vaccine immunity duration (yr)", 
#          pl.y=1/H3vac.sigmaV, grp=signif(H3vac.cov/52*100,3)) %>%
#   # filter(R0==1|R0==1.4|R0==1.8|R0==2.2|R0==2.6|R0==3) %>%
#   # filter(grp==0 | grp==0.577 | grp==1.15 | grp==1.73)
#   filter(grp==1.920)
# 
# simdf.t3.sigmaV <- simdf.t3 %>%
#   # filter(vac.sigmaV==unique(simdf.t3$vac.sigmaV)[3]) %>%
#   mutate(subpl="vaccination rate (%) per week", 
#          pl.y=H3vac.cov/52*100, grp=1/H3vac.sigmaV) %>%
#   # filter(grp==1 | grp==4 | grp==16 | grp==64)
#   filter(grp==16)

# simdf.pl <- simdf.t3.rho %>%
#   rbind(simdf.t3.sigmaV)

# simdf.pl.R0small <- simdf.pl %>%
#   filter(R0_1==1.3)
# simdf.pl.R0big <- simdf.pl %>%
#   filter(R0_1==3.2)

# library(paletteer)
# paletteer_dynamic("cartography::green.pal", 5)
# 
# pl <- ggplot(simdf.pl.R0big,aes(x=H3vac.tau2,y=pl.y,z=growth.rate.t3)) +
#   facet_wrap(subpl~.,scales='free',nrow=2,strip.position = "left")+
#   theme_classic()+theme(strip.placement="outside",strip.background = element_blank(),
#                         strip.text.y = element_text(size = 10, angle = -90,colour = 'black'))+
#   geom_contour_filled(bins=9,breaks = c(-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.4,0.8,1.6))+
#   geom_contour(bins=9,breaks = c(-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.4,0.8,1.6),
#                aes(colour = factor(..level..==0,levels = c(T, F), 
#                                   labels = c("Pandemic emergence","")))) + 
#   scale_colour_manual(values = c("darkred", "#00000000")) + 
#   labs(fill='Invading strain initial\ngrowth rate',color='')+
#   ylab(NULL)+
#   scale_x_continuous(expand = c(0,0),labels=c('0','0.25','0.5','0.75','1'))+
#   scale_y_continuous(expand = c(0,0))+
#   # scale_fill_manual(values=rev(heat.colors(12)))+
#   scale_fill_manual(values=c('#5D9D52FF','#8DBC80FF','#B8D9A9FF',
#                              '#fed976','#feb24c',
#                              '#fc4e2a','#e31a1c','#bd0026','#800026'))+
#                              # '#CCBAD7FF', '#BA9FC7FF','#A783B6FF','#8F5EA1FF',
#                              # '#773A8BFF','#572872FF','#361959FF','black'))+
#                              # '#e31a1c'
#                              # ,'#bd0026','#800026'))+
#   # scale_color_manual(values=c('#1a9850','#ffffcc','#ffeda0','#fed976','#feb24c',
#   #                             'black','#fc4e2a','#e31a1c','#bd0026','#800026'))+
#   labs(x=expression('susceptibility reduction by vaccines'~tau))
# pl
# ggsave('H3endemicVac.pandStrR0_3.2_Invad.NoNatImm.png',pl, width=6,height=6)
# 
# pl2 <- ggplot(simdf.pl.R0small,aes(x=H3vac.tau2,y=pl.y,z=growth.rate.t3)) +
#   facet_wrap(subpl~.,scales='free',nrow=2,strip.position = "left")+
#   theme_classic()+theme(strip.placement="outside",strip.background = element_blank(),
#                         strip.text.y = element_text(size = 10, angle = -90,colour = 'black'))+
#   geom_contour_filled(bins=9,breaks = c(-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.4,0.8,1.6))+
#   geom_contour(bins=9,breaks = c(-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.4,0.8,1.6),
#                aes(colour = factor(..level..==0,levels = c(T, F), 
#                                    labels = c("Pandemic emergence","")))) + 
#   scale_colour_manual(values = c("darkred", "#00000000")) + 
#   labs(fill='Invading strain initial\ngrowth rate',color='')+
#   ylab(NULL)+
#   scale_x_continuous(expand = c(0,0),labels=c('0','0.25','0.5','0.75','1'))+
#   scale_y_continuous(expand = c(0,0))+
#   # scale_fill_manual(values=rev(heat.colors(12)))+
#   scale_fill_manual(values=c('#287A22FF','#5D9D52FF','#8DBC80FF','#B8D9A9FF',
#                              '#fed976','#feb24c',
#                              '#fc4e2a','#e31a1c','#bd0026','#800026'))+
#   # '#CCBAD7FF', '#BA9FC7FF','#A783B6FF','#8F5EA1FF',
#   # '#773A8BFF','#572872FF','#361959FF','black'))+
#   # '#e31a1c'
#   # ,'#bd0026','#800026'))+
#   # scale_color_manual(values=c('#1a9850','#ffffcc','#ffeda0','#fed976','#feb24c',
#   #                             'black','#fc4e2a','#e31a1c','#bd0026','#800026'))+
#   labs(x=expression('susceptibility reduction by vaccines'~tau))
# 
# pl2
# ggsave('H3endemicVac.pandStrR0_1.3_Invad.NoNatImm.png',pl, width=6,height=6)


# simdf.t3.rho <- simdf.t3 %>%
#   mutate(subpl="vaccine immunity duration (yr)", 
#          pl.y=1/H3vac.sigmaV, grp=signif(H3vac.cov/52*100,3)) %>%
#   filter(grp==0 | grp==0.577 | grp==1.15 | grp==1.73)
# 
# simdf.t3.rho <- simdf.t3 %>%
#   mutate(subpl="vaccine immunity duration (yr)", 
#          pl.y=1/H3vac.sigmaV, grp=signif(H3vac.cov/52*100,3)) %>%
#   filter(H3vac.tau2==0 | H3vac.tau2==0.25 |H3vac.tau2==0.5 | H3vac.tau2==0.75 | 
#            H3vac.tau2==1)
# ggplot(simdf.t3.rho,aes(x=H3vac.cov/52*100,y=pl.y,z=growth.rate.t3*52)) +
#   facet_grid(H3vac.tau2~R0_1,scales='free',labeller = label_bquote(rows=tau:.(H3vac.tau2), 
#                                                             cols=R[0]^1:.(R0_1)))+
#   theme_classic()+
#   geom_contour_filled(bins=6,breaks = 52*c(-0.1,0,0.1,0.2,0.4,0.8,1.6))+
#   geom_contour(bins=6,breaks = 52*c(-0.1,0,0.1,0.2,0.4,0.8,1.6),
#                aes(colour = factor(..level..==0,levels = c(T, F), 
#                                    labels = c("Pandemic emergence","")))) + 
#   scale_colour_manual(values = c("darkred", "#00000000")) + 
#   labs(fill='Invading strain initial\ngrowth rate (per week)',color='')+
#   ylab('Vaccine immunity duration (yr)')+
#   # scale_x_continuous(expand = c(0,0),labels=c('0','0.25','0.5','0.75','1'))+
#   scale_y_continuous(expand = c(0,0))+
#   # scale_fill_manual(values=rev(heat.colors(12)))+
#   scale_fill_manual(values=c('#5D9D52FF',
#                              '#fed976','#feb24c','#fc4e2a','#e31a1c','#bd0026'))
# 
# pl.rho.R0 <- ggplot(simdf.t3.rho,aes(x=H3vac.tau2,y=pl.y,z=growth.rate.t3*52)) +
#   facet_grid(grp~R0_1,scales='free',labeller = label_bquote(rows=rho:.(grp), 
#                                                           cols=R[0]^1:.(R0_1)))+
#   theme_classic()+
#   geom_contour_filled(bins=6,breaks = 52*c(-0.1,0,0.1,0.2,0.4,0.8,1.6))+
#   geom_contour(bins=6,breaks = 52*c(-0.1,0,0.1,0.2,0.4,0.8,1.6),
#                aes(colour = factor(..level..==0,levels = c(T, F), 
#                                    labels = c("Pandemic emergence","")))) + 
#   scale_colour_manual(values = c("darkred", "#00000000")) + 
#   labs(fill='Invading strain initial\ngrowth rate (per week)',color='')+
#   ylab('Vaccine immunity duration (yr)')+
#   scale_x_continuous(expand = c(0,0),labels=c('0','0.25','0.5','0.75','1'))+
#   scale_y_continuous(expand = c(0,0))+
#   # scale_fill_manual(values=rev(heat.colors(12)))+
#   scale_fill_manual(values=c('#5D9D52FF',
#                              '#fed976','#feb24c','#fc4e2a','#e31a1c','#bd0026'))+
#                              # ,'#800026'))+
#   # '#CCBAD7FF', '#BA9FC7FF','#A783B6FF','#8F5EA1FF',
#   # '#773A8BFF','#572872FF','#361959FF','black'))+
#   # '#e31a1c'
#   # ,'#bd0026','#800026'))+
#   # scale_color_manual(values=c('#1a9850','#ffffcc','#ffeda0','#fed976','#feb24c',
#   #                             'black','#fc4e2a','#e31a1c','#bd0026','#800026'))+
#   labs(x=expression('susceptibility reduction by vaccines'~tau))
# 
# pl.rho.R0
# ggsave('H3endemicVac.pandStrR0varyRhoInvad.NoNatImm.png',pl.rho.R0, width=9,height=7)

# trajectory
## trough / 5 years trajectory
rm(list=ls())
library(deSolve)
load('H3Vac.invad.RData')
times.invad <- seq(0,5,by=1/(52*7))

simdf.strongVac <- simdf.list[[958]]
simdf.weakVac <- simdf.list[[886]]
simdf.midVac <- simdf.list[[926]]
simdf.list.sh <- list(simdf.weakVac,simdf.midVac,simdf.strongVac)
invad.list.R0.sh <- lapply(c(1,3.2),function (r){
  parms.invad <- parms3
  # natural immunity aginst seasonal variant
  # parms.invad['theta2']=0.5 
  # pandemic strain
  parms.invad['theta2']=0
  parms.invad['R0_1'] <- r
  parms.invad['beta0_1']<-r*(parms.invad['gamma2']+
                               parms.invad['mu'])
  parms.invad['sigma1'] <- parms.invad['sigma2']
  invad.list.sh <- lapply(simdf.list.sh, function(simdf) {
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
  return(list(invad.list.sh,r))
})

simdf.list.sh <- lapply(invad.list.R0.sh, function(x){
  invad.list <- x[[1]]
  trough <- ldply(invad.list,rbind)
  trough$R0 <- x[[2]]
  return (trough)
})

trough <- ldply(simdf.list.sh,rbind)

trough_smallR0 <- simdf.list.sh[[1]]
trough_midR0 <- simdf.list.sh[[2]]
trough_bigR0 <- simdf.list.sh[[3]]

ggplot(data=trough_smallR0,mapping = aes(x=time))+
  geom_line(aes(y=log(inci1),colour=factor(vac.tau2)))

ggplot(data=trough_midR0,mapping = aes(x=time))+
  geom_line(aes(y=log(inci1),colour=factor(vac.tau2)))

ggplot(data=trough_bigR0,mapping = aes(x=time))+
  geom_line(aes(y=log(inci1),colour=factor(vac.tau2)))

simdf.rho <- simdf %>%
  mutate(subpl="vaccine immunity duration (yr)", 
         pl.y=1/vac.sigmaV, grp=signif(vac.cov/52*100,3)) %>%
  filter(R0==1) %>%
  filter(grp==1.920)

test <- simdf.rho %>%
  filter(vac.tau2==0 & pl.y==4)

ggplot(data=test,aes(x=time,y=inci1))+
  geom_point()

trough.rho <- trough %>%
  filter(vac.cov==1) %>%
  mutate(subpl="vaccine immunity duration (yr)",pl.y=vac.eta)

trough.sigmaV <- trough %>%
  filter(vac.sigmaV==0) %>%
  mutate(subpl="vaccination rate",pl.y=vac.cov)

trough.pl <- trough.rho %>%
  bind_rows(trough.sigmaV)

pl.trough <- ggplot(trough.pl,aes(x=vac.tau,y=pl.y,z=min.inci2)) +
  facet_grid(subpl~R0,scales='free',labeller = labeller(R0=label_both))+
  geom_contour_filled(color='black',bins=10)+
  theme_classic()+
  labs(fill='Depth of trough')+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_brewer(palette="RdGy",direction=-1)+
  labs(x=expression('susceptibility reduction by vaccines'~tau))

ggsave(file='trough.png',plot=pl.trough,width=14,height=4)