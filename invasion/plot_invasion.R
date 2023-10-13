rm(list=ls())
load('H3vac.invad.ana_noNatImm.RData')
library(RColorBrewer)
library(dplyr)
library(ggplot2)
simdf <- data.frame()
for(k in 1:length(invad.list)){
  df <- as.data.frame(invad.list[[k]])
  df['rho'] <- vac.parms.list[[k]][3]
  df['sigmaV'] <- vac.parms.list[[k]][2]
  df['tau2'] <- vac.parms.list[[k]][1]
  simdf <- simdf %>%
    bind_rows(df)
}

# first 4 weeks after invasion
simdf.t0 <- simdf %>%
  filter(time==0)
simdf.t1 <- simdf %>%
  filter(time==1/(52*7))
simdf.t2 <- simdf %>%
  filter(time==2/(52*7))
simdf.t3 <- simdf %>%
  filter(time==3/(52*7))

simdf.t3$growth.rate1 <- (simdf.t3$inci1-simdf.t2$inci1)/simdf.t2$inci1
simdf.t3$growth.rate <- (simdf.t2$inci1-simdf.t1$inci1)/simdf.t1$inci1

simdf.t3.rho <- simdf.t3 %>%
  # filter(rho==0) %>%
  mutate(subpl="vaccine immunity duration (yr)", 
         pl.y=1/sigmaV, grp=signif(rho/52*100,3)) %>%
  filter(grp==0 | grp==0.577 | grp==1.15 | grp==1.73)

simdf.t3.sigmaV <- simdf.t3 %>%
  # filter(sigmaV==2) %>%
  mutate(subpl="vaccination rate", 
         pl.y=rho/52*100, grp=1/sigmaV) %>%
  filter(grp==1 | grp==4 | grp==16 | grp==64)

pl.sigmaV <- ggplot(simdf.t3.sigmaV,aes(x=tau2,y=pl.y,z=growth.rate1)) +
  facet_wrap(~grp,scales='free') + 
  theme_classic()+ theme(strip.background = element_blank(),
                         strip.placement = "outside")+
  # geom_contour_filled(breaks=c(-100,0,1,2))+
  geom_contour_filled(bins=6,breaks = seq(-0.4,0.2,by=0.1))+
  labs(fill='Invading strain\ninitial growth rate')+
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), 
                     labels = c('0','0.25','0.5','0.75','1'))+
  scale_y_continuous(breaks = c(0,0.5,1,1.5,2),label=
                       c('0','0.5','1','1.5','2'))+
  # scale_fill_manual(values=c('#1a9850','#ffffcc','#feb24c','#fc4e2a','#800026'),
  #                   labels=c('<0','0-1','1-2','>2'))+
  scale_fill_manual(values=c('#1a9850','#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c'))+
                             # '#fc4e2a','#e31a1c','#bd0026','#800026'))+
  labs(x=expression('Susceptibility reduction by vaccines'~tau),
       y=expression('Vaccination rate (%) per week'~rho))
ggsave('invasion/H3vac.H1invad.rho.png',pl.sigmaV, width=6,height=4)

pl.rho <- ggplot(simdf.t3.rho,aes(x=tau2,y=pl.y,z=growth.rate1)) +
  facet_wrap(~grp,scales='free') + 
  theme_classic()+ theme(strip.background = element_blank(),
                         strip.placement = "outside")+
  # geom_contour_filled(breaks = c(-100,0,1,2)) +
  geom_contour_filled(bins=6,breaks = seq(-0.4,0.2,by=0.1))+
  labs(fill='Invading strain\ninitial growth rate') +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), 
                     labels = c('0','0.25','0.5','0.75','1'))+
  scale_y_continuous(expand = c(0,0))+
  # scale_fill_manual(values=c('#1a9850','#ffffcc','#feb24c','#fc4e2a','#800026'),
  #                   labels=c('<0','0-1','1-2','>2'))+
  scale_fill_manual(values=c('#1a9850','#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c'))+
  # ,'#fc4e2a','#e31a1c','#bd0026','#800026'))+
  labs(x=expression('Susceptibility reduction by vaccines'~tau),
       y='Vaccine immunity duration (yr)')

ggsave('invasion/H3vac.H1invad.vacdur.png',pl.rho, width=5,height=3)


rm(list=ls())
load('H1vac.invad.ana.RData')
library(RColorBrewer)
library(dplyr)
library(ggplot2)
simdf <- data.frame()
for(k in 1:length(invad.list)){
  df <- as.data.frame(invad.list[[k]])
  df['rho'] <- vac.parms.list[[k]][3]
  df['sigmaV'] <- vac.parms.list[[k]][2]
  df['tau'] <- vac.parms.list[[k]][1]
  # df <- df %>%
  #   filter(time<50)
  simdf <- simdf %>%
    bind_rows(df)
}

# first 4 weeks after invasion
simdf.t0 <- simdf %>%
  filter(time==0)
simdf.t1 <- simdf %>%
  filter(time==1/(52*7))
simdf.t2 <- simdf %>%
  filter(time==2/(52*7))
simdf.t3 <- simdf %>%
  filter(time==3/(52*7))

simdf.t3$growth.rate2 <- (simdf.t3$inci2-simdf.t2$inci2)/simdf.t2$inci2
max(simdf.t3$growth.rate2)
min(simdf.t3$growth.rate2)
simdf.t3$growth.rate <- (simdf.t2$inci2-simdf.t1$inci2)/simdf.t1$inci2
max(simdf.t3$growth.rate)
min(simdf.t3$growth.rate)

simdf.t3.rho <- simdf.t3 %>%
  # filter(rho==0) %>%
  mutate(subpl="vaccine immunity duration (yr)", 
         pl.y=1/sigmaV, grp=signif(rho/52*100,3)) 
# %>%
  # filter(grp==0 | grp==0.577 | grp==1.15 | grp==1.73)

simdf.t3.sigmaV <- simdf.t3 %>%
  # filter(sigmaV==2) %>%
  mutate(subpl="vaccination rate", 
         pl.y=rho/52*100, grp=1/sigmaV)
# max(simdf.t3.rho$growth.rate2)
pl.sigmaV <- ggplot(simdf.t3.sigmaV,aes(x=tau,y=pl.y,z=growth.rate2)) +
  facet_wrap(~grp,scales='free') + 
  theme_classic()+ theme(strip.background = element_blank(),
                         strip.placement = "outside")+
  # geom_contour_filled(breaks=c(-100,0,1,2))+
  # geom_contour_filled(bins=10,breaks = seq(-1.2,1.8,by=0.3))+
  geom_contour_filled(bins=6,breaks = seq(-0.3,0.3,by=0.1))+
  labs(fill='Invading strain\ninitial growth rate')+
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), 
                     labels = c('0','0.25','0.5','0.75','1'))+
  scale_y_continuous(breaks = c(0,0.5,1,1.5,2),label=
                       c('0','0.5','1','1.5','2'))+
  # scale_fill_manual(values=c('#1a9850','#ffffcc','#feb24c','#fc4e2a','#800026'),
  #                   labels=c('<0','0-1','1-2','>2'))+
  scale_fill_manual(values=c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a'))+
                             # ,'#e31a1c','#bd0026','#800026'))+
  labs(x=expression('Susceptibility reduction by vaccines'~tau),
       y=expression('Vaccination rate (%) per week'~rho))
ggsave('invasion/H1vac.H3invad.rho.png',pl.sigmaV, width=6,height=4)

pl.rho <- ggplot(simdf.t3.rho,aes(x=tau,y=pl.y,z=growth.rate2)) +
  facet_wrap(~grp,scales='free') + 
  theme_classic()+ theme(strip.background = element_blank(),
                         strip.placement = "outside")+
  geom_contour_filled(bins=6,breaks = seq(-0.3,0.3,by=0.1))+
  labs(fill='Invading strain\ninitial growth rate') +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), 
                     labels = c('0','0.25','0.5','0.75','1'))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values=c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c',
                             '#fc4e2a'))+
                             # '#e31a1c','#bd0026','#800026'))+
  labs(x=expression('Susceptibility reduction by vaccines'~tau),
       y='Vaccine immunity duration (yr)')

ggsave('invasion/H1vac.H3invad.vacdur.png',pl.rho, width=6,height=4)

