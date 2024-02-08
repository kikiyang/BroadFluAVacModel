rm(list=ls())
library(dplyr)
library(deSolve)
library(lubridate)
library(ggplot2)
theme_set(theme_bw())

# load('USA.inci.RData')
inci <- read.csv('data/USA/no.cases_inci.allsubtypes.csv')
inci <- inci %>%
  select(YEAR,WEEK,month,inci.A.H1,inci.A.H3,cases.A.H1,cases.A.H3)

start.season <- 2004
end.season <- 2018
no.season <- c(2008,2009,2010)
## incidence = multiplying
## bunch of things
data <- inci %>%
  mutate(Season= if_else(month<7,as.integer(YEAR-1),YEAR)) %>%
  filter(Season<end.season+1 & Season>start.season-1) %>%
  filter(Season!=no.season[1] & Season!=no.season[2]&Season!=no.season[3])

stats.season <- data %>%
  group_by(Season) %>%
  summarise(inci.total1 = sum(inci.A.H1,na.rm=TRUE),
            inci.total2 = sum(inci.A.H3,na.rm=TRUE))

stats.season <- stats.season %>%
  mutate(predom=if_else(inci.total1>inci.total2,'H1',if_else(inci.total1==inci.total2,'H1 or H3','H3')))

epi.size1.data <- stats.season[['inci.total1']]
epi.size2.data <- stats.season[['inci.total2']]
sumry.stats.data <- data.frame()
sumry.stats.data[1,'corr'] <- cor(epi.size1.data,epi.size2.data)
sumry.stats.data[1,'auto_corr_h1'] <- acf(epi.size1.data, plot=FALSE)[[1]][2,,1]
sumry.stats.data[1,'auto_corr_h3'] <- acf(epi.size2.data, plot=FALSE)[[1]][2,,1]
sumry.stats.data[1,'coeff_var_h1'] <- sd(epi.size1.data)/mean(epi.size1.data)
sumry.stats.data[1,'coeff_var_h3'] <- sd(epi.size2.data)/mean(epi.size2.data)

sir.mod.fit <- function(t,y,parms){
  with(as.list(c(y, parms)), {
    beta0_1 <- R0_1*(gamma1+mu)
    beta0_2 <- R0_2*(gamma2+mu)
    
    beta1 <- beta0_1*(1+a1*cos(2*pi*t)) #seasonality forcing
    beta2 <- beta0_2*(1+a2*cos(2*pi*t)) #seasonality forcing
    
    dS <- mu-(rho+mu)*S-beta1*S*(I1+I21+IV1)-beta2*S*(I2+I12+IV2)+
      sigma1*R1+sigma2*R2+sigmaV*V #susceptible
    
    dI1 <- beta1*S*(I1+I21+IV1)-(gamma1+mu)*I1 #primary Strain1 infection
    dRT1 <- gamma1*I1+sigmaV*RVT1-(mu+sigmaT+rho2)*RT1 #strain-transcending immunized + infection-induced immunized against strain1
    dR1 <- sigma2*R12 +sigmaT*RT1 +sigmaV*RV1 -
      beta2*(1-theta1)*R1*(I2+I12+IV2) -(mu+sigma1+rho2)*R1 #infection-induced immunized against strain1
    dI12 <- beta2*(1-theta1)*R1*(I2+I12+IV2)-(mu+gamma2/(1-nu1))*I12 #recovered from primary strain1 infection and infected by strain2 
    
    dI2 <- beta2*S*(I2+I12+IV2)-(gamma2+mu)*I2 #primary Strain2 infection
    dRT2 <- gamma2*I2 + sigmaV*RVT2 - (mu+sigmaT+rho1)*RT2 #strain-transcending immunized + infection-induced immunized against strain2
    dR2 <- sigma1*R12 +sigmaT*RT2 +sigmaV*RV2 - 
      beta1*(1-theta2)*R2*(I1+I21+IV1) -(mu+sigma2+rho1)*R2 #infection-induced immunized against strain2
    dI21 <- beta1*(1-theta2)*R2*(I1+I21+IV1)-(mu+gamma1/(1-nu2))*I21 #recovered from primary strain2 infection and infected by strain1
    
    dV <- rho*S - (mu+sigmaV)*V + sigma2*RV2 + sigma1*RV1 - 
      beta1*(1-tau2)*V*(I1+I21+IV1) - beta2*(1-tau)*V*(I2+I12+IV2)
    
    dIV1 <- beta1*(1-tau2)*V*(I1+I21+IV1) - (gamma1/(1-eta2)+mu)*IV1
    dRVT1 <- gamma1/(1-eta2)*IV1 + rho2*RT1 - (sigmaV+sigmaT+mu)*RVT1
    dRV1 <- sigmaT*RVT1 + rho2*R1-(sigma1+sigmaV+mu)*RV1
    
    dIV2 <- beta2*(1-tau)*V*(I2+I12+IV2)-(gamma2/(1-eta)+mu)*IV2 #vaccinated individuals infected by strain2
    dRVT2 <- gamma2/(1-eta)*IV2 + rho1*RT2 - (sigmaV+sigmaT+mu)*RVT2 #strain-transcending immunized - recovered form strain2 infection
    dRV2 <- sigmaT*RVT2 +rho1*R2 -(sigma2+sigmaV+mu)*RV2 #vaccinated individuals recovered from strain2 infection 
    
    dR12 <- gamma2/(1-nu1)*I12 +gamma1/(1-nu2)*I21 -
      (mu+sigma1+sigma2)*R12 #recovered from natural infection of 1 and 2
    
    res <- c(dS,dI1,dI2,dV,dR12,dR1,dR2,dI12,dI21,dIV2,dRT1,dRT2,
             dRV1,dRV2,dRVT2,dIV1,dRVT1)
    
    return(list(res,inci1=beta1*S*(I1+I21+IV1)+beta1*(1-theta2)*R2*(I1+I21+IV1)+
                  beta1*(1-tau2)*V*(I1+I21+IV1),
                inci2=beta2*S*(I2+I12+IV2)+beta2*(1-theta1)*R1*(I2+I12+IV2)+
                  beta2*(1-tau)*V*(I2+I12+IV2)))
  })
}

# stats.all <- data.frame()
# for(k in 1:nrow(sumry.stats.df_a04_top)){
const_parms <- c(R0_1=1.44,R0_2=1.60,sigma1=1/3.12,sigma2=1/2.28,mu=1/75, #from Yang et al., 2020
                 # theta1=sumry.stats.df_a04_top[k,'theta1'],
                 # theta2=sumry.stats.df_a04_top[k,'theta2'],
                 theta1=0.55,theta2=0.35,
                 gamma1=365/2.64,gamma2=365/3.03, #from Yang et al., 2020
                 a1=0.04,a2=0.04,sigmaV=0,
                 # sigmaT=SIGMAT[idx],
                 sigmaT=365/(30*2), #from Ferguson et al., 2003
                 nu1=0,nu2=0,rho=0,rho1=0,rho2=0,eta=0,eta2=0,tau=0,tau2=0)

sim <- rk4(func=sir.mod.fit,
           times=seq(0,2000,by=1/104),
           y=c(S=0.998,I1=0.001,I2=0.001,V=0,R12=0,R1=0,R2=0,I12=0,I21=0,IV2=0,
               RT1=0,RT2=0,RV1=0,RV2=0,RVT2=0,IV1=0,RVT1=0),
           parms=c(const_parms['R0_1'],const_parms['R0_2'],
                   const_parms['sigma1'],const_parms['sigma2'],
                   const_parms['sigmaV'],const_parms['sigmaT'],
                   # theta1=parms2fit[1],theta2=parms2fit[2], 
                   const_parms['theta1'],const_parms['theta2'],
                   const_parms['a1'], const_parms['a2'],
                   const_parms['gamma1'],const_parms['gamma2'],
                   const_parms['mu'],const_parms['rho'],const_parms['rho1'],
                   const_parms['rho2'],const_parms['nu1'],const_parms['nu2'],
                   const_parms['tau'],const_parms['tau2'],
                   const_parms['eta'],const_parms['eta2']))
sim.output <- as.data.frame(sim) %>%
  filter(time>800)

sim.stats <- sim.output %>%
  mutate(year=year(date_decimal(time)),month=month(date_decimal(time))) %>%
  # filter(time>decimal_date(ymd("2004-06-30")) &
  #          time<decimal_date(ymd("2019-07-01"))) %>%
  mutate(Season=if_else(month<7,year-1,year)) %>%
  # filter(Season!=2009 & Season!=2010) %>%
  group_by(Season) %>%
  # summarise(inci.total1 = parms2fit[3]*sum(inci1,na.rm=TRUE),
  #           inci.total2 = parms2fit[3]*sum(inci2,na.rm=TRUE))%>%
  summarise(inci.total1 = sum(inci1,na.rm=TRUE),
            inci.total2 = sum(inci2,na.rm=TRUE))

sumry.stats.peryr <- data.frame()

for(i in 1:((nrow(sim.stats)-1)/12)){
  epi.size1.sim <- sim.stats[((i-1)*12+1):(i*12),]$inci.total1
  epi.size2.sim <- sim.stats[((i-1)*12+1):(i*12),]$inci.total2
  
  sumry.stats.peryr[i,'corr'] <- cor(epi.size1.sim,epi.size2.sim)
  sumry.stats.peryr[i,'auto_corr_h1'] <- acf(epi.size1.sim, plot=FALSE)[[1]][2,,1]
  sumry.stats.peryr[i,'auto_corr_h3'] <- acf(epi.size2.sim, plot=FALSE)[[1]][2,,1]
  sumry.stats.peryr[i,'coeff_var_h1'] <- sd(epi.size1.sim)/mean(epi.size1.sim)
  sumry.stats.peryr[i,'coeff_var_h3'] <- sd(epi.size2.sim)/mean(epi.size2.sim)
  sumry.stats.peryr[i,'theta1'] <- sumry.stats.df_a04_top[k,'theta1']
  sumry.stats.peryr[i,'theta2'] <- sumry.stats.df_a04_top[k,'theta2']
}
  
  # stats.all <- stats.all %>%
  #   bind_rows(sumry.stats.peryr)
# }
# for(idx in 1:40){
#   stats.all[(100*idx-99):(100*idx),'param_grp'] <- idx
# }

## plot results (horizental line to indicate data?)
# ggplot(stats.all,aes(corr)) +
#   facet_wrap(~param_grp) + 
ggplot(sumry.stats.peryr,aes(corr)) +
  geom_histogram() +
  geom_vline(xintercept=sumry.stats.data$corr,color='grey',linetype=2)+
  labs(x='Correlation of H1 and H3 incidence')

# ggplot(stats.all) + 
#   facet_wrap(~param_grp) + 
ggplot(sumry.stats.peryr) +
  geom_histogram(aes(auto_corr_h1),fill='darkred')+
  geom_histogram(aes(auto_corr_h3),fill='steelblue')+
  geom_vline(xintercept=sumry.stats.data$auto_corr_h1,color='darkred',linetype=2)+
  geom_vline(xintercept=sumry.stats.data$auto_corr_h3,color='steelblue',linetype=2)+
  labs(x='Auto-correlation of the incidence')
  
# ggplot(stats.all) + 
#   facet_wrap(~param_grp) +
ggplot(sumry.stats.peryr) +
  geom_histogram(aes(coeff_var_h1),fill='darkred')+
  geom_histogram(aes(coeff_var_h3),fill='steelblue')+
  geom_vline(xintercept=sumry.stats.data$coeff_var_h1,color='darkred',linetype=2)+
  geom_vline(xintercept=sumry.stats.data$coeff_var_h3,color='steelblue',linetype=2)+
  labs(x='Coefficient of variation of the incidence')


theta1 <- seq(0.05,1,by=0.05)
theta2 <- seq(0.05,1,by=0.05)
sumry.stats.ls <- list()
sumry.stats.df <- data.frame()
for(idx1 in 1:length(theta1)){
  for(idx2 in 1:length(theta2)){
    const_parms <- c(R0_1=1.44,R0_2=1.60,sigma1=1/3.12,sigma2=1/2.28,mu=1/75, #from Yang et al., 2020
                     theta1=theta1[idx1],theta2=theta2[idx2],gamma1=365/2.64,gamma2=365/3.03, #from Yang et al., 2020
                     a1=0.1,a2=0.1,sigmaV=0,
                     sigmaT=365/(30*2), #from Ferguson et al., 2003
                     nu1=0,nu2=0,rho=0,rho1=0,rho2=0,eta=0,eta2=0,
                     tau=0,tau2=0)
    sim <- rk4(func=sir.mod.fit,
               times=seq(0,1650,by=1/104),
               y=c(S=0.998,I1=0.001,I2=0.001,V=0,R12=0,R1=0,R2=0,I12=0,I21=0,IV2=0,
                   RT1=0,RT2=0,RV1=0,RV2=0,RVT2=0,IV1=0,RVT1=0),
               parms=c(const_parms['R0_1'],const_parms['R0_2'],
                       const_parms['sigma1'],const_parms['sigma2'],
                       const_parms['sigmaV'],const_parms['sigmaT'],
                       const_parms['theta1'],const_parms['theta2'],
                       const_parms['a1'], const_parms['a2'],
                       const_parms['gamma1'],const_parms['gamma2'],
                       const_parms['mu'],const_parms['rho'],const_parms['rho1'],
                       const_parms['rho2'],const_parms['nu1'],const_parms['nu2'],
                       const_parms['tau'],const_parms['tau2'],
                       const_parms['eta'],const_parms['eta2']))
    sim.output <- as.data.frame(sim) %>%
      filter(time>150)
    sim.stats <- sim.output %>%
      mutate(year=year(date_decimal(time)),month=month(date_decimal(time))) %>%
      mutate(Season=if_else(month<7,year-1,year)) %>%
      group_by(Season) %>%
      summarise(inci.total1 = sum(inci1,na.rm=TRUE),
                inci.total2 = sum(inci2,na.rm=TRUE))
    sumry.stats.per15yrs <- data.frame()
    for(i in 1:((nrow(sim.stats)-1)/15)){
      epi.size1.sim <- sim.stats[((i-1)*15+1):(i*15),]$inci.total1
      epi.size2.sim <- sim.stats[((i-1)*15+1):(i*15),]$inci.total2
      sumry.stats.per15yrs[i,'corr'] <- cor(epi.size1.sim,epi.size2.sim)
      sumry.stats.per15yrs[i,'auto_corr_h1'] <- acf(epi.size1.sim, plot=FALSE)[[1]][2,,1]
      sumry.stats.per15yrs[i,'auto_corr_h3'] <- acf(epi.size2.sim, plot=FALSE)[[1]][2,,1]
      sumry.stats.per15yrs[i,'coeff_var_h1'] <- sd(epi.size1.sim)/mean(epi.size1.sim)
      sumry.stats.per15yrs[i,'coeff_var_h3'] <- sd(epi.size2.sim)/mean(epi.size2.sim)
    }
    
    row.num = length(theta1)*(idx1-1)+idx2
    sumry.stats.df[row.num,'mean_corr'] <- mean(sumry.stats.per15yrs$corr)
    sumry.stats.df[row.num,'mean_auto_corr_h1'] <- mean(sumry.stats.per15yrs$auto_corr_h1)
    sumry.stats.df[row.num,'mean_auto_corr_h3'] <- mean(sumry.stats.per15yrs$auto_corr_h3)
    sumry.stats.df[row.num,'mean_coeff_var_h1'] <- mean(sumry.stats.per15yrs$coeff_var_h1)
    sumry.stats.df[row.num,'mean_coeff_var_h3'] <- mean(sumry.stats.per15yrs$coeff_var_h3)
    sumry.stats.df[row.num,'theta1'] <- theta1[idx1]
    sumry.stats.df[row.num,'theta2'] <- theta2[idx2]
    
    sumry.stats.ls[[row.num]] <- sumry.stats.per15yrs
  }
}

ggplot(data=sumry.stats.df,aes(theta1,theta2,fill=mean_corr))+
  geom_raster()+
  scale_fill_viridis_c()

ggplot(data=sumry.stats.df,aes(theta1,theta2,fill=mean_auto_corr_h1))+
  geom_raster()+
  scale_fill_viridis_c()

ggplot(data=sumry.stats.df,aes(theta1,theta2,fill=mean_auto_corr_h3))+
  geom_raster()+
  scale_fill_viridis_c()

ggplot(data=sumry.stats.df,aes(theta1,theta2,fill=mean_coeff_var_h1))+
  geom_raster()+
  scale_fill_viridis_c()

ggplot(data=sumry.stats.df,aes(theta1,theta2,fill=mean_coeff_var_h3))+
  geom_raster()+
  scale_fill_viridis_c()
load('sim_mean_dist_cos_a0.rdata')
for(i in 1:nrow(sumry.stats.df)){
  sumry.stats.df[i,'corr'] <- ifelse(is.na(sumry.stats.df[i,'mean_corr']),0,
                                           sumry.stats.df[i,'mean_corr'])
  sumry.stats.df[i,'ac1'] <- ifelse(is.na(sumry.stats.df[i,'mean_auto_corr_h1']),0,sumry.stats.df[i,'mean_auto_corr_h1'])
  sumry.stats.df[i,'ac3'] <- ifelse(is.na(sumry.stats.df[i,'mean_auto_corr_h3']),0,sumry.stats.df[i,'mean_auto_corr_h3'])
  
}

sumry.stats.df_a0 <- sumry.stats.df %>%
  mutate(a=0,dist=abs(corr+sumry.stats.data$corr)/abs(sumry.stats.data$corr)+
           abs(ac1-sumry.stats.data$auto_corr_h1)/abs(sumry.stats.data$auto_corr_h1)+
           abs(ac3-sumry.stats.data$auto_corr_h3)/abs(sumry.stats.data$auto_corr_h3)+
           abs(mean_coeff_var_h1-sumry.stats.data$coeff_var_h1)/abs(sumry.stats.data$coeff_var_h1)+
           abs(mean_coeff_var_h3-sumry.stats.data$coeff_var_h3)/abs(sumry.stats.data$coeff_var_h3)
  )

setwd('parameter/')


load('sim_mean_dist_cos_a002.rdata')
for(i in 1:nrow(sumry.stats.df)){
  sumry.stats.df[i,'corr'] <- ifelse(is.na(sumry.stats.df[i,'mean_corr']),0,
                                     sumry.stats.df[i,'mean_corr'])
  sumry.stats.df[i,'ac1'] <- ifelse(is.na(sumry.stats.df[i,'mean_auto_corr_h1']),0,sumry.stats.df[i,'mean_auto_corr_h1'])
  sumry.stats.df[i,'ac3'] <- ifelse(is.na(sumry.stats.df[i,'mean_auto_corr_h3']),0,sumry.stats.df[i,'mean_auto_corr_h3'])
  
}

sumry.stats.df_a02 <- sumry.stats.df %>%
  mutate(a=0.02,dist=abs(corr+sumry.stats.data$corr)/abs(sumry.stats.data$corr)+
           abs(ac1-sumry.stats.data$auto_corr_h1)/abs(sumry.stats.data$auto_corr_h1)+
           abs(ac3-sumry.stats.data$auto_corr_h3)/abs(sumry.stats.data$auto_corr_h3)+
           abs(mean_coeff_var_h1-sumry.stats.data$coeff_var_h1)/abs(sumry.stats.data$coeff_var_h1)+
           abs(mean_coeff_var_h3-sumry.stats.data$coeff_var_h3)/abs(sumry.stats.data$coeff_var_h3)
  )

load('sim_mean_dist_cos_a004.rdata')
for(i in 1:nrow(sumry.stats.df_a04)){
  sumry.stats.df[i,'corr'] <- ifelse(is.na(sumry.stats.df[i,'mean_corr']),0,
                                     sumry.stats.df[i,'mean_corr'])
  sumry.stats.df[i,'ac1'] <- ifelse(is.na(sumry.stats.df[i,'mean_auto_corr_h1']),0,sumry.stats.df[i,'mean_auto_corr_h1'])
  sumry.stats.df[i,'ac3'] <- ifelse(is.na(sumry.stats.df[i,'mean_auto_corr_h3']),0,sumry.stats.df[i,'mean_auto_corr_h3'])
  
}

sumry.stats.df_a04 <- sumry.stats.df_a04 %>%
  mutate(a=0.04,dist=abs(corr+sumry.stats.data$corr)/abs(sumry.stats.data$corr)+
           abs(ac1-sumry.stats.data$auto_corr_h1)/abs(sumry.stats.data$auto_corr_h1)+
           abs(ac3-sumry.stats.data$auto_corr_h3)/abs(sumry.stats.data$auto_corr_h3)+
           abs(mean_coeff_var_h1-sumry.stats.data$coeff_var_h1)/abs(sumry.stats.data$coeff_var_h1)+
           abs(mean_coeff_var_h3-sumry.stats.data$coeff_var_h3)/abs(sumry.stats.data$coeff_var_h3)
  )
load('sim_mean_dist_cos_a006.rdata')
for(i in 1:nrow(sumry.stats.df)){
  sumry.stats.df[i,'corr'] <- ifelse(is.na(sumry.stats.df[i,'mean_corr']),0,
                                     sumry.stats.df[i,'mean_corr'])
  sumry.stats.df[i,'ac1'] <- ifelse(is.na(sumry.stats.df[i,'mean_auto_corr_h1']),0,sumry.stats.df[i,'mean_auto_corr_h1'])
  sumry.stats.df[i,'ac3'] <- ifelse(is.na(sumry.stats.df[i,'mean_auto_corr_h3']),0,sumry.stats.df[i,'mean_auto_corr_h3'])
  
}

sumry.stats.df_a06 <- sumry.stats.df %>%
  mutate(a=0.06,dist=abs(corr+sumry.stats.data$corr)/abs(sumry.stats.data$corr)+
           abs(ac1-sumry.stats.data$auto_corr_h1)/abs(sumry.stats.data$auto_corr_h1)+
           abs(ac3-sumry.stats.data$auto_corr_h3)/abs(sumry.stats.data$auto_corr_h3)+
           abs(mean_coeff_var_h1-sumry.stats.data$coeff_var_h1)/abs(sumry.stats.data$coeff_var_h1)+
           abs(mean_coeff_var_h3-sumry.stats.data$coeff_var_h3)/abs(sumry.stats.data$coeff_var_h3)
  )
load('sim_mean_dist_cos_a008.rdata')
for(i in 1:nrow(sumry.stats.df)){
  sumry.stats.df[i,'corr'] <- ifelse(is.na(sumry.stats.df[i,'mean_corr']),0,
                                     sumry.stats.df[i,'mean_corr'])
  sumry.stats.df[i,'ac1'] <- ifelse(is.na(sumry.stats.df[i,'mean_auto_corr_h1']),0,sumry.stats.df[i,'mean_auto_corr_h1'])
  sumry.stats.df[i,'ac3'] <- ifelse(is.na(sumry.stats.df[i,'mean_auto_corr_h3']),0,sumry.stats.df[i,'mean_auto_corr_h3'])
  
}

sumry.stats.df_a08 <- sumry.stats.df %>%
  mutate(a=0.08,dist=abs(corr+sumry.stats.data$corr)/abs(sumry.stats.data$corr)+
           abs(ac1-sumry.stats.data$auto_corr_h1)/abs(sumry.stats.data$auto_corr_h1)+
           abs(ac3-sumry.stats.data$auto_corr_h3)/abs(sumry.stats.data$auto_corr_h3)+
           abs(mean_coeff_var_h1-sumry.stats.data$coeff_var_h1)/abs(sumry.stats.data$coeff_var_h1)+
           abs(mean_coeff_var_h3-sumry.stats.data$coeff_var_h3)/abs(sumry.stats.data$coeff_var_h3)
  )
load('sim_mean_dist_cos_a01.rdata')
for(i in 1:nrow(sumry.stats.df)){
  sumry.stats.df[i,'corr'] <- ifelse(is.na(sumry.stats.df[i,'mean_corr']),0,
                                     sumry.stats.df[i,'mean_corr'])
  sumry.stats.df[i,'ac1'] <- ifelse(is.na(sumry.stats.df[i,'mean_auto_corr_h1']),0,sumry.stats.df[i,'mean_auto_corr_h1'])
  sumry.stats.df[i,'ac3'] <- ifelse(is.na(sumry.stats.df[i,'mean_auto_corr_h3']),0,sumry.stats.df[i,'mean_auto_corr_h3'])
  
}

sumry.stats.df_a1 <- sumry.stats.df %>%
  mutate(a=0.1,dist=abs(corr+sumry.stats.data$corr)/abs(sumry.stats.data$corr)+
           abs(ac1-sumry.stats.data$auto_corr_h1)/abs(sumry.stats.data$auto_corr_h1)+
           abs(ac3-sumry.stats.data$auto_corr_h3)/abs(sumry.stats.data$auto_corr_h3)+
           abs(mean_coeff_var_h1-sumry.stats.data$coeff_var_h1)/abs(sumry.stats.data$coeff_var_h1)+
           abs(mean_coeff_var_h3-sumry.stats.data$coeff_var_h3)/abs(sumry.stats.data$coeff_var_h3)
  )

sumry <- bind_rows(sumry.stats.df_a0,sumry.stats.df_a02,sumry.stats.df_a04,
                   sumry.stats.df_a06,sumry.stats.df_a08,sumry.stats.df_a1)  

sumry <- sumry %>%
  mutate(
    dist2=ifelse(dist > 5, 5, dist)
  )

sumry_corrb0 <- sumry %>% filter(mean_corr<0)
ggplot(data=sumry_corrb0,aes(theta1,theta2,fill=dist2))+
  facet_wrap(~a,labeller = label_bquote(cols = a:.(a)))+
  geom_raster()+
  scale_fill_viridis_c(breaks=c(3, 4, 5),
                       labels=c("3", "4", ">5"))+
  labs(fill='Distance',x=expression(''~theta[1]),
       y=expression(''~theta[2]))+
  scale_x_continuous(expand = c(0,0),labels=c('0','0.25','0.5','0.75','1'))
  
ggsave('distance_a_theta.png')

sumry3 <- sumry %>% filter(dist<3 & mean_corr<0)
ggplot(data=sumry3,aes(theta1,theta2,fill=dist))+facet_wrap(~a)+
  geom_raster()+ scale_fill_viridis_c()+
  labs(fill='Distance',x=expression(''~theta[1]),
       y=expression(''~theta[2]))+
  scale_x_continuous(expand = c(0,0),labels=c('0','0.25','0.5','0.75','1'))
ggsave('distance_a_theta_3.png')
ggplot(sumry3,aes(theta1,theta2)) + geom_count() +scale_size_area()
ggsave('count.png')
