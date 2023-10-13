rm(list=ls())
library(dplyr)
load('sim_OldModel_NewParms1.rdata')
rm(simdf.meta.list,simdf.meta.list.vac.dur)
stats.sum.all <- stats.sum %>%
  mutate(param=paste0(strain,'_p1'))
stats.sum.vac.dur.all <- stats.sum.vac.dur %>%
  mutate(param=paste0(strain,'_p1'))
rm(stats.sum,stats.sum.vac.dur)
load('sim_OldModel_NewParms2.rdata')
rm(simdf.meta.list,simdf.meta.list.vac.dur)
stats.sum <- stats.sum %>%
  mutate(param=paste0(strain,'_p2'))
stats.sum.vac.dur <- stats.sum.vac.dur %>%
  mutate(param=paste0(strain,'_p2'))
stats.sum.all <- stats.sum.all %>%
  rbind(stats.sum)
stats.sum.vac.dur.all <- stats.sum.vac.dur.all %>%
  rbind(stats.sum.vac.dur)

rm(stats.sum,stats.sum.vac.dur)
load('sim_OldModel_NewParms3.rdata')
rm(simdf.meta.list,simdf.meta.list.vac.dur)
stats.sum <- stats.sum %>%
  mutate(param=paste0(strain,'_p3'))
stats.sum.vac.dur <- stats.sum.vac.dur %>%
  mutate(param=paste0(strain,'_p3'))
stats.sum.all <- stats.sum.all %>%
  rbind(stats.sum)
stats.sum.vac.dur.all <- stats.sum.vac.dur.all %>%
  rbind(stats.sum.vac.dur)
rm(stats.sum,stats.sum.vac.dur)

# stats.sum.sub <- stats.sum.all %>%
#   filter(param=="H1_p5" | param=="H3_p5" |
#            param=="H1_p9_1" | param=="H3_p9_1") %>%
#   mutate(param2=param)
# 
# stats.sum.sub2 <- stats.sum.all %>%
#   filter(param!="H1_p1" & param!="H3_p1" & param!="H1_p5" & param!="H3_p5"
#          & param!="H1_p9_1" & param!="H3_p9_1" ) %>%
#   mutate(param2=paste0('other'))
# 
# stats.sum.all2 <- rbind(stats.sum.sub,stats.sum.sub2)

stats.sum.all2 <- stats.sum.all %>%
  mutate(vaccineType=case_when(vaccineImmunity=='tau1=b*theta1\ntau2=b'~'H1 vaccine',
                               vaccineImmunity=='tau1=b\ntau2=b*theta2'~'H3 vaccine',
                               vaccineImmunity=='tau1=b\ntau2=b'~'Bivalent vaccines')) %>%
  mutate(susReduct=paste("Relative strength of\nvaccinal vs natural\nimmunity:",scaling))


# stats.sum.dur.sub <- stats.sum.vac.dur.all %>%
#   filter(param=="H1_p5" | param=="H3_p5" |
#            param=="H1_p9_1" | param=="H3_p9_1") %>%
#   mutate(param2=param)
# 
# stats.sum.dur.sub2 <- stats.sum.vac.dur.all %>%
#   filter(param!="H1_p1" & param!="H3_p1" & param!="H1_p5" & param!="H3_p5"
#          & param!="H1_p9_1" & param!="H3_p9_1" ) %>%
#   mutate(param2=paste0('other'))

mycolors <- rgb(c(col2rgb("darkred")[1],col2rgb("steelblue")[1],211),
    c(col2rgb("darkred")[2],col2rgb("steelblue")[2],211),
    c(col2rgb("darkred")[1],col2rgb("steelblue")[3],211),
    max=255,alpha=125,names = c("red50","blue50","grey50"))
library(ggplot2)
# library(tikzDevice)
png('vac_all.png',res=300,units='cm',width=13,height=14)
ggplot(stats.sum.all2,aes(x=vac.cov/52*100,group=param,colour=param,
                          linetype=param))+ geom_line(aes(y=avg.inci*1000))+
  xlab('Vaccination rate (%) per week') +
  ylab('Average weekly incidence in 1000') +
  facet_grid(vaccineType ~ susReduct,
             labeller=labeller(.rows = label_value, .cols = label_value)) +
  # label_bquote(cols = \tau_i/\theta_i: .(scaling), rows = vaccineType)) +
  scale_color_manual(values=c('darkred','steelblue',
                              mycolors['red50'][[1]],mycolors['blue50'][[1]],
                              mycolors['red50'][[1]],mycolors['blue50'][[1]]),
                     labels=c("H1 - best fit",
                              "H3 - best fit",
                              "H1 - moderate fit 1",
                              "H3 - moderate fit 1", 
                              "H1 - moderate fit 2",
                              "H3 - moderate fit 2"),
                     limits=c("H1_p1","H3_p1","H1_p2","H3_p2","H1_p3","H3_p3")) +
  scale_linetype_manual(values=c(1,1,2,2,3,3),
                        labels=c("H1 - best fit",
                                 "H3 - best fit",
                                 "H1 - moderate fit 1",
                                 "H3 - moderate fit 1", 
                                 "H1 - moderate fit 2",
                                 "H3 - moderate fit 2"),
                        limits=c("H1_p1","H3_p1","H1_p2","H3_p2","H1_p3","H3_p3"))+
  scale_x_continuous(breaks = c(0,0.5,1,1.5,2),label=
                       c('0','0.5','1','1.5','2'))+
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = 'bottom',
    # legend.direction="horizontal",
    axis.line = element_line())+
  guides(linetype=guide_legend(nrow=2,byrow=FALSE))+
  guides(color=guide_legend(nrow=2,byrow=FALSE))

dev.off()

stats.sum.vac.dur.all2 <- stats.sum.vac.dur.all %>%
  mutate(vaccineType=case_when(vaccineImmunity=='tau1=b*theta1\ntau2=b\nsigmaV=b*sigma'~'H1 vaccine',
                               vaccineImmunity=='tau1=b\ntau2=b*theta2\nsigmaV=b*sigma'~'H3 vaccine',
                               vaccineImmunity=='tau1=theta1\ntau2=theta2\nsigmaV=b*sigma'~'Bivalent vaccines')) %>%
  mutate(VacDur=paste("Relative duration of\nvaccinal vs natural\nimmunity:",1/scaling))

png('vac_dur_all.png',res=300,units='cm',width=13,height=14)
ggplot(stats.sum.vac.dur.all2,aes(x=vac.cov/52*100,group=param,colour=param,
                                  linetype=param))+
  geom_line(aes(y=avg.inci*1000))+
  xlab('Vaccination rate (%) per week') +
  ylab('Average weekly incidence in 1000')+
  # facet_grid(vaccineType~factor(scaling2,
  #                                   levels=rev(levels(unique(scaling2)))),
             # labeller=labeller(.rows = label_value, .cols = label_parsed))+
  facet_grid(vaccineType~VacDur,labeller=labeller(.rows = label_value, .cols = label_value))+
  scale_color_manual(values=c('darkred','steelblue',
                              mycolors['red50'][[1]],mycolors['blue50'][[1]],
                              mycolors['red50'][[1]],mycolors['blue50'][[1]]),
                     labels=c("H1 - best fit",
                              "H3 - best fit",
                              "H1 - moderate fit 1",
                              "H3 - moderate fit 1", 
                              "H1 - moderate fit 2",
                              "H3 - moderate fit 2"),
                     limits=c("H1_p1","H3_p1","H1_p2","H3_p2","H1_p3","H3_p3")) +
  scale_linetype_manual(values=c(1,1,2,2,3,3),
                        labels=c("H1 - best fit",
                                 "H3 - best fit",
                                 "H1 - moderate fit 1",
                                 "H3 - moderate fit 1", 
                                 "H1 - moderate fit 2",
                                 "H3 - moderate fit 2"),
                        limits=c("H1_p1","H3_p1","H1_p2","H3_p2","H1_p3","H3_p3"))+
  scale_x_continuous(breaks = c(0,0.5,1,1.5,2),label=
                       c('0','0.5','1','1.5','2'))+
  theme_classic()+
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.line = element_line())+
  guides(linetype=guide_legend(nrow=2,byrow=TRUE))+
  guides(color=guide_legend(nrow=2,byrow=TRUE))
dev.off()
