rm(list=ls())
library(dplyr)
load('sim_scale2Clades_1.rdata')
rm(simdf.meta.list,simdf.meta.list.vac.dur)
stats.sum.all <- stats.sum %>%
  mutate(param=paste0(strain,'_p1'))
stats.sum.vac.dur.all <- stats.sum.vac.dur %>%
  mutate(param=paste0(strain,'_p1'))
rm(stats.sum,stats.sum.vac.dur)
load('sim_scale2Clades_2.rdata')
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
# load('sim_scale2Clades_3.rdata')
# rm(simdf.meta.list,simdf.meta.list.vac.dur)
# stats.sum <- stats.sum %>%
#   mutate(param=paste0(strain,'_p3'))
# stats.sum.vac.dur <- stats.sum.vac.dur %>%
#   mutate(param=paste0(strain,'_p3'))
# stats.sum.all <- stats.sum.all %>%
#   rbind(stats.sum)
# stats.sum.vac.dur.all <- stats.sum.vac.dur.all %>%
#   rbind(stats.sum.vac.dur)
# rm(stats.sum,stats.sum.vac.dur)

stats.sum.all2 <- stats.sum.all %>%
  mutate(susReduct=paste("Relative strength of\nvaccinal vs natural\nimmunity:",scaling))
mycolors <- rgb(c(col2rgb("darkred")[1],col2rgb("steelblue")[1],211),
                c(col2rgb("darkred")[2],col2rgb("steelblue")[2],211),
                c(col2rgb("darkred")[1],col2rgb("steelblue")[3],211),
                max=255,alpha=125,names = c("red50","blue50","grey50"))
library(ggplot2)
# png('vac_all_scale2clade.png',res=300,units='cm',width=13,height=14)
# ggplot(stats.sum.all2,aes(x=vac.cov/52*100,group=param,colour=param,
#                           linetype=param))+ geom_line(aes(y=avg.inci*1000))+
#   xlab('Vaccination rate (%) per week') +
#   ylab('Average weekly incidence in 1000') +
#   facet_grid(vaccineImmunity ~ susReduct,
#              labeller=labeller(.rows = label_value, .cols = label_value)) +
#   # label_bquote(cols = \tau_i/\theta_i: .(scaling), rows = vaccineType)) +
#   scale_color_manual(values=c('darkred','steelblue',
#                               mycolors['red50'][[1]],mycolors['blue50'][[1]],
#                               mycolors['red50'][[1]],mycolors['blue50'][[1]]),
#                      labels=c("H1 - best fit",
#                               "H3 - best fit",
#                               "H1 - moderate fit",
#                               "H3 - moderate fit"),
#                      limits=c("H1_p1","H3_p1","H1_p2","H3_p2")) +
#   scale_linetype_manual(values=c(1,1,2,2),
#                         labels=c("H1 - best fit",
#                                  "H3 - best fit",
#                                  "H1 - moderate fit",
#                                  "H3 - moderate fit"),
#                         limits=c("H1_p1","H3_p1","H1_p2","H3_p2"))+
#   scale_x_continuous(breaks = c(0,0.5,1,1.5,2),label=
#                        c('0','0.5','1','1.5','2'))+
#   theme_classic() +
#   theme(
#     legend.title = element_blank(),
#     legend.position = 'bottom',
#     # legend.direction="horizontal",
#     axis.line = element_line())+
#   guides(linetype=guide_legend(nrow=2,byrow=FALSE))+
#   guides(color=guide_legend(nrow=2,byrow=FALSE))
# 
# dev.off()
# Best fit - vaccine strength
png('vac_all_scale2clade_bestFit.png',res=300,units='cm',width=11,height=11.5)
ggplot(stats.sum.all2 %>% filter(param=="H1_p1"|param=="H3_p1"),
       aes(x=vac.cov/52*100,group=param,colour=param))+ 
  geom_line(aes(y=avg.inci*1000))+
  xlab('Vaccination rate (%) per week') +
  ylab('Average weekly incidence in 1000') +
  facet_grid(vaccineImmunity ~ susReduct,
             labeller=labeller(.rows = label_value, .cols = label_value)) +
  # label_bquote(cols = \tau_i/\theta_i: .(scaling), rows = vaccineType)) +
  scale_color_manual(values=c('darkred','steelblue',
                              mycolors['red50'][[1]],mycolors['blue50'][[1]],
                              mycolors['red50'][[1]],mycolors['blue50'][[1]]),
                     labels=c("H1",
                              "H3"),
                     limits=c("H1_p1","H3_p1")) +
  scale_x_continuous(breaks = c(0,0.5,1,1.5,2),label=
                       c('0','0.5','1','1.5','2'))+
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.9,0.88),
    # legend.direction="horizontal",
    axis.line = element_line())+
  guides(color=guide_legend(nrow=2,byrow=FALSE))

dev.off()

# Moderate fit - vaccine strength
png('vac_all_scale2clade_moderateFit.png',res=300,units='cm',width=15,height=11)
ggplot(stats.sum.all2 %>% filter(param=="H1_p2"|param=="H3_p2"),
       aes(x=vac.cov/52*100,group=param,colour=param))+ 
  geom_line(aes(y=avg.inci*1000),linetype=2)+
  xlab('Vaccination rate (%) per week') +
  ylab('Average weekly incidence in 1000') +
  facet_grid(vaccineImmunity ~ susReduct,
             labeller=labeller(.rows = label_value, .cols = label_value)) +
  # label_bquote(cols = \tau_i/\theta_i: .(scaling), rows = vaccineType)) +
  scale_color_manual(values=c('darkred','steelblue',
                              mycolors['red50'][[1]],mycolors['blue50'][[1]],
                              mycolors['red50'][[1]],mycolors['blue50'][[1]]),
                     labels=c("H1 - moderate fit",
                              "H3 - moderate fit"),
                     limits=c("H1_p2","H3_p2")) +
  scale_x_continuous(breaks = c(0,0.5,1,1.5,2),label=
                       c('0','0.5','1','1.5','2'))+
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    # legend.direction="horizontal",
    axis.line = element_line())+
  guides(color=guide_legend(nrow=2,byrow=FALSE))

dev.off()


stats.sum.vac.dur.all2 <- stats.sum.vac.dur.all %>%
  mutate(VacDur=paste("Relative duration of\nvaccinal vs natural\nimmunity:",1/scaling))

# png('vac_dur_scale2clade.png',res=300,units='cm',width=13,height=14)
# ggplot(stats.sum.vac.dur.all2,aes(x=vac.cov/52*100,group=param,colour=param,
#                                   linetype=param))+
#   geom_line(aes(y=avg.inci*1000))+
#   xlab('Vaccination rate (%) per week') +
#   ylab('Average weekly incidence in 1000')+
#   # facet_grid(vaccineType~factor(scaling2,
#   #                                   levels=rev(levels(unique(scaling2)))),
#   # labeller=labeller(.rows = label_value, .cols = label_parsed))+
#   facet_grid(vaccineImmunity~VacDur,labeller=labeller(.rows = label_value, .cols = label_value))+
#   scale_color_manual(values=c('darkred','steelblue',
#                               mycolors['red50'][[1]],mycolors['blue50'][[1]],
#                               mycolors['red50'][[1]],mycolors['blue50'][[1]]),
#                      labels=c("H1 - best fit",
#                               "H3 - best fit",
#                               "H1 - moderate fit",
#                               "H3 - moderate fit"),
#                      limits=c("H1_p1","H3_p1","H1_p2","H3_p2")) +
#   scale_linetype_manual(values=c(1,1,2,2,3,3),
#                         labels=c("H1 - best fit",
#                                  "H3 - best fit",
#                                  "H1 - moderate fit",
#                                  "H3 - moderate fit"),
#                         limits=c("H1_p1","H3_p1","H1_p2","H3_p2"))+
#   scale_x_continuous(breaks = c(0,0.5,1,1.5,2),label=
#                        c('0','0.5','1','1.5','2'))+
#   theme_classic()+
#   theme(
#     legend.title = element_blank(),
#     legend.position = "bottom",
#     axis.line = element_line())+
#   guides(linetype=guide_legend(nrow=2,byrow=TRUE))+
#   guides(color=guide_legend(nrow=2,byrow=TRUE))
# dev.off()

# Best fit - vaccine duration
png('vac_dur_scale2clade_bestFit.png',res=300,units='cm',width=11,height=11.5)
ggplot(stats.sum.vac.dur.all2 %>% filter(param=='H1_p1'|param=='H3_p1'),
       aes(x=vac.cov/52*100,group=param,colour=param))+
  geom_line(aes(y=avg.inci*1000))+
  xlab('Vaccination rate (%) per week') +
  ylab('Average weekly incidence in 1000')+
  # facet_grid(vaccineType~factor(scaling2,
  #                                   levels=rev(levels(unique(scaling2)))),
  # labeller=labeller(.rows = label_value, .cols = label_parsed))+
  facet_grid(vaccineImmunity~VacDur,labeller=labeller(.rows = label_value, .cols = label_value))+
  scale_color_manual(values=c('darkred','steelblue',
                              mycolors['red50'][[1]],mycolors['blue50'][[1]],
                              mycolors['red50'][[1]],mycolors['blue50'][[1]]),
                     labels=c("H1",
                              "H3"),
                     limits=c("H1_p1","H3_p1")) +
  scale_linetype_manual(values=c(1,1,2,2,3,3),
                        labels=c("H1",
                                 "H3"),
                        limits=c("H1_p1","H3_p1"))+
  scale_x_continuous(breaks = c(0,0.5,1,1.5,2),label=
                       c('0','0.5','1','1.5','2'))+
  theme_classic()+
  theme(
    legend.title = element_blank(),
    legend.position = c(0.88,0.85),
    axis.line = element_line())+
  guides(color=guide_legend(nrow=2,byrow=TRUE))
dev.off()

# Moderate fit - vaccine strength
png('vac_dur_scale2clade_moderateFit.png',res=300,units='cm',width=15,height=11)
ggplot(stats.sum.vac.dur.all2 %>% filter(param=="H1_p2"|param=="H3_p2"),
       aes(x=vac.cov/52*100,group=param,colour=param))+ 
  geom_line(aes(y=avg.inci*1000),linetype=2)+
  xlab('Vaccination rate (%) per week') +
  ylab('Average weekly incidence in 1000') +
  facet_grid(vaccineImmunity~VacDur,labeller=labeller(.rows = label_value, .cols = label_value))+
  # label_bquote(cols = \tau_i/\theta_i: .(scaling), rows = vaccineType)) +
  scale_color_manual(values=c('darkred','steelblue',
                              mycolors['red50'][[1]],mycolors['blue50'][[1]],
                              mycolors['red50'][[1]],mycolors['blue50'][[1]]),
                     labels=c("H1 - moderate fit",
                              "H3 - moderate fit"),
                     limits=c("H1_p2","H3_p2")) +
  scale_x_continuous(breaks = c(0,0.5,1,1.5,2),label=
                       c('0','0.5','1','1.5','2'))+
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    # legend.direction="horizontal",
    axis.line = element_line())+
  guides(color=guide_legend(nrow=2,byrow=FALSE))

dev.off()
