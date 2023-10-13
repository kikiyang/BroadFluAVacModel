setwd('~/Documents/P3_UIV/BroadFluVacModel/data/USA')
pop = read.csv(file=c('Monthly_population/natmonthl_2000_2020.csv'))
library(lubridate)
pop$year <- year(mdy(pop$Month))
pop$month <- month(mdy(pop$Month))
pop = pop[pop$year>2003,]

ILI = read.csv('ILINet_CLEAN.csv',header = TRUE)
ILI = ILI[ILI$YEAR>2003 & ILI$YEAR<2021,]
ILI$prop <- ILI$ILITOTAL/ILI$TOTAL.PATIENTS

# virus surveillance
viral15 = read.csv('WHO_NREVSS_Combined_prior_to_2015_16_CLEAN.csv')
viral_ph = read.csv("WHO_NREVSS_Public_Health_Labs_CLEAN.csv")
viral_cl = read.csv("WHO_NREVSS_Clinical_Labs_CLEAN.csv")

viral04_15 <- viral15[viral15$YEAR>2003,]
viral04_15$pos_prop <- viral04_15$PERCENT.POSITIVE/100
# library(dplyr)
# viral_ph %>%
#   filter(.$H3N2v>.$A..H3.)

library(dplyr)
viral_ph <- viral_ph %>% select(!(REGION.TYPE:REGION))
viral_cl <- viral_cl %>% select(!(REGION.TYPE:REGION))
viral_cl$pos_prop <- viral_cl$PERCENT.POSITIVE/100

viral_ph <- viral_ph %>%
  mutate(pos_specimens=A..2009.H1N1.+A..H3.+
           A..Subtyping.not.Performed.+B+BVic+BYam+H3N2v,
         pos_prop=pos_specimens/TOTAL.SPECIMENS,
         subtyped.A.H1=A..2009.H1N1.,subtyped.A.H3=A..H3.+H3N2v,subtyped.B=B+BVic+BYam) %>%
  mutate(ratio.A.H1=subtyped.A.H1/(subtyped.A.H1+subtyped.A.H3),
         ratio.A.H3=subtyped.A.H3/(subtyped.A.H1+subtyped.A.H3))

viral15_21 <- viral_ph %>% 
  inner_join(viral_cl,by=c('YEAR','WEEK')) %>%
  mutate(TOTAL.SPECIMENS=TOTAL.SPECIMENS.x+TOTAL.SPECIMENS.y,
         pos_specimens=pos_specimens+TOTAL.A+TOTAL.B,
         UNSUBTYPED.A.H1.EST=(A..Subtyping.not.Performed.+TOTAL.A)*ratio.A.H1,
         UNSUBTYPED.A.H3.EST=(A..Subtyping.not.Performed.+TOTAL.A)*ratio.A.H3) %>%
  mutate(PROP.A.H1.EST=(UNSUBTYPED.A.H1.EST+subtyped.A.H1)/pos_specimens,
         PROP.A.H3.EST=(UNSUBTYPED.A.H3.EST+subtyped.A.H3)/pos_specimens,
         PROP.B.EST=(TOTAL.B+subtyped.B)/pos_specimens,
         pos_prop=pos_specimens/TOTAL.SPECIMENS)
  

viral04_15<- viral04_15 %>% 
  select(!(REGION.TYPE:REGION)) %>%
  mutate(pos_specimens=A..2009.H1N1.+A..H1.+A..H3.+A..Subtyping.not.Performed.+
           A..Unable.to.Subtype.+B+H3N2v,pos_prop=pos_specimens/TOTAL.SPECIMENS,
         subtyped.A.H1=A..2009.H1N1.+A..H1.,subtyped.A.H3=A..H3.+H3N2v,
         ratio.A.H1=subtyped.A.H1/(subtyped.A.H1+subtyped.A.H3),
         ratio.A.H3=subtyped.A.H3/(subtyped.A.H1+subtyped.A.H3),
         UNSUBTYPED.A.H1.EST=(A..Subtyping.not.Performed.+A..Unable.to.Subtype.)*ratio.A.H1,
         UNSUBTYPED.A.H3.EST=(A..Subtyping.not.Performed.+A..Unable.to.Subtype.)*ratio.A.H3,
         PROP.A.H1.EST=(UNSUBTYPED.A.H1.EST+subtyped.A.H1)/pos_specimens,
         PROP.A.H3.EST=(UNSUBTYPED.A.H3.EST+subtyped.A.H3)/pos_specimens,
         PROP.B.EST=B/pos_specimens)




viral04_15 <- viral04_15 %>%
  mutate(unsubtyped.A=A..Subtyping.not.Performed.+A..Unable.to.Subtype.) %>%
  select(!(A..2009.H1N1.:A..Unable.to.Subtype.))

cols <- unlist(colnames(viral04_15))

viral15_21 <- viral15_21 %>%
  mutate(unsubtyped.A=A..Subtyping.not.Performed.) %>%
  filter(YEAR<2021) %>%
  select(cols)
  
rm(list = c('viral_cl','viral_ph','viral15'))

# Incidence = Population size * Outpatient data (Ratio of ILI patients among all hospital visits) * 
# Ratio of influenza positive sample among all testings * Ratio of specific subtype among all positive samples
viral04_20 <- viral04_15 %>% 
  bind_rows(viral15_21) %>%
  inner_join(ILI,by=c('YEAR','WEEK')) %>%
  select(!(REGION.TYPE:REGION))

inci <- viral04_20 %>%
  mutate(month=month(ymd(paste0(YEAR,'-01-01'))+weeks(WEEK-1))) %>%
  left_join(pop,by=c('YEAR'='year','month')) %>%
  select(!Month) %>%
  rename(ILI.prop=X..WEIGHTED.ILI) %>%
  mutate(cases.A.H1=Resident.population*ILI.prop*pos_prop*PROP.A.H1.EST,
         cases.A.H3=Resident.population*ILI.prop*pos_prop*PROP.A.H3.EST,
         cases.B=Resident.population*ILI.prop*pos_prop*PROP.B.EST,
         inci.A.H1=ILI.prop*pos_prop*PROP.A.H1.EST,
         inci.A.H3=ILI.prop*pos_prop*PROP.A.H3.EST,
         inci.B=ILI.prop*pos_prop*PROP.B.EST)

# write.csv(inci,file='no.cases_inci.allsubtypes.csv')

library(tidyverse)
inci4pl <- inci %>%
  select(YEAR,WEEK,inci.A.H1,inci.A.H3,inci.B) %>%
  pivot_longer(cols=c(inci.A.H1,inci.A.H3,inci.B),names_to = 'subtype',values_to='inci')

# incidence A H1
incih1 <- inci %>%
  select(YEAR,WEEK,inci.A.H1)
# incidence A H3
incih3 <- inci %>%
  select(YEAR,WEEK,inci.A.H3)
# incidence B
inciB <- inci %>%
  select(YEAR,WEEK,inci.B)

# write.csv(inci4pl,'incidence_allsubtypes.csv')
# write.csv(incih1,'incidence_USA_h1.csv')
# write.csv(incih3,'incidence_USA_h3.csv')
# write.csv(inciB,'incidence_USA_B.csv')

ggplot(data=inci4pl,mapping = aes(x=ymd(paste0(YEAR,'-01-01'))+weeks(WEEK-1),y=inci,color=subtype)) +
       geom_line() + 
  scale_color_manual(labels= c('H1','H3','B'), values = c('darkred','steelblue','darkgrey'))+
  xlab('Time (weekly)') + 
  ylab('Incidence') + 
  theme_bw()
ggsave('USA.Inci.png',width = 8,height = 4)

# Subtype A
inciA<-inci4pl %>%
  filter(subtype!='inci.B')
ggplot(data=inciA,mapping = aes(x=ymd(paste0(YEAR,'-01-01'))+weeks(WEEK-1),
                                y=inci,color=subtype)) + 
  annotate("rect", xmin = as.Date('2008-07-01') , xmax = as.Date("2011-07-01"), 
           ymin = 0, ymax = 3.1, alpha = .5,fill = "grey") + 
  geom_line() + 
  scale_color_manual(labels= c('H1','H3'), values = c('darkred','steelblue','darkgrey'))+
  xlab('Time (weekly)') +
  ylab('Incidence') + 
  theme_bw()
ggsave('USA.InciA.png',width = 8,height = 4)

       