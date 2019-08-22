
# Expected relative proportion of breakthrough hospitalized dengue by serostatus in populations vaccinated with the CYD- tetravalent dengue vaccine in the Philippines
# written by Stefan Flasche using R version 3.5.3

# load packages
require(tidyverse)


# Data --------------------------------------------------------------------

# data: cumulative incidence over 5yrs observed in the trial and its stratification into 
#       seropositive and seronegative using PRNT50 with multiple imputation (in per 100). from Shridar Table S10
df <- tibble(no=1:8,
             randomisation = c("vacc","vacc","control","control","vacc","vacc","control","control"),
             serostatus = c("pos","neg","pos","neg","pos","neg","pos","neg"),
             outcome = c(rep("hosp",4),rep("severe",4)),
             incidence.ph.mid = c(.375,1.571,1.883,1.093,0.075,0.404,0.48,0.174),
             incidence.ph.lo = c(.263,1.125,1.536,.526,0.034,0.218,0.335,0.036),
             incidence.ph.hi = c(.535,2.193,2.307,2.265,0.165,0.749,0.688,0.834))

#seroprevalence Phil
sp=.70

# Cumulative proportion of 9-16y olds hospitalised with Dengue over time (0m - 60m). Digitalised from Sridhar Figure 3
df.time <- tibble(time = rep(c(0,6,12,18,24,30,36,42,48,54,60),4),
                  randomisation = c(rep("vacc",11*2),rep("control",11*2)),
                  serostatus = rep(c(rep("pos",11),rep("neg",11)),2),
                  incidence = c(0,0.031,0.051,0.061,0.113,0.153,0.184,0.235,0.276,0.307,0.368,
                                0,0.061,0.143,0.205,0.348,0.512,0.757,0.982,1.289,1.432,1.535,
                                0,0.143,0.307,0.573,0.982,1.125,1.350,1.514,1.678,1.739,1.862,
                                0,0.082,0.164,0.399,0.471,0.522,0.583,0.696,0.757,0.757,1.033),
                  No.at.Risk =c(1503,1490,1462,1453,1447,1434,1422,1403,1351,1284,1226,
                                375,374,368,366,364,361,358,350,339,322,311,
                                730,725,705,699,695,689,688,686,658,615,592,
                                207,207,201,199,199,195,194,192,185,174,167))
                  
# plot data
p.data <- df %>% ggplot(aes(x= serostatus, y = incidence.ph.mid, ymin = incidence.ph.lo, ymax = incidence.ph.hi,color=randomisation)) +
  geom_linerange(position=position_dodge(width = 0.5)) +
  geom_point(position=position_dodge(width = 0.5)) +
  facet_grid(.~outcome, scales = "free") +
  coord_flip() + ylab("incidence per 100") +
  labs(y = "Incidence per 100", x = "Serostatus", color = "Randomisation")
ggsave(filename = "Pics\\Fig_Data.tiff",p.data ,unit="cm", width = 14, height = 5, compression = "lzw", dpi = 300)


### analyses on cumulative 5y estimates
vph = subset(df, randomisation=="vacc" & serostatus == "pos" & outcome =="hosp")$incidence.ph.mid
vnh = subset(df, randomisation=="vacc" & serostatus == "neg" & outcome =="hosp")$incidence.ph.mid
cph = subset(df, randomisation=="control" & serostatus == "pos" & outcome =="hosp")$incidence.ph.mid
cnh = subset(df, randomisation=="control" & serostatus == "neg" & outcome =="hosp")$incidence.ph.mid
vps = subset(df, randomisation=="vacc" & serostatus == "pos" & outcome =="severe")$incidence.ph.mid
vns = subset(df, randomisation=="vacc" & serostatus == "neg" & outcome =="severe")$incidence.ph.mid
cps = subset(df, randomisation=="control" & serostatus == "pos" & outcome =="severe")$incidence.ph.mid
cns = subset(df, randomisation=="control" & serostatus == "neg" & outcome =="severe")$incidence.ph.mid

#increase in seroneg vaccinees
vnh/cnh
vns/cns #severe

#RR in seropo / seroneg controls
cph/cnh
cps/cns #severe

#reduction in seropos vaccinees
1-vph/cph
1-vps/cps#severe

#cases averted in seropos vacc per 1 seroneg vacc case
sp*(cph-vph) / ((1-sp)*(vnh-cnh))
sp*(cps-vps) / ((1-sp)*(vns-cns))#severe

# breakthrough among all hospitalised
sp*vph  / (sp*vph + (1-sp)*vnh)
sp*vps  / (sp*vps + (1-sp)*vns)#severe

#proportion of excess cases among the cases in seronegatives
1-cnh/vnh
1-cns/vns#severe

#proportion of excess cases among all cases
(1-(sp*vph  / (sp*vph + (1-sp)*vnh))) * (1-cnh/vnh)
(1-(sp*vps  / (sp*vps + (1-sp)*vns))) * (1-cns/vns)#severe

# proportion of cases averted through vacc
1 - (sp*vph+(1-sp)*vnh) / (sp*cph+(1-sp)*cnh)
1 - (sp*vps+(1-sp)*vns) / (sp*cps+(1-sp)*cns)


#calculate  time dependent proportion of cases attributeable to breakthrough in seropositive, and excess cases n seronegative vaccinees
  #proportion of excess cases among the cases in seronegatives by time
  excess.prop = 1 - subset(df.time, randomisation=="control" & serostatus=="neg")$incidence / subset(df.time, randomisation=="vacc" & serostatus=="neg")$incidence
  excess.prop = excess.prop * (excess.prop>0)
  
df.prop = tibble(serostatus = c("neg","pos"),
                 prop = c(1-sp,sp))
data.time = df.time %>% merge(df.prop) %>% 
  mutate(cases = incidence*prop*830000/100) %>% 
  select(serostatus:randomisation, cases) %>%
  spread(serostatus, cases) %>%
  mutate(propSeroNeg = neg/(pos+neg)) %>%
  filter(randomisation=="vacc") %>%
  mutate(propVaccINduced = propSeroNeg * excess.prop) %>%
  mutate(propSeroNeg = propSeroNeg - propVaccINduced) %>%
  mutate(propSeroPos = 1-propSeroNeg-propVaccINduced) %>%
  gather(key = "key", value="prop", -(time:pos)) 
data.time$key = factor(data.time$key, levels = c("propSeroPos","propSeroNeg","propVaccINduced"))
data.time %>%
  ggplot(aes(x=time, y=prop, fill=key))+ 
  geom_area(alpha=.7) + 
  coord_cartesian(ylim=c(0,1)) + 
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#009E73","#E69F00","#D55E00"), name = "", labels = c("breakthrough cases in\nseropositive vaccinees\n", "cases in\nseronegative vaccinees\n", "vaccine attributable cases in\nseronegative vaccinees\n"))+
  theme_bw() +
  xlab("month since first vaccine dose") + ylab("percentage of all\nhospitalised cases")
ggsave(filename = "Pics\\Fig_PropCases.tiff",unit="cm", width = 17, height = 7, compression = "lzw", dpi = 300)



