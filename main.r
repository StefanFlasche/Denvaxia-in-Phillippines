
# Expected relative proportion of breakthrough hospitalized dengue by serostatus in populations vaccinated with the CYD- tetravalent dengue vaccine in the Philippines
# written by Stefan Flasche using R version 3.5.3

# load packages
require(tidyverse)


# Data --------------------------------------------------------------------

# data: cumulative incidence over 5yrs observed in the trial and its stratification into 
#       seropositive and seronegative using PRNT50 with multiple imputation (in per 100). from Table S10
df <- tibble(no=1:8,
             randomisation = c("vacc","vacc","control","control","vacc","vacc","control","control"),
             serostatus = c("pos","neg","pos","neg","pos","neg","pos","neg"),
             outcome = c(rep("hosp",4),rep("severe",4)),
             incidence.ph.mid = c(.375,1.571,1.883,1.093,0.075,0.404,0.48,0.174),
             incidence.ph.lo = c(.263,1.125,1.536,.526,0.034,0.218,0.335,0.036),
             incidence.ph.hi = c(.535,2.193,2.307,2.265,0.165,0.749,0.688,0.834))

#seroprevalence Phil
sp=.85

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
ggsave(filename = "Pics\\Fig1_Data.tiff",p.data ,unit="cm", width = 14, height = 5, compression = "lzw", dpi = 300)

# analyses

#increase in seroneg vaccinees
1.57/1.09
.404/.174#severe
#RR in seropo / seroneg controls
1.88/1.09
.48/.174#severe
#reduction in seropos vaccinees
1-0.375/1.88
1-.075/.48#severe

#cases averted in seropos vacc per 1 seroneg vacc case
.85*(1.88-.375) / (.15*(1.57-1.09))
.85*(.48-.075) / (.15*(.404-.174))#severe
# breakthrough among all hospitalised
0.85*.375  / (.85*.375 + .15*1.57)
0.85*.075  / (.85*.075 + .15*.404)#severe
#proportion of excess cases among the cases in seronegatives
1-1.09/1.57
1-.174/.404#severe

#calculate  proportion of cases attributeable to seronegative vaccinees
df.prop = tibble(serostatus = c("neg","pos"),
                 prop = c(0.15,.85))
df.time %>% merge(df.prop) %>% 
  mutate(cases = incidence*prop*830000/100) %>% 
  select(serostatus:randomisation, cases) %>%
  spread(serostatus, cases) %>%
  mutate(propSeroNeg = neg/(pos+neg)) %>%
  filter(randomisation=="vacc") %>%
  mutate(propVaccINduced = propSeroNeg * (1-1.09/1.57)) %>%
  gather(key = "key", value="prop", -(time:pos)) %>%
  ggplot(aes(x=time, y=prop, color=key))+ 
  geom_line() + 
  coord_cartesian(ylim=c(0,1)) + 
  scale_y_continuous(labels = scales::percent) +
    scale_color_discrete(name = "", labels = c("seronegative vaccinees", "vaccine induced"))+
  theme_bw() +
  xlab("Month since vaccination") + ylab("proportion of all\nhospitalised cases")
ggsave(filename = "Pics\\Fig3_Data.tiff",unit="cm", width = 17, height = 7, compression = "lzw", dpi = 300)



##################### end of numbers in manuscript # start of previous work








# sample from lognormal distribution in which the 50%, 2.5% and 97.5% quanitles fit the observed mean, CI.lo and CI.hi respectively
LnfitSample <- function(input= c(mean=1, lo=.5, hi=2), N=10000){
  mean = input[1] %>% as.numeric()
  lo = input[2] %>% as.numeric()
  hi = input[3] %>% as.numeric()
  rlnorm(N, meanlog = log(mean), sdlog = (log(mean/lo) +  log(hi/mean))/2/1.96) %>%
    return()
}


# Analyses ----------------------------------------------------------

# calculate population impact of test and vaccinate strategy (no test = 100% sens, 0% spec)
CasesAverted <- function(seroPrevalence = .7, sensitivity = 1, specificity = 0, 
                         cohortSize=100000, df.tmp = df, outcm = "hosp"){
  df.tmp = df.tmp %>% filter(outcome == outcm)
  Inc.SeroPos.Vacc <- df.tmp %>% filter(randomisation=="vacc" & serostatus =="pos") %>% 
    select(incidence.ph.mid:incidence.ph.hi) %>% LnfitSample()
  Inc.SeroPos.Cont <- df.tmp %>% filter(randomisation=="control" & serostatus =="pos") %>% 
    select(incidence.ph.mid:incidence.ph.hi) %>% LnfitSample()
  Inc.SeroNeg.Vacc <- df.tmp %>% filter(randomisation=="vacc" & serostatus =="neg") %>% 
    select(incidence.ph.mid:incidence.ph.hi) %>% LnfitSample()
  Inc.SeroNeg.Cont <- df.tmp %>% filter(randomisation=="control" & serostatus =="neg") %>% 
    select(incidence.ph.mid:incidence.ph.hi) %>% LnfitSample()
  
  CasesNoVaccSeroPos = cohortSize * seroPrevalence * Inc.SeroPos.Cont / 100
  CasesNoVaccSeroNeg = cohortSize * (1-seroPrevalence) * Inc.SeroNeg.Cont / 100
  CasesVaccSeroPos = cohortSize * seroPrevalence * Inc.SeroPos.Vacc / 100
  CasesVaccSeroNeg = cohortSize * (1-seroPrevalence) * Inc.SeroNeg.Vacc / 100
  
  CasesAvertedSeroPos = cohortSize * seroPrevalence * sensitivity * (Inc.SeroPos.Cont - Inc.SeroPos.Vacc) / 100
  CasesAvertedSeroNeg = cohortSize * (1-seroPrevalence) * (1 - specificity) * (Inc.SeroNeg.Cont - Inc.SeroNeg.Vacc) /100
  
  df_res <- tibble(seroPrevalence, sensitivity, specificity, outcm, 
                   CasesNoVaccSeroPos, CasesNoVaccSeroNeg,
                   CasesNoVaccTotal = CasesNoVaccSeroPos + CasesNoVaccSeroNeg,
                   CasesVaccSeroPos, CasesVaccSeroNeg,
                   CasesVaccTotal = CasesVaccSeroPos + CasesVaccSeroNeg,
                   CasesAvertedSeroPos, CasesAvertedSeroNeg, 
                   CasesAvertedTotal = CasesAvertedSeroPos + CasesAvertedSeroNeg)
  return(df_res)
}

# calculate cases exspected and averted in Phillippines
res <- CasesAverted(seroPrevalence = .85, sensitivity = 1, specificity = 0, cohortSize=830000, df.tmp = df, outcm = "hosp") %>%
  gather(key, value, -(seroPrevalence:outcm)) %>% 
  group_by(key) %>% summarise(mid = median(value), 
                              lo = quantile(value, probs = 0.025),
                              hi = quantile(value, probs = 0.975)) %>% 
  mutate(outcome = rep(c("Averted","Control","Vacc"),each=3)) %>%
  mutate(serostatus = rep(c("Negative", "Positive", "Total"), 3))
res
res %>% filter(outcome != "Averted") %>%
  ggplot(aes(x=serostatus, y=mid, ymin=lo, ymax=hi, color = outcome)) +
    geom_pointrange(position=position_dodge(width = 0.5)) + 
    coord_flip() + ylab("Number of cases") + 
    theme_bw()
ggsave(filename = "Pics\\Fig2_Data.tiff",unit="cm", width = 17, height = 7, compression = "lzw", dpi = 300)



# proportion seropos
CasesAverted(seroPrevalence = .85, sensitivity = 1, specificity = 0, cohortSize=830000, df.tmp = df, outcm = "hosp") %>%
  mutate(prop_Neg_vac = CasesVaccSeroNeg/CasesVaccTotal) %>%
  mutate(prop_Neg_Novac = CasesNoVaccSeroNeg/CasesNoVaccTotal) %>%
  select("prop_Neg_vac","prop_Neg_Novac") %>%
  gather(key, value) %>% group_by(key) %>% summarise(mid = median(value),
                                                     lo = quantile(value, probs = 0.025),
                                                     hi = quantile(value, probs = 0.975)) 



