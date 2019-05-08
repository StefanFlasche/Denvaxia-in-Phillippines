
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

# data : Hazard Ratio in seronegatives by time period (table S23), all studies, MI MO used.
df.t <- tibble(Year = rep(c("Active Phase","Year 1 of Hospital Phase","Year 2 of Hospital Phase", "Beyond Year 2 of Hospital Phase"),2),
               randomisation = rep(c("vacc","control"),each=4),
               serostatus = "neg",
               outcome = "hosp",
               n = c(14, 15.1, 16.7, 18.4,9.7,2.9,4.2,8.5),
               N = c(rep(375.1,4),rep(207.2,4)))

# plot data
p.data <- df %>% ggplot(aes(x= serostatus, y = incidence.ph.mid, ymin = incidence.ph.lo, ymax = incidence.ph.hi,color=randomisation)) +
  geom_linerange(position=position_dodge(width = 0.5)) +
  geom_point(position=position_dodge(width = 0.5)) +
  facet_grid(.~outcome, scales = "free") +
  coord_flip() + ylab("incidence per 100") +
  labs(y = "Incidence per 100", x = "Serostatus", color = "Randomisation")
ggsave(filename = "Pics\\Fig1_Data.tiff",p.data ,unit="cm", width = 14, height = 5, compression = "lzw", dpi = 300)

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

# proportion seropos
CasesAverted(seroPrevalence = .85, sensitivity = 1, specificity = 0, cohortSize=830000, df.tmp = df, outcm = "hosp") %>%
  mutate(prop_Neg_vac = CasesVaccSeroNeg/CasesVaccTotal) %>%
  mutate(prop_Neg_Novac = CasesNoVaccSeroNeg/CasesNoVaccTotal) %>%
  select("prop_Neg_vac","prop_Neg_Novac") %>%
  gather(key, value) %>% group_by(key) %>% summarise(mid = median(value),
                                                     lo = quantile(value, probs = 0.025),
                                                     hi = quantile(value, probs = 0.975)) 


#calculate excess cases over time

