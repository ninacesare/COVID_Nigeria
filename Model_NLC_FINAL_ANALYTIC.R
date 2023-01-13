########################################################################################
# PROJECT       :	Leveraging AI
# SPONSOR/PI    : Elaine Nsoesie
# PROGRAM NAME  : 
# DESCRIPTION   : Modeling COVID-19 mortality across states in Nigeria
#                 This code loads analytic data and builds models for the paper
#                 Note that this code does not generate figures or load raw COVID/demographic data due to unclear data sharing restrictions
#                 
# PROGRAMMER    : Nina Cesare
# DATE WRITTEN  : 2/24/2022
########################################################################################
# INPUT FILES   : 	
# OUTPUT FILES  : 	
#######################################################################################
# MODIFICATIONS : 
#
# DATE          :
# PROGRAMMER    : 
# DESCRIPTION   : 	
#######################################################################################

rm(list = ls())

library(glm)
library(sandwich) #for robust standard errors 
#https://stats.idre.ucla.edu/r/dae/poisson-regression/
library(readxl)
library(MASS)
library(fitdistrplus)
library(dplyr)
library(ggplot2)
library(msm)
library(corrplot)
library(stargazer)
library(lmtest)

#dir <- "~/Desktop/Nigeria COVID/"
dir <- "//ad.bu.edu/bumcfiles/SPH/DCC/Dept/LeveragingAI/Spring2022/Nigeria COVID/"

### Sociodemographic data for 2019
dat3_agg <- read_excel(paste0(dir, "analytic_data/dat3_agg.xlsx"))

### COVID analytic data
time_dat_merge1_agg <- read_excel(paste0(dir, "analytic_data/time_dat_merge1_agg.xlsx"))
time_dat_merge2_agg <- read_excel(paste0(dir, "analytic_data/time_dat_merge2_agg.xlsx"))
time_dat_mergeA_agg <- read_excel(paste0(dir, "analytic_data/time_dat_mergeA_agg.xlsx"))


### generate model data for three time periods
d3_m <- merge(dat3_agg, time_dat_mergeA_agg[,-2], by = "state", all =TRUE)  # time dat merge from above (remove political region)
d3_m1 <- merge(dat3_agg, time_dat_merge1_agg[,-2], by = "state", all =TRUE)  
d3_m2 <- merge(dat3_agg, time_dat_merge2_agg[,-2], by = "state", all =TRUE)  

d3_m$mortality_rate <- d3_m$period_deaths/d3_m$population_2016 * 100000
d3_m1$mortality_rate <- d3_m1$period_deaths/d3_m1$population_2016 * 100000
d3_m2$mortality_rate <- d3_m2$period_deaths/d3_m2$population_2016 * 100000

d3_m$infection_fatality_rate <- d3_m$period_cases/d3_m$population_2016 * 1000
d3_m1$infection_fatality_rate <- d3_m1$period_cases/d3_m1$population_2016 * 1000
d3_m2$infection_fatality_rate <- d3_m2$period_cases/d3_m2$population_2016 * 1000



### Plot mortality rates ####

p_all <- ggplot(d3_m, aes(reorder(state, -mortality_rate), mortality_rate, color = political_region, fill = political_region)) + 
  geom_bar(stat = "identity") +
  xlab("State") + 
  ylab("Mortality Rate \n (Deaths per 100,000 residents)") + 
  ggtitle("Full time period: February 27th 2020 to July 25th 2021") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p_1 <- ggplot(d3_m, aes(reorder(state, -mortality_rate), mortality_rate, color = political_region, fill = political_region)) + 
  geom_bar(stat = "identity") +
  xlab("State") + 
  ylab("Mortality Rate \n (Deaths per 100,000 residents)") + 
  ggtitle("Period one: February 27th 2020 to October 24th 2020") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p_2 <- ggplot(d3_m, aes(reorder(state, -mortality_rate), mortality_rate, color = political_region, fill = political_region)) + 
  geom_bar(stat = "identity") +
  xlab("State") + 
  ylab("Mortality Rate \n (Deaths per 100,000 residents)") + 
  ggtitle("Period two: October 25th 2020 to July 25th 2021") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("//ad.bu.edu/bumcfiles/SPH/DCC/Dept/LeveragingAI/Spring2022/Nigeria COVID/Figure_2.pdf", width=6, height=10)
Rmisc::multiplot(p_all, p_1, p_2, cols = 1)
dev.off()


d3_m[order(-d3_m$population_2016),c("population_2016","state")]

head(d3_m[order(-d3_m$mortality_rate),])
head(d3_m1[order(-d3_m1$mortality_rate),])
head(d3_m2[order(-d3_m2$mortality_rate),])

###### PERIOD DEATHS ######


### check for overdisperesion #####

mean(d3_m$period_deaths)
(d3_m$period_deaths)

poisson <- glm(period_deaths ~ wealth_index_poorest_per + insured_per + median_age + urban_residence_per + population_2016_t, family = poisson(link = "log"), data=d3_m)
poisson1 <- glm(period_deaths ~ wealth_index_poorest_per + insured_per + median_age + urban_residence_per + population_2016_t, family = poisson(link = "log"), data=d3_m1)
poisson2 <- glm(period_deaths ~ wealth_index_poorest_per + insured_per + median_age + urban_residence_per + population_2016_t, family = poisson(link = "log"), data=d3_m2)

AER::dispersiontest(poisson) ## way over 7... data may need quasipoisson 
AER::dispersiontest(poisson1) ## way over 7... data may need quasipoisson 
AER::dispersiontest(poisson2) ## way over 7... data may need quasipoisson 


## negative binomial model
nb <- glm.nb(period_deaths ~ wealth_index_poorest_per + insured_per + median_age + urban_residence_per + population_2016_t, data=d3_m)
nb1 <- glm.nb(period_deaths ~ wealth_index_poorest_per + insured_per + median_age + urban_residence_per + population_2016_t, data=d3_m1)
nb2 <- glm.nb(period_deaths ~ wealth_index_poorest_per + insured_per + median_age + urban_residence_per + population_2016_t, data=d3_m1)



### Confirm that the negative binomial model is an improved fit over the poisson model
pchisq(2 * (logLik(nb) - logLik(poisson)), df = 4, lower.tail = FALSE)
pchisq(2 * (logLik(nb1) - logLik(poisson1)), df = 4, lower.tail = FALSE)
pchisq(2 * (logLik(nb2) - logLik(poisson2)), df = 4, lower.tail = FALSE)


### Incidence rate ratios with 95% confidence intervals
out <- cbind(Estimate = coef(nb), confint(nb))
out1 <- cbind(Estimate = coef(nb1), confint(nb1))
out2 <- cbind(Estimate = coef(nb2), confint(nb2))

irr <- exp(out)
irr1 <- exp(out1)
irr2 <- exp(out2)



#### Descriptive tables #####


### Table 1: Sociodemographic characteristics ####

tab1 <- d3_m[,c("political_region","state","wealth_index_poorest_per","insured_per","median_age","urban_residence_per","population_2016_t")]
#write.csv(tab1, paste0(dir, "Manuscript/table_1.csv"), row.names = FALSE)


### Table 2: Infections and deaths #####
time_dat_merge1_agg <- read_excel(paste0(dir, "analytic_data/time_dat_merge1_agg.xlsx"))
time_dat_merge2_agg <- read_excel(paste0(dir, "analytic_data/time_dat_merge2_agg.xlsx"))
time_dat_mergeA_agg <- read_excel(paste0(dir, "analytic_data/time_dat_mergeA_agg.xlsx"))

names(time_dat_merge1_agg)[c(3,4)] <- paste0(colnames(time_dat_merge1_agg)[c(3,4)], "_period1")
names(time_dat_merge2_agg)[c(3,4)] <- paste0(colnames(time_dat_merge2_agg)[c(3,4)], "_period2")
names(time_dat_mergeA_agg)[c(3,4)] <- paste0(colnames(time_dat_mergeA_agg)[c(3,4)], "_periodFull")


tab2 <- merge(time_dat_merge1_agg, time_dat_merge2_agg[,-2], by = "state")
tab2 <- merge(tab2, time_dat_mergeA_agg[,-2], by = "state")
tab2 <- merge(tab2, d3_m[,c("state","population_2016")], by = "state")

tab2$infection_fatality_rate_period1 <- tab2$period_deaths_period1/tab2$period_cases_period1 * 1000
tab2$infection_fatality_rate_period2 <- tab2$period_deaths_period2/tab2$period_cases_period2 * 1000
tab2$infection_fatality_rate_periodFull <- tab2$period_deaths_period2/tab2$period_cases_periodFull * 1000

tab2$mortality_rate_period1 <- tab2$period_deaths_period1/tab2$population_2016 * 100000
tab2$mortality_rate_period2 <- tab2$period_deaths_period2/tab2$population_2016 * 100000
tab2$mortality_rate_periodFull <- tab2$period_deaths_periodFull/tab2$population_2016 * 100000


tab2_sort <- tab2[,c("political_region",
                     "state",
                     "period_cases_period1",
                     "period_deaths_period1",
                     "infection_fatality_rate_period1",
                     "mortality_rate_period1",
                     "period_cases_period2",
                     "period_deaths_period2",
                     "infection_fatality_rate_period2",
                     "mortality_rate_period2",
                     "period_cases_periodFull",
                     "period_deaths_periodFull",
                     "infection_fatality_rate_periodFull",
                     "mortality_rate_periodFull")]

write.csv(tab2_sort, paste0(dir, "Manuscript/table_2.csv"), row.names = FALSE)



