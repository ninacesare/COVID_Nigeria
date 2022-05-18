########################################################################################
# PROJECT       :	Leveraging AI
# SPONSOR/PI    : Elaine Nsoesie
# PROGRAM NAME  : 
# DESCRIPTION   : Modeling COVID-19 mortality across states in Nigeria
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

#dir <- "~/Desktop/Nigeria COVID/"
dir <- "//ad.bu.edu/bumcfiles/SPH/DCC/Dept/LeveragingAI/Spring2022/Nigeria COVID/"

### Export analytic datasets 
time_dat_merge1_agg <- read_excel(paste0(dir, "analytic_data/time_dat_merge1_agg.xlsx"))
time_dat_merge2_agg <- read_excel(paste0(dir, "analytic_data/time_dat_merge2_agg.xlsx"))
time_dat_mergeA_agg <- read_excel(paste0(dir, "analytic_data/time_dat_mergeA_agg.xlsx"))


### generate model data for three time periods
d3_m <- merge(dat3_agg, time_dat_mergeA_agg, by = "state", all =TRUE)  # time dat merge from above
d3_m1 <- merge(dat3_agg, time_dat_merge1_agg, by = "state", all =TRUE)  # time dat merge from above
d3_m2 <- merge(dat3_agg, time_dat_merge2_agg, by = "state", all =TRUE)  # time dat merge from above

## Build models
m3_m_all <- glm(period_deaths ~ wealth_index_poorest_per + insured_per + median_age + urban_residence_per + period_cases + population_2016_t, family = "poisson", data=d3_m)

m3_m_t1 <- glm(period_deaths ~ wealth_index_poorest_per + insured_per + median_age + urban_residence_per + period_cases + population_2016_t, family = "poisson", data=d3_m1)
m3_m_t2 <- glm(period_deaths ~ wealth_index_poorest_per + insured_per + median_age + urban_residence_per + period_cases + population_2016_t, family = "poisson", data=d3_m2)


stargazer::stargazer(m3_m_all, type = "text")
stargazer::stargazer(m3_m_t1, type = "text")
stargazer::stargazer(m3_m_t2, type = "text")


## Robust standard errors and P values
cov.mA <- vcovHC(m3_m_all, type="HC0")
std.errA <- sqrt(diag(cov.mA))
r.estA <- cbind(Estimate= coef(m3_m_all), "Robust SE" = std.errA,
                "Pr(>|z|)" = 2 * pnorm(abs(coef(m3_m_all)/std.errA), lower.tail=FALSE),
                LL = coef(m3_m_all) - 1.96 * std.errA,
                UL = coef(m3_m_all) + 1.96 * std.errA)

#openxlsx::write.xlsx(r.estA, paste0(dir, "Models/tA_m_mod.xlsx"), row.names = TRUE)


cov.m1 <- vcovHC(m3_m_t1, type="HC0")
std.err1 <- sqrt(diag(cov.m1))
r.est1 <- cbind(Estimate= coef(m3_m_t1), "Robust SE" = std.err1,
                "Pr(>|z|)" = 2 * pnorm(abs(coef(m3_m_t1)/std.err1), lower.tail=FALSE),
                LL = coef(m3_m_t1) - 1.96 * std.err1,
                UL = coef(m3_m_t1) + 1.96 * std.err1)

#openxlsx::write.xlsx(r.est1, paste0(dir, "Models/t1_m_mod.xlsx"), row.names = TRUE)

cov.m2 <- vcovHC(m3_m_t2, type="HC0")
std.err2 <- sqrt(diag(cov.m2))
r.est2 <- cbind(Estimate= coef(m3_m_t2), 
                "Robust SE" = std.err2,
                "Pr(>|z|)" = 2 * pnorm(abs(coef(m3_m_t2)/std.err2), lower.tail=FALSE),
                LL = coef(m3_m_t2) - 1.96 * std.err2,
                UL = coef(m3_m_t2) + 1.96 * std.err2)

#openxlsx::write.xlsx(r.est2, paste0(dir, "Models/t2_m_mod.xlsx"), row.names = TRUE)




### cfs for incident rate ratios
s1 <- deltamethod(list(~ exp(x1), ~ exp(x2), ~ exp(x3), ~ exp(x4), ~ exp(x5), ~ exp(x6), ~ exp(x7)), 
                  coef(m3_m_t1), cov.m1)

rexp.est1 <- exp(r.est1[, -3])
rexp.est1[, "Robust SE"] <- s1

rexp.est1


s2 <- deltamethod(list(~ exp(x1), ~ exp(x2), ~ exp(x3), ~ exp(x4), ~ exp(x5), ~ exp(x6), ~ exp(x7)), 
                  coef(m3_m_t2), cov.m2)

rexp.est2 <- exp(r.est2[, -3])
rexp.est2[, "Robust SE"] <- s2


rexp.est2


sA <- deltamethod(list(~ exp(x1), ~ exp(x2), ~ exp(x3), ~ exp(x4), ~ exp(x5), ~ exp(x6), ~ exp(x7)),  
                  coef(m3_m_all), cov.mA)

rexp.estA <- exp(r.estA[, -3])
rexp.estA[, "Robust SE"] <- sA

rexp.estA



### Descriptive table #######
out <- d3_m[,c("political_region","state","wealth_index_poorest_per","insured_per","median_age","urban_residence_per","period_deaths")]
#write.csv(out, paste0(dir, "Manuscript/table_S1.csv"), row.names = FALSE)



