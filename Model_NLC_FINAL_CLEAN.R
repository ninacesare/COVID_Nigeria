########################################################################################
# PROJECT       :	Leveraging AI
# SPONSOR/PI    : Elaine Nsoesie
# PROGRAM NAME  : 
# DESCRIPTION   : Modeling COVID-19 mortality across states in Nigeria
#                 This file creates analytic datasets, generates figure 1, and calculates descriptive statistics
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


### Load original dataset for states and political regions, which are referenced in death count and sociodemographic variable sections
d <- read_excel(paste0(dir, "updated_nigeria_risk_factor_june_24.xlsx"), sheet = 1)
## linking political region to state (from https://en.wikipedia.org/wiki/Geopolitical_zones_of_Nigeria)

d$political_region <- NA

d$political_region[which(d$state == "Benue")] <- "North Central"
d$political_region[which(d$state == "Kogi")] <- "North Central"
d$political_region[which(d$state == "Kwara")] <- "North Central"
d$political_region[which(d$state == "Nasarawa")] <- "North Central"
d$political_region[which(d$state == "Niger")] <- "North Central"
d$political_region[which(d$state == "Plateau")] <- "North Central"
d$political_region[which(d$state == "FCT Abuja")] <- "North Central"

d$political_region[which(d$state == "Adamawa")] <- "North East"
d$political_region[which(d$state == "Bauchi")] <- "North East"
d$political_region[which(d$state == "Borno")] <- "North East"
d$political_region[which(d$state == "Gombe")] <- "North East"
d$political_region[which(d$state == "Taraba")] <- "North East"
d$political_region[which(d$state == "Yobe")] <- "North East"

d$political_region[which(d$state == "Jigawa")] <- "North West"
d$political_region[which(d$state == "Kaduna")] <- "North West"
d$political_region[which(d$state == "Kano")] <- "North West"
d$political_region[which(d$state == "Katsina")] <- "North West"
d$political_region[which(d$state == "Kebbi")] <- "North West"
d$political_region[which(d$state == "Sokoto")] <- "North West"
d$political_region[which(d$state == "Zamfara")] <- "North West"

d$political_region[which(d$state == "Abia")] <- "South East"
d$political_region[which(d$state == "Anambra")] <- "South East"
d$political_region[which(d$state == "Ebonyi")] <- "South East"
d$political_region[which(d$state == "Enugu")] <- "South East"
d$political_region[which(d$state == "Imo")] <- "South East"

d$political_region[which(d$state == "Akwa Ibom")] <- "South South"
d$political_region[which(d$state == "Ibom")] <- "South South"
d$political_region[which(d$state == "Bayelsa")] <- "South South"
d$political_region[which(d$state == "Cross River")] <- "South South"
d$political_region[which(d$state == "Delta")] <- "South South"
d$political_region[which(d$state == "Edo")] <- "South South"
d$political_region[which(d$state == "Rivers")] <- "South South"

d$political_region[which(d$state == "Ekiti")] <- "South West"
d$political_region[which(d$state == "Lagos")] <- "South West"
d$political_region[which(d$state == "Ogun")] <- "South West"
d$political_region[which(d$state == "Ondo")] <- "South West"
d$political_region[which(d$state == "Osun")] <- "South West"
d$political_region[which(d$state == "Oyo")] <- "South West"




#### PROCESSING COVID DATA ####

#### Plotting deaths and cases over time #########
pop <- unique(d[,c("state","political_region","population_2016")])
names(pop) <- c("state","political_region","Population")

time_dat <- read.csv("//ad.bu.edu/bumcfiles/SPH/DCC/Dept/LeveragingAI/Spring2022/Nigeria COVID/DSN/covid_19_daily_data.csv") ## using this because it's the most recent daily data we have

time_dat$region_name <- as.character(time_dat$region_name)
time_dat$region_name[which(time_dat$region_name == "Federal Capital Territory")] <- "FCT Abuja"
names(time_dat)[grep("region", colnames(time_dat))] <- "state"

time_dat_merge <- merge(time_dat, pop, by = "state", all =TRUE) # Merge pop and political region with COVID daily data

time_dat_merge$date_new <- as.Date(as.character(time_dat_merge$date), "%Y-%m-%d")
time_dat_merge <- time_dat_merge[which(!is.na(time_dat_merge)),]
time_dat_merge <- time_dat_merge[which(time_dat_merge$state != "Non spécifié"),]

time_dat_merge$state <- as.factor(time_dat_merge$state)
time_dat_merge$per.100000 <- c(time_dat_merge$affected_infected/time_dat_merge$Population) * 100000  
time_dat_merge$per.100000_deaths <- c(time_dat_merge$affected_killed/time_dat_merge$Population) * 100000

time_dat_merge$week <- lubridate::week(time_dat_merge$date_new)
time_dat_merge$year <- unlist(lapply(time_dat_merge$date, function(x) unlist(strsplit(as.character(x), "-"))[1]))

time_dat_merge$week2 <- time_dat_merge$week
time_dat_merge$week2[which(time_dat_merge$year == 2021)] <- time_dat_merge$week[which(time_dat_merge$year == 2021)] + 53

time_dat_merge <- time_dat_merge[order(time_dat_merge$date_new),]

time_dat_merge <- time_dat_merge %>% group_by(week2) %>% mutate(date_new2 = first(date_new))

time_dat_merge_collapse <- time_dat_merge %>% dplyr::group_by(state, date_new2) %>% dplyr::summarize(affected_killed = sum(affected_killed, na.rm =TRUE), affected_infected = sum(affected_infected, na.rm =TRUE))
time_dat_merge_collapse <- merge(time_dat_merge_collapse, pop, by = "state")

time_dat_merge_collapse$per.100000 <- time_dat_merge_collapse$affected_infected/time_dat_merge_collapse$Population * 100000
time_dat_merge_collapse$per.100000_deaths<- time_dat_merge_collapse$affected_killed/time_dat_merge_collapse$Population * 100000

time_dat_merge_collapse$date_new2 <- as.Date(time_dat_merge_collapse$date_new2, by = "%Y-%m-%d")



#### break this into individual plots so we can see country lines
p1 <- ggplot(time_dat_merge_collapse[which(time_dat_merge_collapse$political_region == "North Central"),], aes(date_new2, per.100000_deaths, fill = state, color = state)) + 
  geom_line() + 
  xlab("Date") +
  ylab("Deaths per 100,000 residents") +
  ggtitle("North Central") +
  facet_wrap(.~political_region, ncol = 3) + 
  #theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


p2 <- ggplot(time_dat_merge_collapse[which(time_dat_merge_collapse$political_region == "North East"),], aes(date_new2, per.100000_deaths, fill = state, color = state)) + 
  geom_line() + 
  xlab("Date") +
  ylab("Deaths per 100,000 residents") +
  ggtitle("North East") +
  facet_wrap(.~political_region, ncol = 3) + 
  #theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


p3 <- ggplot(time_dat_merge_collapse[which(time_dat_merge_collapse$political_region == "North West"),], aes(date_new2, per.100000_deaths, fill = state, color = state)) + 
  geom_line() + 
  xlab("Date") +
  ylab("Deaths per 100,000 residents") +
  ggtitle("North West") +
  facet_wrap(.~political_region, ncol = 3) + 
  #theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


p4 <- ggplot(time_dat_merge_collapse[which(time_dat_merge_collapse$political_region == "South East"),], aes(date_new2, per.100000_deaths, fill = state, color = state)) + 
  geom_line() + 
  xlab("Date") +
  ylab("Deaths per 100,000 residents") +
  ggtitle("South East") +
  facet_wrap(.~political_region, ncol = 3) + 
  #theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


p5 <- ggplot(time_dat_merge_collapse[which(time_dat_merge_collapse$political_region == "South South"),], aes(date_new2, per.100000_deaths, fill = state, color = state)) + 
  geom_line() + 
  xlab("Date") +
  ylab("Deaths per 100,000 residents") +
  ggtitle("South South") +
  facet_wrap(.~political_region, ncol = 3) + 
  #theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


p6 <- ggplot(time_dat_merge_collapse[which(time_dat_merge_collapse$political_region == "South West"),], aes(date_new2, per.100000_deaths, fill = state, color = state)) + 
  geom_line() + 
  xlab("Date") +
  ylab("Deaths per 100,000 residents") +
  ggtitle("South West") +
  facet_wrap(.~political_region, ncol = 3) + 
  #theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#pdf("//ad.bu.edu/bumcfiles/SPH/DCC/Dept/LeveragingAI/Spring2022/Nigeria COVID/politicalRegion_deathsPer100K_time_long.pdf", width=10, height=12)
Rmisc::multiplot(p1, p2, p3, p4, p5, p6, cols = 2)
#dev.off()



#### break this into individual plots so we can see country lines
p1 <- ggplot(time_dat_merge_collapse[which(time_dat_merge_collapse$political_region == "North Central"),], aes(date_new2, affected_killed, fill = state, color = state)) + 
  geom_line() + 
  xlab("Date") +
  ylab("COVID-19 Deaths") +
  #ggtitle("North Central") +
  facet_wrap(.~political_region, ncol = 3) + 
  #theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme_light()


p2 <- ggplot(time_dat_merge_collapse[which(time_dat_merge_collapse$political_region == "North East"),], aes(date_new2, affected_killed, fill = state, color = state)) + 
  geom_line() + 
  xlab("Date") +
  ylab("COVID-19 Deaths") +
  #ggtitle("North East") +
  facet_wrap(.~political_region, ncol = 3) + 
  #theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme_light()


p3 <- ggplot(time_dat_merge_collapse[which(time_dat_merge_collapse$political_region == "North West"),], aes(date_new2, affected_killed, fill = state, color = state)) + 
  geom_line() + 
  xlab("Date") +
  ylab("COVID-19 Deaths") +
  #ggtitle("North West") +
  facet_wrap(.~political_region, ncol = 3) + 
  #theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme_light()


p4 <- ggplot(time_dat_merge_collapse[which(time_dat_merge_collapse$political_region == "South East"),], aes(date_new2, affected_killed, fill = state, color = state)) + 
  geom_line() + 
  xlab("Date") +
  ylab("COVID-19 Deaths") +
  #ggtitle("South East") +
  facet_wrap(.~political_region, ncol = 3) + 
  #theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme_light()


p5 <- ggplot(time_dat_merge_collapse[which(time_dat_merge_collapse$political_region == "South South"),], aes(date_new2, affected_killed, fill = state, color = state)) + 
  geom_line() + 
  xlab("Date") +
  ylab("COVID-19 Deaths") +
  #ggtitle("South South") +
  facet_wrap(.~political_region, ncol = 3) + 
  #theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme_light()


p6 <- ggplot(time_dat_merge_collapse[which(time_dat_merge_collapse$political_region == "South West"),], aes(date_new2, affected_killed, fill = state, color = state)) + 
  geom_line() + 
  xlab("Date") +
  ylab("COVID-19 Deaths") +
  #ggtitle("South West") +
  facet_wrap(.~political_region, ncol = 3) + 
  #theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme_light()

pdf("//ad.bu.edu/bumcfiles/SPH/DCC/Dept/LeveragingAI/Spring2022/Nigeria COVID/Figure_1.pdf", width=10, height=12)
Rmisc::multiplot(p1, p2, p3, p4, p5, p6, cols = 2)
dev.off()




#### For this analysis, we'll break the time data into two sections #### 
#https://gh.bmj.com/content/6/11/e007076

time_dat_merge1 <- time_dat_merge[which(time_dat_merge$date_new < "2020-10-25"),]
time_dat_merge2 <- time_dat_merge[which(time_dat_merge$date_new >= "2020-10-25"),]



#### Export analytic data for COVID infections and deaths by state ####
time_dat_mergeA_agg <- time_dat_merge %>% group_by(state) %>% summarize(political_region = political_region[1],
                                                                        period_deaths = sum(affected_killed, na.rm =TRUE), 
                                                                        period_cases = sum(affected_infected, na.rm =TRUE))

time_dat_merge1_agg <- time_dat_merge1 %>% group_by(state) %>% summarize(political_region = political_region[1], 
                                                                         period_deaths = sum(affected_killed, na.rm =TRUE), 
                                                                         period_cases = sum(affected_infected, na.rm =TRUE))

time_dat_merge2_agg <- time_dat_merge2 %>% group_by(state) %>% summarize(political_region = political_region[1], 
                                                                         period_deaths = sum(affected_killed, na.rm =TRUE), 
                                                                         period_cases = sum(affected_infected, na.rm =TRUE))

#openxlsx::write.xlsx(time_dat_mergeA_agg, "//ad.bu.edu/bumcfiles/SPH/DCC/Dept/LeveragingAI/Spring2022/Nigeria COVID/analytic_data/time_dat_mergeA_agg.xlsx", row.names = FALSE)
#openxlsx::write.xlsx(time_dat_merge1_agg, "//ad.bu.edu/bumcfiles/SPH/DCC/Dept/LeveragingAI/Spring2022/Nigeria COVID/analytic_data/time_dat_merge1_agg.xlsx", row.names = FALSE)
#openxlsx::write.xlsx(time_dat_merge2_agg, "//ad.bu.edu/bumcfiles/SPH/DCC/Dept/LeveragingAI/Spring2022/Nigeria COVID/analytic_data/time_dat_merge2_agg.xlsx", row.names = FALSE)



##### PROCESSING MALE AND FEMALE DEMOGRAPHIC DATA ######## 

#### Load and clean data ####

dat1 <- data.table::fread("//ad.bu.edu/bumcfiles/SPH/DCC/Dept/LeveragingAI/Spring2022/Nigeria COVID/NDHS/new_NGIR7AFL.csv") ## female data 
dat2 <- data.table::fread("//ad.bu.edu/bumcfiles/SPH/DCC/Dept/LeveragingAI/Spring2022/Nigeria COVID/NDHS/new_NGMR7AFL.csv") ## male data

names(dat1) <- gsub(" ", ".", colnames(dat1))
names(dat2) <- gsub(" ", ".", colnames(dat2))

names(dat1) <- tolower(colnames(dat1))
names(dat2) <- tolower(colnames(dat2))


summary(dat1$`respondent's.current.age`)
summary(dat2$current.age)

dat1 <- dat1 %>% rename(state = state.of.residence,
                        current.age = "respondent's.current.age",
                        type.of.place.of.residence = type.of.place.of.residence.b)

dat2 <- dat2 %>% rename(education.in.single.years = total.number.of.years.of.education)


#### Combine male and female data ####


dat3 <- rbind(dat1[,c("state",
                      "wealth.index.combined",
                      "current.age",
                      "covered.by.health.insurance",
                      "type.of.place.of.residence",
                      "frequency.of.reading.newspaper.or.magazine",
                      "frequency.of.listening.to.radio",
                      "frequency.of.watching.television",
                      "education.in.single.years")], 
              dat2[,c("state",
                      "wealth.index.combined",
                      "current.age",
                      "covered.by.health.insurance",
                      "type.of.place.of.residence",
                      "frequency.of.reading.newspaper.or.magazine",
                      "frequency.of.listening.to.radio",
                      "frequency.of.watching.television",
                      "education.in.single.years")])


#### Create measures within combined dataset ####
dat3$count <- 1


## Check for normality in age by state
dat3$state <- as.character(dat3$state)

par(mfrow = c(1,2))
plot(density(dat3$current.age, na.rm = TRUE))
qqnorm(dat3$current.age, main='Normal')
qqline(dat3$current.age)


shapiro.test(dat3$current.age[sample(1:nrow(dat3), 5000)])

normality <- dat3 %>% group_by(state) %>% summarise(shapiro.test(current.age)$p.value)



dat3_agg <- dat3 %>% dplyr::group_by(state) %>% dplyr::summarize(no_of_participants = sum(count), 
                                                                 wealth_index_poorest = length(which(wealth.index.combined == "poorest")),
                                                                 median_age = median(current.age, na.rm =TRUE),
                                                                 mean_age = mean(current.age, na.rm = TRUE),
                                                                 insured = length(which(covered.by.health.insurance == "yes")),
                                                                 urban_residence = length(which(type.of.place.of.residence == "urban")),
                                                                 paper_at_least_once_week = length(which(frequency.of.reading.newspaper.or.magazine == "at least once a week")),
                                                                 radio_at_least_once_week = length(which(frequency.of.listening.to.radio == "at least once a week")),
                                                                 tv_at_least_once_week = length(which(frequency.of.watching.television == "at least once a week")),
                                                                 avg_education = mean(education.in.single.years, na.rm = TRUE))


dat3_agg$wealth_index_poorest_per <- dat3_agg$wealth_index_poorest/dat3_agg$no_of_participants * 100
dat3_agg$insured_per <- dat3_agg$insured/dat3_agg$no_of_participants * 100
dat3_agg$urban_residence_per <- dat3_agg$urban_residence/dat3_agg$no_of_participants * 100
dat3_agg$paper_at_least_once_week_pct <- dat3_agg$paper_at_least_once_week/dat3_agg$no_of_participants * 100
dat3_agg$radio_at_least_once_week_pct <- dat3_agg$radio_at_least_once_week/dat3_agg$no_of_participants * 100
dat3_agg$tv_at_least_once_week_pct <- dat3_agg$tv_at_least_once_week/dat3_agg$no_of_participants * 100


d$state_lower <- tolower(d$state)
dat3_agg$state_lower <- dat3_agg$state

dat3_agg <- merge(dat3_agg[,-1], d[,c("political_region","state","state_lower","cases_confirmed_in_lab","population_2016")], by = "state_lower")

dat3_agg$population_2016_t <- dat3_agg$population_2016/10000


#### Export analytic set

openxlsx::write.xlsx(dat3_agg, "//ad.bu.edu/bumcfiles/SPH/DCC/Dept/LeveragingAI/Spring2022/Nigeria COVID/analytic_data/dat3_agg.xlsx")


############ Calculate variables and export for male, female dataset for robustness checks (dat1 = female; dat2 = male) ###############

#### For female (dat1) ####

dat1 <- dat1[,c("state",
               "wealth.index.combined",
               "current.age",
               "covered.by.health.insurance",
               "type.of.place.of.residence",
               "frequency.of.reading.newspaper.or.magazine",
               "frequency.of.listening.to.radio",
               "frequency.of.watching.television",
               "education.in.single.years")]             

dat1$count <- 1


## Check for normality in age by state
dat1$state <- as.character(dat1$state)

par(mfrow = c(1,2))
plot(density(dat1$current.age, na.rm = TRUE))
qqnorm(dat1$current.age, main='Normal')
qqline(dat1$current.age)


shapiro.test(dat1$current.age[sample(1:nrow(dat1), 5000)])

normality <- dat1 %>% group_by(state) %>% summarise(shapiro.test(current.age)$p.value)



dat1_agg <- dat1 %>% dplyr::group_by(state) %>% dplyr::summarize(no_of_participants = sum(count), 
                                                                 wealth_index_poorest = length(which(wealth.index.combined == "poorest")),
                                                                 median_age = median(current.age, na.rm =TRUE),
                                                                 mean_age = mean(current.age, na.rm = TRUE),
                                                                 insured = length(which(covered.by.health.insurance == "yes")),
                                                                 urban_residence = length(which(type.of.place.of.residence == "urban")),
                                                                 paper_at_least_once_week = length(which(frequency.of.reading.newspaper.or.magazine == "at least once a week")),
                                                                 radio_at_least_once_week = length(which(frequency.of.listening.to.radio == "at least once a week")),
                                                                 tv_at_least_once_week = length(which(frequency.of.watching.television == "at least once a week")),
                                                                 avg_education = mean(education.in.single.years, na.rm = TRUE))


dat1_agg$wealth_index_poorest_per <- dat1_agg$wealth_index_poorest/dat1_agg$no_of_participants * 100
dat1_agg$insured_per <- dat1_agg$insured/dat1_agg$no_of_participants * 100
dat1_agg$urban_residence_per <- dat1_agg$urban_residence/dat1_agg$no_of_participants * 100
dat1_agg$paper_at_least_once_week_pct <- dat1_agg$paper_at_least_once_week/dat1_agg$no_of_participants * 100
dat1_agg$radio_at_least_once_week_pct <- dat1_agg$radio_at_least_once_week/dat1_agg$no_of_participants * 100
dat1_agg$tv_at_least_once_week_pct <- dat1_agg$tv_at_least_once_week/dat1_agg$no_of_participants * 100


d$state_lower <- tolower(d$state)
dat1_agg$state_lower <- dat1_agg$state

dat1_agg <- merge(dat1_agg[,-1], d[,c("political_region","state","state_lower","cases_confirmed_in_lab","population_2016")], by = "state_lower")

dat1_agg$population_2016_t <- dat1_agg$population_2016/10000


#### Export analytic set

openxlsx::write.xlsx(dat1_agg, "//ad.bu.edu/bumcfiles/SPH/DCC/Dept/LeveragingAI/Spring2022/Nigeria COVID/analytic_data/dat1_agg.xlsx")


#### For male (dat2) ####


dat2 <- dat2[,c("state",
                "wealth.index.combined",
                "current.age",
                "covered.by.health.insurance",
                "type.of.place.of.residence",
                "frequency.of.reading.newspaper.or.magazine",
                "frequency.of.listening.to.radio",
                "frequency.of.watching.television",
                "education.in.single.years")]   

dat2$count <- 1


## Check for normality in age by state
dat2$state <- as.character(dat2$state)

par(mfrow = c(1,2))
plot(density(dat2$current.age, na.rm = TRUE))
qqnorm(dat2$current.age, main='Normal')
qqline(dat2$current.age)


shapiro.test(dat2$current.age[sample(1:nrow(dat2), 5000)])

normality <- dat2 %>% group_by(state) %>% summarise(shapiro.test(current.age)$p.value)



dat2_agg <- dat2 %>% dplyr::group_by(state) %>% dplyr::summarize(no_of_participants = sum(count), 
                                                                 wealth_index_poorest = length(which(wealth.index.combined == "poorest")),
                                                                 median_age = median(current.age, na.rm =TRUE),
                                                                 mean_age = mean(current.age, na.rm = TRUE),
                                                                 insured = length(which(covered.by.health.insurance == "yes")),
                                                                 urban_residence = length(which(type.of.place.of.residence == "urban")),
                                                                 paper_at_least_once_week = length(which(frequency.of.reading.newspaper.or.magazine == "at least once a week")),
                                                                 radio_at_least_once_week = length(which(frequency.of.listening.to.radio == "at least once a week")),
                                                                 tv_at_least_once_week = length(which(frequency.of.watching.television == "at least once a week")),
                                                                 avg_education = mean(education.in.single.years, na.rm = TRUE))


dat2_agg$wealth_index_poorest_per <- dat2_agg$wealth_index_poorest/dat2_agg$no_of_participants * 100
dat2_agg$insured_per <- dat2_agg$insured/dat2_agg$no_of_participants * 100
dat2_agg$urban_residence_per <- dat2_agg$urban_residence/dat2_agg$no_of_participants * 100
dat2_agg$paper_at_least_once_week_pct <- dat2_agg$paper_at_least_once_week/dat2_agg$no_of_participants * 100
dat2_agg$radio_at_least_once_week_pct <- dat2_agg$radio_at_least_once_week/dat2_agg$no_of_participants * 100
dat2_agg$tv_at_least_once_week_pct <- dat2_agg$tv_at_least_once_week/dat2_agg$no_of_participants * 100


d$state_lower <- tolower(d$state)
dat2_agg$state_lower <- dat2_agg$state

dat2_agg <- merge(dat2_agg[,-1], d[,c("political_region","state","state_lower","cases_confirmed_in_lab","population_2016")], by = "state_lower")

dat2_agg$population_2016_t <- dat2_agg$population_2016/10000


#### Export analytic set 

openxlsx::write.xlsx(dat2_agg, "//ad.bu.edu/bumcfiles/SPH/DCC/Dept/LeveragingAI/Spring2022/Nigeria COVID/analytic_data/dat2_agg.xlsx")

