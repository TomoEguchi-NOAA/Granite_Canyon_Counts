#SplineJags_AllYears.R
#Runs the spline-only model on all years. Jags data are combined between 
#LaakeData2SplineJags.R and SplineJags_2006to2023.R
#

rm(list = ls())
library(tidyverse)
library(jagsUI)
library(ggplot2)
library(loo)
library(bayesplot)
library(ERAnalysis)

jags.model <- paste0("models/model_Nmix_Spline_v2_JAGS.txt")


data(ERSurveyData)
data("Observer")

# The data in PrimarySightings are all southbound sightings for all years in which visibility and beaufort
# are less than or equal to 4. Below the counts are shown for the 2 dataframes for
# recent surveys since 1987/88.
data(PrimaryEffort)

# Filter effort and sightings and store in dataframes Effort and Sightings
Effort.1 = PrimaryEffort[PrimaryEffort$Use,]  

data("SecondaryEffort")
Effort.2 <- SecondaryEffort[SecondaryEffort$Use,]

# For jags and WinBugs code, what I need are
# 1. observed number of whales per day n[d, s, y], where d = # days since 12/1,
# s = station (1 = primary, 2 = secondary), y = year. For Laake's data, maximum 
# number of days per year was 94. 
# 2. Beaufort sea state bf[d,y]
# 3. Visibility code vs[d,y]
# 4. observer code obs[d,s,y]
# 5. the proportion of watch duration per day (out of 540 minutes or 9 hrs) watch.prop[d,y]
# 6. index of survey day, i.e., the number of days since 12/1 day[d,y]

# Need to count the number of days since 12-01 for each year. But 12-01 is 1.
# Then... 
# Count the number of whales per day and daily effort
# In early years, surveys were conducted 10 hrs. So, the watch proportion
# can be > 1.0, because we have used 9 hrs as maximum. 

# Summarizing by day worked fine but the model requires counts per observation
# period. Needs to be redone. 2023-09-15 DONE.

Effort.1 %>% 
  mutate(Day1 = as.Date(paste0(Start.year, "-11-30")),
         dt = as.numeric(as.Date(Date) - Day1),
         obs = Observer) %>%
  dplyr::select(Start.year, nwhales, effort, vis, beaufort, obs, dt) %>%
  group_by(Start.year) %>%
  mutate(effort.min = effort * 24 * 60,
         watch.prop = effort.min/540) -> Effort.by.period.1

Effort.2 %>% 
  mutate(Day1 = as.Date(paste0(Start.year, "-11-30")),
         dt = as.numeric(as.Date(Date) - Day1) + 1,
         obs = Observer) %>%
  dplyr::select(Start.year, nwhales, effort, vis, beaufort, obs, dt) %>%
  group_by(Start.year) %>%
  mutate(effort.min = effort * 24 * 60,
         watch.prop = effort.min/540) -> Effort.by.period.2

# Create matrices that contain necessary data
all.years <- unique(Effort.1$Start.year)
n.year <- length(all.years)

# Primary station
Effort.by.period.1 %>% 
  summarize(n = n()) -> periods.1

n.1 <- matrix(data = 0, nrow = max(periods.1$n)+2, ncol = n.year)
Watch.Length.1 <- matrix(data = 0, nrow = max(periods.1$n)+2, ncol = n.year)
day.1 <- matrix(nrow = max(periods.1$n)+2, ncol = n.year)
bf.1 <- vs.1 <- matrix(data = 0, nrow = max(periods.1$n), ncol = n.year)
obs.1 <-  matrix(nrow = max(periods.1$n), ncol = n.year)
u.1 <- matrix(data = 0, nrow = max(periods.1$n)+2, ncol = n.year)

for (k in 1:n.year){

  n.1[1:periods.1$n[k], k] <- Effort.by.period.1 %>% 
    ungroup() %>%
    filter(Start.year == all.years[k]) %>%
    dplyr::select(nwhales) %>%
    pull()
  n.1[(periods.1$n[k]+1):(periods.1$n[k]+2), k] <- 0
  
  Watch.Length.1[1:periods.1$n[k], k] <- Effort.by.period.1 %>% 
    ungroup() %>%
    filter(Start.year == all.years[k]) %>%
    dplyr::select(watch.prop) %>%
    pull()
  Watch.Length.1[(periods.1$n[k]+1):(periods.1$n[k]+2), k] <- 0
  
  bf.1[1:periods.1$n[k], k] <- Effort.by.period.1 %>% 
    ungroup() %>%
    filter(Start.year == all.years[k]) %>%
    dplyr::select(beaufort) %>%
    pull()
  
  vs.1[1:periods.1$n[k], k] <- Effort.by.period.1 %>% 
    ungroup() %>%
    filter(Start.year == all.years[k]) %>%
    dplyr::select(vis) %>%
    pull()
  
  obs.1[1:periods.1$n[k], k] <- Effort.by.period.1 %>% 
    ungroup() %>%
    filter(Start.year == all.years[k]) %>%
    dplyr::select(obs) %>%
    pull()
  
  u.1[1:periods.1$n[k], k] <- 1
  
  day.1[1:periods.1$n[k], k] <- Effort.by.period.1 %>% 
    ungroup() %>%
    filter(Start.year == all.years[k]) %>%
    dplyr::select(dt) %>%
    pull()
  day.1[(periods.1$n[k]+1):(periods.1$n[k]+2), k] <- c(1,90)
  
}

# Secondary station
Effort.by.period.2 %>% 
  summarize(n = n()) -> periods.2

n.2 <- matrix(data = 0, nrow = max(periods.2$n)+2, ncol = n.year)
Watch.Length.2 <- matrix(data = 0, nrow = max(periods.2$n)+2, ncol = n.year)
day.2 <- matrix(data = 1, nrow = max(periods.2$n)+2, ncol = n.year)
bf.2 <- vs.2 <- matrix(data = 0, nrow = max(periods.2$n), ncol = n.year)
obs.2 <- matrix(data = 36, nrow = max(periods.2$n), ncol = n.year)
u.2 <- matrix(data = 0, nrow = max(periods.2$n)+2, ncol = n.year)

c <- 1
for (k in 1:n.year){
  
  if (all.years[k] %in% periods.2$Start.year){
    n.2[1:periods.2$n[c], c] <- Effort.by.period.2 %>% 
      ungroup() %>%
      filter(Start.year == all.years[k]) %>%
      dplyr::select(nwhales) %>%
      pull()
    n.2[(periods.2$n[c]+1):(periods.2$n[c]+2), k] <- 0
    
    Watch.Length.2[1:periods.2$n[c], c] <- Effort.by.period.2 %>% 
      ungroup() %>%
      filter(Start.year == all.years[k]) %>%
      dplyr::select(watch.prop) %>%
      pull()
    Watch.Length.1[(periods.2$n[c]+1):(periods.2$n[c]+2), k] <- 0
    
    bf.2[1:periods.2$n[c], c] <- Effort.by.period.2 %>% 
      ungroup() %>%
      filter(Start.year == all.years[k]) %>%
      dplyr::select(beaufort) %>%
      pull()
    
    vs.2[1:periods.2$n[c], c] <- Effort.by.period.2 %>% 
      ungroup() %>%
      filter(Start.year == all.years[k]) %>%
      dplyr::select(vis) %>%
      pull()
    
    obs.2[1:periods.2$n[c], c] <- Effort.by.period.2 %>% 
      ungroup() %>%
      filter(Start.year == all.years[k]) %>%
      dplyr::select(obs) %>%
      pull()
    
    u.2[1:periods.2$n[c], c] <- 1
    
    day.2[1:periods.2$n[c], c] <- Effort.by.period.2 %>% 
      ungroup() %>%
      filter(Start.year == all.years[k]) %>%
      dplyr::select(dt) %>%
      pull()
    day.2[(periods.2$n[c]+1):(periods.2$n[c]+2), c] <- c(1,90)
    
    c <- c + 1
  }
}

periods.1 %>%
  left_join(periods.2, by = "Start.year") -> all.periods

all.periods[is.na(all.periods)] <- 0

# combine input data for Laake's and Durban's approaches
# 
# More recent data
# Bring in the output from a successful run of WinBUGS (a rare event...)
# BUGS input and output
x <- 10
min.duration <- 85
in.file.name <- paste0("RData/WinBUGS_", x, "yr_v2_min", min.duration, ".rds")
BUGS.results <- readRDS(in.file.name)
BUGS.data <- BUGS.results$BUGS.data

# The first year is a duplicate from Laake's so remove.
n.1 <- cbind(n.1, 
             rbind(as.matrix(BUGS.data$n[,1,2:x]), 
                   matrix(data = 0, 
                          nrow = nrow(n.1) - dim(BUGS.data$n)[1], 
                          ncol = dim(BUGS.data$n)[3]-1)))

n.2 <- cbind(n.2, 
             rbind(as.matrix(BUGS.data$n[,2,2:x]), 
                   matrix(data = 0, 
                          nrow = nrow(n.2) - dim(BUGS.data$n)[1], 
                          ncol = dim(BUGS.data$n)[3]-1)))

Watch.Length.1 <- cbind(Watch.Length.1,
                        rbind(BUGS.data$Watch.Length[, 2:x], 
                              matrix(data = NA, 
                                     nrow = nrow(Watch.Length.1) - dim(BUGS.data$Watch.Length)[1], 
                                     ncol = dim(BUGS.data$Watch.Length)[2]-1)))

Watch.Length.2 <- cbind(Watch.Length.2,
                        rbind(BUGS.data$Watch.Length[, 2:x], 
                              matrix(data = NA, 
                                     nrow = nrow(Watch.Length.2) - dim(BUGS.data$Watch.Length)[1], 
                                     ncol = dim(BUGS.data$Watch.Length)[2]-1)))
bf.1 <- cbind(bf.1,
              rbind(BUGS.data$bf[, 2:x], 
                    matrix(data = NA, 
                           nrow = nrow(bf.1) - dim(BUGS.data$bf)[1], 
                           ncol = dim(BUGS.data$bf)[2]-1)))

vs.1 <- cbind(vs.1,
              rbind(BUGS.data$vs[, 2:x], 
                    matrix(data = NA, 
                           nrow = nrow(vs.1) - dim(BUGS.data$vs)[1], 
                           ncol = dim(BUGS.data$vs)[2]-1)))

bf.2 <- cbind(bf.2,
              rbind(BUGS.data$bf[, 2:x], 
                    matrix(data = NA, 
                           nrow = nrow(bf.2) - dim(BUGS.data$bf)[1], 
                           ncol = dim(BUGS.data$bf)[2]-1)))

vs.2 <- cbind(vs.2,
              rbind(BUGS.data$vs[, 2:x], 
                    matrix(data = NA, 
                           nrow = nrow(vs.2) - dim(BUGS.data$vs)[1], 
                           ncol = dim(BUGS.data$vs)[2]-1)))


obs.1 <- cbind(obs.1, rbind(as.matrix(BUGS.data$obs[,1,2:x]), 
                            matrix(data = 36, 
                                   nrow = nrow(obs.1) - dim(BUGS.data$obs)[1], 
                                   ncol = dim(BUGS.data$obs)[3]-1)))

obs.2 <- cbind(obs.2, rbind(as.matrix(BUGS.data$obs[,2,2:x]), 
                            matrix(data = 36, 
                                   nrow = nrow(obs.2) - dim(BUGS.data$obs)[1], 
                                   ncol = dim(BUGS.data$obs)[3]-1)))

day.1 <- cbind(day.1,
               rbind(BUGS.data$day[, 2:x], 
                     matrix(data = NA, 
                            nrow = nrow(day.1) - dim(BUGS.data$day)[1], 
                            ncol = dim(BUGS.data$day)[2]-1)))

day.2 <- cbind(day.2,
               rbind(BUGS.data$day[, 2:x], 
                     matrix(data = NA, 
                            nrow = nrow(day.2) - dim(BUGS.data$day)[1], 
                            ncol = dim(BUGS.data$day)[2]-1)))

u.1 <- cbind(u.1, 
             rbind(as.matrix(BUGS.data$u[,1,2:x]), 
                   matrix(data = 0, 
                          nrow = nrow(u.1) - dim(BUGS.data$u)[1], 
                          ncol = dim(BUGS.data$u)[3]-1)))

u.2 <- cbind(u.2, 
             rbind(as.matrix(BUGS.data$u[,1,2:x]), 
                   matrix(data = 0, 
                          nrow = nrow(u.2) - dim(BUGS.data$u)[1], 
                          ncol = dim(BUGS.data$u)[3]-1)))

periods.1 <- c(all.periods$n.x, BUGS.data$periods[2:length(BUGS.data$periods)])
periods.2 <- c(all.periods$n.y, 0, BUGS.data$periods[3:4], rep(0, x-3))

# Create jags input
N.1.inits <- matrix(nrow = nrow(n.1), ncol = ncol(n.1))
N.1.inits[day.1 == 1 | day.1 > 89] <- 0

N.2.inits <- matrix(nrow = nrow(n.2), ncol = ncol(n.2))
N.2.inits[day.2 == 1 | day.2 > 89] <- 0

jags.data <- list(n.1 = n.1,
                  n.2 = n.2,
                  N.1 = N.1.inits,
                  N.2 = N.2.inits,
                  n.station = ifelse(colSums(n.2) > 0, 1, 0),
                  n.year = ncol(n.1),
                  n.obs = max(obs.1, na.rm = T),
                  periods.1 = periods.1,
                  periods.2 = periods.2,
                  obs.1 = obs.1,
                  obs.2 = obs.2,
                  vs.1 = vs.1,
                  vs.2 = vs.2,
                  bf.1 = bf.1,
                  bf.2 = bf.2,
                  Watch.Length.1 = Watch.Length.1,
                  Watch.Length.2 = Watch.Length.2,
                  day.1 = day.1,
                  day.2 = day.2,
                  u.1 = u.1,
                  u.2 = u.2,
                  n.days = apply(day.1, MARGIN = 2, FUN = max, na.rm = T),
                  knot = c(-1.46,-1.26,-1.02,-0.78,
                           -0.58,-0.34,-0.10,0.10,
                           0.34,0.57,0.78,1.02,1.26,1.46),
                  n.knots = 14)

jags.params <- c("lambda.sp",
                 "OBS.RF.sp",
                 "OBS.Switch.sp",
                 "BF.Switch.sp",
                 "BF.Fixed.sp",
                 "VS.Switch.sp",
                 "VS.Fixed.sp",
                 "mean.prob.sp",
                 "BF.Fixed.sp",
                 "VS.Fixed.sp",
                 "Corrected.Est",
                 "Raw.Est",
                 "z",
                 "sp",
                 "Daily.Est",
                 "beta.sp",
                 "b.sp",
                 "sd.b.sp",
                 "log.lkhd")

MCMC.params <- list(n.samples = 250000,
                    n.thin = 100,
                    n.burnin = 200000,
                    n.chains = 5)

# Needs initial values for N to avoid errors. N has to be significantly larger
# than observed n.
jags.inits <- function() list(N.1 = jags.data$n.1[,] * 5,
                              N.2 = jags.data$n.2[,] * 5)

out.file.name <- "RData/JAGS_Spline_results_All_Data.rds"

if (!file.exists(out.file.name)){
  Start_Time<-Sys.time()
  
  jm <- jagsUI::jags(jags.data,
                     inits = function() list(N.1 = jags.data$n.1[,] * 5,
                                             N.2 = jags.data$n.2[,] * 5),
                     parameters.to.save= jags.params,
                     model.file = jags.model,
                     n.chains = MCMC.params$n.chains,
                     n.burnin = MCMC.params$n.burnin,
                     n.thin = MCMC.params$n.thin,
                     n.iter = MCMC.params$n.samples,
                     DIC = T, 
                     parallel=T)
  
  Run_Time <- Sys.time() - Start_Time
  jm.out <- list(jm = jm,
                 jags.data = jags.data,
                 jags.params = jags.params,
                 jags.model = jags.model,
                 MCMC.params = MCMC.params,
                 Run_Time = Run_Time,
                 System = Sys.getenv(),
                 Run_Date = Sys.Date())
  
  saveRDS(jm.out,
          file = out.file.name)
  
} else {
  jm.out <- readRDS(out.file.name)
}
# 
# seasons <- c("2006/2007", "2007/2008", "2009/2010", "2010/2011", 
#              "2014/2015", "2015/2016", "2019/2020", "2021/2022",
#              "2022/2023")
# 
