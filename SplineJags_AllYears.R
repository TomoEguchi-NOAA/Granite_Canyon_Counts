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

# I want to convert observers to integers as soon as possible
# create a new observer list from all analysis after Laake's
# output from Ver2.0 extraction
years <- c("2010", "2011", "2015", "2016", "2020", "2022", "2023", "2024")

# I don't have raw data for 2006/2007 and 2007/2008. Laake's data include
# 2006/2007. So, missing 2007/2008.

# Extracted data using v2 method (mine)
out.v2 <- lapply(years, 
                 FUN = function(x) readRDS(paste0("RData/V2.1_Feb2024/out_", x,
                                                  "_min85_Tomo_v2.rds")))
# extract unique observers
all.new.obs <- lapply(out.v2, FUN = function(x) x$Complete_Data %>% 
                        select(obs) %>% 
                        unique()) %>% 
  unlist() %>% 
  unique() %>% 
  unname()

#old.obs.list <- read.csv(file = "Data/Observer list.csv")
Laake.obs <- data.frame(Observer = unique(c(unique(Effort.1$Observer), 
                                            unique(Effort.2$Observer)))) %>%
  transmute(Observer.char = as.character(Observer))

# The last two rows are duplicates
Observer <- Observer[1:67,]

Laake.obs %>% left_join(Observer %>% 
                          mutate(Observer.char = as.character(Observer)), 
                        by = "Observer.char") %>%
  select(Observer.char, Observer) -> Laake.obs.1

Laake.obs.1$Observer[is.na(Laake.obs.1$Observer)] <- (max(Laake.obs.1$Observer, na.rm = T)+1) : 
  ((sum(is.na(Laake.obs.1$Observer))) + max(Laake.obs.1$Observer, na.rm = T))

# Combine Laake's observer list to the new one
all.obs <- unique(c(select(Laake.obs.1, Observer.char) %>% pull(), all.new.obs))

all.obs.df <- data.frame(Observer.char = all.obs)

Laake.obs.1 %>% right_join(all.obs.df, by = "Observer.char") %>%
  transmute(Observer.char = Observer.char,
            obs.int = Observer) -> updated.obs.list
updated.obs.list$obs.int[is.na(updated.obs.list$obs.int)] <- (max(updated.obs.list$obs.int, 
                                                                  na.rm = T)+1) : 
  ((sum(is.na(updated.obs.list$obs.int))) + max(updated.obs.list$obs.int, na.rm = T))


# Finally, I need to convert 2007/2008 observer data into the new observer code.
# More recent data
# Bring in the output from a successful run of WinBUGS 
# BUGS input and output
x <- 10

in.file.name <- paste0("RData/WinBUGS_", x, "yr_v2_min85.rds")
BUGS.results <- readRDS(in.file.name)
BUGS.data <- BUGS.results$BUGS.data

# observers in the winbugs data need to be all replaced with the new numbering system.
# So... I need to re-do the numbering. 

all.obs.WinBUGS <- lapply(out.v2, FUN = function(x){
  return(x$Final_Data %>%
           #filter(station == "P") %>%
           select(obs) %>%
           unique()) 
})

# I don't have raw data for 2007/2008 so I have to back convert using the observer list csv file
# I'm assuming that IDs in the observer list csv file match with what were found in WinBUGS
# input. first find all unique observer initials and create matching new IDs.
obs.old.list <- read.csv("Data/Observer list.csv") %>%
  filter(ID != 36) %>%   # 36 was used for no observer
  transmute(Observer.char = Observer,
            ID = ID)

obs.2008.1 <- data.frame(ID = unique(BUGS.data$obs[,1,2])) %>%
  left_join(obs.old.list, by = "ID") %>%
  filter(ID != 36) %>%   # 36 was used as no observer.
  pull(Observer.char)

all.obs.uniq <- unique(c(unique(unlist(all.obs.WinBUGS)), 
                         updated.obs.list$Observer.char,
                         obs.2008.1))

obs.df <- data.frame(Observer.char = all.obs.uniq,
                     obs.int = seq(1, length(all.obs.uniq)))

no.obs <- max(obs.df$obs.int) + 1

# Replace old IDs with new IDs for 2007/2008
obs.2008.df <- data.frame(ID = BUGS.data$obs[,1,2]) %>%
  filter(ID != 36) %>%
  left_join(obs.old.list, by = "ID") %>%
  left_join(obs.df, by = "Observer.char")

# Find all observers for the primary stations
all.Final_Data.1 <- lapply(out.v2, FUN = function(x){
  return(x$Final_Data %>%
           filter(station == "P") %>%
           mutate(Observer.char = obs) %>%
           left_join(obs.df, by = "Observer.char"))
})

# assign new IDs (obs.int)
obs.1.1.list <- lapply(all.Final_Data.1, FUN = function(x) x$obs.int)
max.length.1 <- lapply(obs.1.1.list, FUN = length) %>% unlist() 

obs.1.1 <- matrix(data = no.obs, nrow = max(max.length.1), ncol = length(obs.1.1.list))

for (k in 1:length(obs.1.1.list)){
  obs.1.1[1:max.length.1[k], k] <- obs.1.1.list[[k]]
}

# Add 2008 observer IDs (with new numbering system)
obs.1.1 <- cbind(c(obs.2008.df$obs.int, 
                   rep(no.obs, 
                       times = (nrow(obs.1.1) - (nrow(obs.2008.df))))), 
                 obs.1.1)

# Do the same with the secondary station
all.Final_Data.2 <- lapply(out.v2, FUN = function(x){
  return(x$Final_Data %>%
           filter(station == "S") %>%
           mutate(Observer.char = obs) %>%
           left_join(obs.df, by = "Observer.char"))
})

obs.2.1.list <- lapply(all.Final_Data.2, FUN = function(x) x$obs.int)
max.length.2 <- lapply(obs.2.1.list, FUN = length) %>% unlist() 

obs.2.1 <- matrix(data = no.obs, 
                  nrow = max(max.length.2), 
                  ncol = length(obs.2.1.list))

for (k in 1:length(obs.2.1.list)){
  if (length(obs.2.1.list[[k]]) > 0)
    obs.2.1[1:max.length.2[k], k] <- obs.2.1.list[[k]]
}

# Add 2008 observer IDs (with new numbering system)
# There was no secondary station for 2007/2008
obs.2.1 <- cbind( rep(no.obs, nrow(obs.2.1)), obs.2.1)

# For new jags code, what I need are
# 1. observed number of whales per sampling period n[t, y], where t = sequential period,
# n.1 and n.2 for primary and secondary, respectively. y = year. 
#  For Laake's data, maximum # number of days per year was 94. 
# 2. Beaufort sea state bf[t,y]
# 3. Visibility code vs[t,y]
# 4. observer code obs.1[t,y] and obs.2
# 5. the proportion of watch duration per day (out of 540 minutes or 9 hrs) watch.prop.1[t,y]
# and watch.prop.2[t,y]
# 6. index of survey day, i.e., the number of days since 12/1 day[t,y]

# Need to count the number of days since 12-01 for each year. But 12-01 is 1.
# Then... 
# Count the number of whales per day and daily effort
# In early years, surveys were conducted 10 hrs. So, the watch proportion
# can be > 1.0, because we have used 9 hrs as maximum. 

Effort.1 %>% 
  mutate(Day1 = as.Date(paste0(Start.year, "-11-30")),
         dt = as.numeric(as.Date(Date) - Day1),
         obs = Observer) %>%
  dplyr::select(Start.year, nwhales, effort, vis, beaufort, obs, dt) %>%
  group_by(Start.year) %>%
  mutate(effort.min = effort * 24 * 60,
         watch.prop = effort.min/540,
         Observer.char = as.character(obs)) %>%
  left_join(updated.obs.list, by = "Observer.char")-> Effort.by.period.1

Effort.2 %>% 
  mutate(Day1 = as.Date(paste0(Start.year, "-11-30")),
         dt = as.numeric(as.Date(Date) - Day1) + 1,
         obs = Observer) %>%
  dplyr::select(Start.year, nwhales, effort, vis, beaufort, obs, dt) %>%
  group_by(Start.year) %>%
  mutate(effort.min = effort * 24 * 60,
         watch.prop = effort.min/540,
         Observer.char = as.character(obs)) %>%
  left_join(updated.obs.list, by = "Observer.char")-> Effort.by.period.2

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
obs.1 <-  matrix(data = no.obs, nrow = max(periods.1$n), ncol = n.year)
u.1 <- matrix(data = 0, nrow = max(periods.1$n)+2, ncol = n.year)

k <- 17
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
    dplyr::select(obs.int) %>%
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
obs.2 <- matrix(data = no.obs, nrow = max(periods.2$n), ncol = n.year)
u.2 <- matrix(data = 0, nrow = max(periods.2$n)+2, ncol = n.year)

k <- c <- 1
for (k in 1:n.year){
  
  if (all.years[k] %in% periods.2$Start.year){
    n.2[1:periods.2$n[c], k] <- Effort.by.period.2 %>% 
      ungroup() %>%
      filter(Start.year == all.years[k]) %>%
      dplyr::select(nwhales) %>%
      pull()
    n.2[(periods.2$n[c]+1):(periods.2$n[c]+2), k] <- 0
    
    Watch.Length.2[1:periods.2$n[c], k] <- Effort.by.period.2 %>% 
      ungroup() %>%
      filter(Start.year == all.years[k]) %>%
      dplyr::select(watch.prop) %>%
      pull()
    Watch.Length.2[(periods.2$n[c]+1):(periods.2$n[c]+2), k] <- 0
    
    bf.2[1:periods.2$n[c], k] <- Effort.by.period.2 %>% 
      ungroup() %>%
      filter(Start.year == all.years[k]) %>%
      dplyr::select(beaufort) %>%
      pull()
    
    vs.2[1:periods.2$n[c], k] <- Effort.by.period.2 %>% 
      ungroup() %>%
      filter(Start.year == all.years[k]) %>%
      dplyr::select(vis) %>%
      pull()
    
    obs.2[1:periods.2$n[c], k] <- Effort.by.period.2 %>% 
      ungroup() %>%
      filter(Start.year == all.years[k]) %>%
      dplyr::select(obs.int) %>%
      pull()
    
    u.2[1:periods.2$n[c], k] <- 1
    
    day.2[1:periods.2$n[c], k] <- Effort.by.period.2 %>% 
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


obs.1 <- cbind(obs.1, rbind(obs.1.1, 
                            matrix(data = no.obs, 
                                   nrow = nrow(obs.1) - nrow(obs.1.1), 
                                   ncol = ncol(obs.1.1))))

obs.2 <- cbind(obs.2, rbind(as.matrix(obs.2.1), 
                            matrix(data = no.obs, 
                                   nrow = nrow(obs.2) - nrow(obs.2.1), 
                                   ncol = ncol(obs.2.1))))

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
             rbind(as.matrix(BUGS.data$u[,2,2:x]), 
                   matrix(data = 0, 
                          nrow = nrow(u.2) - dim(BUGS.data$u)[1], 
                          ncol = dim(BUGS.data$u)[3]-1)))

periods.1 <- c(all.periods$n.x, BUGS.data$periods[2:length(BUGS.data$periods)])
periods.2 <- c(all.periods$n.y, 0, BUGS.data$periods[3:4], rep(0, x-4))

Watch.Length.1[day.1 == 1 | day.1 > 89] <- 1
Watch.Length.2[day.2 == 1 | day.2 > 89] <- 1

# Day 1 and 90 should be zeros
n.1[day.1 == 1 | day.1 > 89] <- 0
n.2[day.2 == 1 | day.2 > 89] <- 0

N.1.inits <- n.1 * 3 + 2
for (k in 1:length(periods.1))
  N.1.inits[(periods.1[k]+1):nrow(N.1.inits), k] <- NA

N.1.obs <- matrix(nrow = nrow(N.1.inits),  ncol = ncol(N.1.inits))

N.1.inits[day.1 == 1 | day.1 > 89] <- NA
N.1.obs[day.1 == 1 | day.1 > 89] <- 0    # "partially observed" as in assumed zeros

N.2.inits <- n.2 * 3 + 2
for (k in 1:length(periods.2))
  N.2.inits[(periods.2[k]+1):nrow(N.2.inits), k] <- NA

N.2.obs <- matrix(nrow = nrow(N.2.inits),  ncol = ncol(N.2.inits))

N.2.inits[day.2 == 1 | day.2 > 89] <- NA
N.2.obs[day.2 == 1 | day.2 > 89] <- 0

# function to create initial values for N.1 and N.2. Without providing initial 
# values, the binomial likelihood doesn't work so well. Results in too small of
# N values.
N.inits <- function(N.1.inits, N.2.inits,n.chains){
  out.list <- vector(mode = "list", length = n.chains)
  return(lapply(out.list, FUN = function(x) list(N.1 = N.1.inits, N.2 = N.2.inits)))
}

jags.data <- list(n.1 = n.1,
                  n.2 = n.2,
                  N.1 = N.1.obs,
                  N.2 = N.2.obs,
                  n.station = ifelse(colSums(n.2, na.rm = T) > 0, 1, 0),
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

jags.params <- c("lambda.1",
                 "lambda.2",
                 "N.1", "N.2",
                 "OBS.RF.sp",
                 "OBS.Switch.sp.1",
                 "OBS.Switch.sp.2",
                 "BF.Switch.sp.1",
                 "BF.Switch.sp.2",
                 "BF.Fixed.sp.1",
                 "BF.Fixed.sp.2",
                 "VS.Switch.sp.1",
                 "VS.Switch.sp.2",
                 "VS.Fixed.sp.1",
                 "VS.Fixed.sp.2",
                 "mean.prob.sp.1",
                 "mean.prob.sp.2",
                 "BF.Fixed.sp.1",
                 "BF.Fixed.sp.2",
                 "VS.Fixed.sp.1",
                 "VS.Fixed.sp.2",
                 "Corrected.Est",
                 "Raw.Est",
                 "sp",
                 "Daily.Est",
                 "beta.sp",
                 "b.sp",
                 "sd.b.sp",
                 "log.lkhd.1",
                 "log.lkhd.2")

MCMC.params <- list(n.samples = 125000,
                    n.thin = 10,
                    n.burnin = 25000,
                    n.chains = 5)


out.file.name <- "RData/JAGS_Spline_results_All_Data_2024-04-22.rds"

if (!file.exists(out.file.name)){
  Start_Time<-Sys.time()
  
  jm <- jagsUI::jags(jags.data,
                     inits = N.inits(N.1.inits, N.2.inits, MCMC.params$n.chains),
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

# Some errors 2024-03-21:
# 
# Caught error when creating stat array for 'BF.Fixed.sp':
#   Error in array(NA, dim = apply(indices, 2, max)): negative length vectors are not allowed
# 
# 
# Caught error when creating stat array for 'VS.Fixed.sp':
#   Error in array(NA, dim = apply(indices, 2, max)): negative length vectors are not allowed
#   
# These have been fixed but forgot to record how.

# Summarize the output
Laake.seasons <- lapply(all.years, FUN = function(x) paste0(x, "/", x+1)) %>%
  unlist()

seasons <- c("2007/2008", "2009/2010", "2010/2011", 
             "2014/2015", "2015/2016", "2019/2020", 
             "2021/2022", "2022/2023", "2023/2024")

all.seasons <- c(Laake.seasons, seasons)

# Extract estimated counts
Daily.Est <- exp(jm.out$jm$sims.list$sp)
Corrected.Est <- jm.out$jm$sims.list$Corrected.Est

stats.list <- vector(mode = "list",
                     length = dim(Daily.Est)[3])

for (k in 1:dim(Daily.Est)[3]){
  # Daily.Est.list[[k]] <- Daily.Est[,,k]
  # Daily.Est.UCIs[[k]] <- apply(Daily.Est[,,k],2,quantile,0.975)
  # Daily.Est.LCIs[[k]] <- apply(Daily.Est[,,k],2,quantile,0.275)
  # 
  # sp.list[[k]] <- sp[,,k]
  # com.list[[k]] <- com[,,k]
  
  stats.list[[k]] <- data.frame(Daily.Est.median = apply(Daily.Est[,,k], 2,
                                                         median, na.rm = T),
                                Daily.Est.LCL = apply(Daily.Est[,,k], 2,
                                                      quantile,0.275, na.rm = T),
                                Daily.Est.UCL = apply(Daily.Est[,,k], 2,
                                                      quantile,0.975, na.rm = T),
                                #total.median = apply(exp(sp[,,k]), 1, sum),
                                days = 1:dim(Daily.Est)[2],
                                year = all.seasons[k])
}

all.stats <- do.call("rbind", stats.list) %>% group_by(year)

obs.day <- jm.out$jags.data$day.1
Nhats <- list()
k <- 1
for (k in 1:ncol(obs.day)){
  tmp.day <- jm.out$jags.data$day.1[,k] 
  tmp.n <- jm.out$jags.data$n.1[tmp.day > 1 & tmp.day < 90 & !is.na(tmp.day), k]
  tmp.watch <- jm.out$jags.data$Watch.Length.1[tmp.day > 1 & tmp.day < 90 & !is.na(tmp.day), k] %>% na.omit()
  
  Nhats[[k]] <- data.frame(day = tmp.day[tmp.day > 1 & tmp.day < 90 & !is.na(tmp.day)],
                           n = tmp.n,
                           watch = tmp.watch) %>%
    group_by(day) %>%
    summarize(day = first(day),
              n = sum(n),
              watch = sum(watch),
              Nhat = n/watch,
              year = all.seasons[k])

}

Nhats.df <- do.call("rbind", Nhats) %>% group_by(year)

ggplot(data = all.stats) + 
  geom_line(aes(x = days, y = Daily.Est.median),
            color = "darkorange") + 
  geom_ribbon(aes(x = days, 
                  ymin = Daily.Est.LCL, 
                  ymax = Daily.Est.UCL),
              fill = "orange", 
              alpha = 0.5) +
  geom_point(data = Nhats.df,
             aes(x = day, y = Nhat),
             shape = 4, color = "darkgreen", alpha = 0.5) +
  facet_wrap(vars(year)) +
  xlab("Days since December 1") + 
  ylab("Whales per day")

abundance.df <- data.frame(Year = lapply(str_split(all.seasons, "/"), 
                                         FUN = function(x) x[2]) %>% 
                             unlist() %>% 
                             as.numeric(),
                           Nhat = apply(Corrected.Est,
                                        FUN = mean,
                                        MARGIN = 2),
                           # Nhat = apply(Corrected.Est,
                           #              FUN = median,
                           #              MARGIN = 2),
                           SE = apply(Corrected.Est,
                                      FUN = function(x) sqrt(var(x)),
                                      MARGIN = 2),
                           CV = apply(Corrected.Est,
                                      FUN = function(x) 100*sqrt(var(x))/mean(x),
                                      MARGIN = 2),
                           Season = all.seasons,
                           CL.low = apply(Corrected.Est, 
                                          MARGIN = 2, 
                                          FUN = quantile, 0.025),
                           CL.high = apply(Corrected.Est, 
                                           MARGIN = 2, 
                                           FUN = quantile, 0.975),
                           Method = "Spline")

# Bring in the estimates from the other analyses (Laake's and Durban's) and compare:

Laake.estimates <- read_csv("Data/all_estimates_Laake_2024.csv",
                            col_types = cols(Year = col_integer(),
                                             Nhat = col_double(),
                                             SE = col_double(),
                                             CV = col_double(),
                                             Season = col_character(),
                                             CL.low = col_double(),
                                             CL.high = col_double())) %>%
  rename(Start.Year = Year) %>%
  mutate(Year = Start.Year + 1,
         Method = "Laake") %>%
  select(Year, Nhat, CV, SE, Season, CL.low, CL.high, Method)

# The most recent estimates (2023/2024) are here:
Durban.estimates <- read_csv(file = "Data/abundance_2024_85min.csv",
                             col_types = cols(Season = col_character(),
                                              total.mean = col_double(),
                                              total.CV = col_double(),
                                              total.median = col_double(),
                                              total.LCL = col_double(),
                                              total.UCL = col_double())) %>%
  transmute(Year = lapply(str_split(Season, "/"), 
                          FUN = function(x) x[2]) %>% 
              unlist() %>% 
              as.numeric(),
            Nhat = total.mean,
            CV = total.CV/100,
            SE = CV * Nhat,   
            Season = Season,
            CL.low = total.LCL,
            CL.high = total.UCL,
            Method = "Durban") %>%
  relocate(SE, .before = CV)


all.estimates <- rbind(abundance.df, Laake.estimates, Durban.estimates) 

ggplot(all.estimates) +
  geom_point(aes(x = Season, y = Nhat, color = Method)) +
  geom_errorbar(aes(x = Season, ymin = CL.low, ymax = CL.high, color = Method)) 

# Strange that estimates were lower than Laake's estimates but a lot higher for
# the last 9 years. I have the feeling data are not treated equally...
# 
# Take a look at watch length

jm.out$jags.data$Watch.Length.1 %>%
  as.data.frame() -> tmp

colnames(tmp) <- all.seasons

tmp %>% pivot_longer(cols = everything(), 
                     names_to = "year",
                     values_to = "Watch.Length")  %>%
  arrange(rev(desc(year))) -> watch.length.1.long

watch.length.1.long %>%
  group_by(year) %>%
  summarize(min = min(Watch.Length, na.rm = T),
            max = max(Watch.Length, na.rm = T),
            mean = mean(Watch.Length, na.rm = T)) -> watch.length.stats

ggplot(watch.length.stats)+
  geom_point(aes(x = year, y = mean)) 
  #geom_errorbar(aes(x = year, ymin = min, ymax = max))
