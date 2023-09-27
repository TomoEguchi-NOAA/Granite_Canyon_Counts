# LaakeData2SplineJags
# 
# Creates JAGS input data from Laake's ERAnalysis library and runs JAGS on 
# Model_Nmix_Spline.txt. This model is identical to Durban's Nmixture model, except
# the normal (common) model was removed (consequently the cut function was removed).
# 
# The spline needs anchor points (days 1 and 90 in Durban's model) to work. In some years,
# data were collected after day 90. So, the end point needs to be redefined. Perhaps,
# it is 90, if survey ends before 90, but it is x days after the last survey day if
# the last day is after 90. x = 5, 10, ...?

# Code chunks with var1 = expressions are by Laake. I use var1 <- expressions

rm(list = ls())

library(ERAnalysis)
library(tidyverse)
library(ggplot2)

# From example code in the ERAnalysis library
# 
# The recent survey data 1987 and after are stored in ERSurveyData and those data
# are processed by the ERAbund program to produce files of sightings and effort.
# The sightings files are split into Primary observer and Secondary observer sightings.
# Primary observer sightings are whales that are not travelling North and are defined by
# those when EXPERIMENT==1 (single observer) or a designated LOCATION when EXPERIMENT==2.
#  For surveys 2000/2001 and 2001/2002, the primary observer was at LOCATION=="N"
# and for all other years, LOCATION=="S".
#
# Based on the projected timing of the passage of the whale (t241) perpendicular to the 
# watch station, the sighting was either contained in the watch (on effort) or not (off effort).
# The dataframe Primary contains all of the on effort sightings and PrimaryOff contains all
# of the off-effort sightings.  
#
#data(PrimaryOff)   # off-effort sightings
data(Primary)      # on-effort sightings
data(ERSurveyData)
data("Observer")

# The data in PrimarySightings are all southbound sightings for all years in which visibility and beaufort
# are less than or equal to 4. Below the counts are shown for the 2 dataframes for
# recent surveys since 1987/88.
#table(Primary$Start.year[Primary$vis<=4 & Primary$beaufort<=4])
data(PrimarySightings)
data(PrimaryEffort)

# Likewise, the secondary sightings are those with EXPERIMENT==2 but the LOCATION that
# is not designated as primary.  there is no effort data for the secondary sightings... 
# so, can't use it for BUGS/jags - ignore it for now.
# data(SecondarySightings)

# Effort and sightings prior to 1987 were filtered for an entire watch if vis or beaufort 
# exceeded 4 at any time during the watch.  This is done for surveys starting in 1987 with the
# Use variable which is set to FALSE for all effort records in a watch if at any time the vis or
# beaufort exceeded 4 during the watch.
# Here are the hours of effort that are excluded (FALSE) and included (TRUE) by each year
# Note that for most years <1987 there are no records with Use==FALSE because the filtered records
# were excluded at the time the dataframe was constructed. The only exception is for 1978 in which  
# one watch (5 hours) was missing a beaufort value so it was excluded.
# tapply(PrimaryEffort$effort,
#        list(PrimaryEffort$Use,
#             PrimaryEffort$Start.year),
#        sum)*24
#        

# Filter effort and sightings and store in dataframes Effort and Sightings
Effort = PrimaryEffort[PrimaryEffort$Use,]  

Sightings = PrimarySightings
Sightings$seq = 1:nrow(Sightings)
Sightings = merge(Sightings, subset(Effort, select=c("key")))
Sightings = Sightings[order(Sightings$seq),]

# filter off-effort sightings and high Beaufort/vis lines from secondary sightings
# but... there is no effort data for the secondary sightings... so, can't use it
# for BUGS/jags - ignore it for now.
# SecondarySightings %>% 
#   filter(vis < 5, beaufort < 5, is.na(off)) -> secondary.sightings

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

Effort %>% 
  mutate(Day1 = as.Date(paste0(Start.year, "-12-01")),
         dt = as.numeric(as.Date(Date) - Day1) + 1,
         obs = Observer) %>%
  select(Start.year, nwhales, effort, vis, beaufort, obs, dt) %>%
  group_by(Start.year) %>%
  mutate(effort.min = effort * 24 * 60,
         watch.prop = effort.min/540) -> Effort.by.period

  # group_by(Start.year, dt) %>%
  # summarise(Start.year = first(Start.year),
  #           dt = first(dt),
  #           vs = max(vis),
  #           bf = max(beaufort),
  #           obs = first(Observer),
  #           effort = sum(effort),
  #           n = sum(nwhales)) %>%
  # mutate(effort.min = effort * 24 * 60,
  #        watch.prop = effort.min/540) -> Effort.by.day

# Need to give numeric IDs to observers
Observer %>%
  mutate(ID.char = as.character(ID)) -> Observer

# Lines 68 and 69 are duplicates. 
Observer.1 <- Observer[1:67,]

#Effort.by.day %>%
Effort.by.period %>%
  mutate(Initials = obs) %>%
  left_join(Observer, by = "Initials") %>%
  dplyr::select(-c(Initials, Observer, Name, Sex)) %>%
  rename(ID.1 = ID) %>%
  mutate(ID.char = obs) %>%
  left_join(Observer.1, by = "ID.char") %>%
  dplyr::select(-c(Initials, Observer, Name, Sex, ID.char)) -> Effort.by.period.1

#Effort.by.day.1$ID[is.na(Effort.by.day.1$ID)] <- Effort.by.day.1$ID.1[is.na(Effort.by.day.1$ID)]
Effort.by.period.1$ID[is.na(Effort.by.period.1$ID)] <- Effort.by.period.1$ID.1[is.na(Effort.by.period.1$ID)]

create.jags.data <- function(Effort.by.period.1){
  # the number of years in the dataset. A lot! 
  all.years <- unique(Effort.by.period.1$Start.year)
  
  Effort.by.period.1 %>% 
    select(Start.year) %>% 
    summarise(n = n()) -> n.year

  # re-index observers
  obs.df <- data.frame(ID = unique(Effort.by.period.1$ID %>% sort),
                       seq.ID = seq(1, length(unique(Effort.by.period.1$ID))))
  
  Effort.by.period.1 %>% 
    left_join(obs.df, by = "ID") -> Effort.by.period.1
  
  # create matrices - don't know how to do this in one line...  
  bf <- vs <- matrix(nrow = max(n.year$n), ncol = length(all.years))

  day <- effort <- matrix(nrow = (max(n.year$n) + 2), 
                          ncol = length(all.years))
  
  n <- obs <- array(dim = c(max(n.year$n), 2, length(all.years)))
  
   periods <- vector(mode = "numeric", length = length(all.years))
  k <- 1
  for (k in 1:length(all.years)){
    Effort.by.period.1 %>% 
      filter(Start.year == all.years[k]) -> tmp
    
    n[1:nrow(tmp), 1, k] <- tmp$nwhales
    day[1:nrow(tmp), k] <- tmp$dt
    
    day[(nrow(tmp)+1):(nrow(tmp)+2), k] <- c(1,100)
    bf[1:nrow(tmp), k] <- tmp$beaufort
    vs[1:nrow(tmp), k] <- tmp$vis
    #watch.prop[1:nrow(tmp), k] <- tmp$watch.prop
    effort[1:nrow(tmp), k] <- tmp$effort
    effort[(nrow(tmp)+1):(nrow(tmp)+2), k] <- c(1,1)
    obs[1:nrow(tmp), 1, k] <- tmp$seq.ID

    periods[k] <- nrow(tmp)
  }
  
  Effort.by.period.1 %>% 
    group_by(Start.year) %>%
    summarise(max.dt = max(dt)) %>%
    pull(max.dt) -> max.dt
  
  jags.data <- list(n.sp = abind::abind(n, 
                                        array(NA, dim = c(2, 2, dim(n)[3])), 
                                        along = 1), 
                    n.station = rep(1, length(all.years)),
                    n.year = length(all.years),
                    n.obs = length(unique(Effort.by.period.1$seq.ID)),
                    periods = periods,
                    obs = obs,
                    vs = scale(vs),
                    bf = scale(bf),
                    vs.raw = vs,
                    bf.raw = bf,
                    Watch.Length = effort,
                    day = day,
                    n.days = rep(100, times = dim(n)[3])) #max(Effort.by.period.1$dt))
  
  jags.data$knot <-  c(-1.46,-1.26,-1.02,-0.78,
                       -0.58,-0.34,-0.10,0.10,
                       0.34,0.57,0.78,1.02,1.26,1.46)
  
  jags.data$n.knots <- 14
  
  u <- array(data = 0, dim = dim(jags.data$n.sp))
  
  for (y in 1:jags.data$n.year){
    for (s in 1:jags.data$n.station[y]){
      for (p in 1:jags.data$periods[y]){
        
        u[p,s,y] <- 1
      }
    }
  }
  
  jags.data$u <- u
  
  return(jags.data)
}

jags.data <- create.jags.data(Effort.by.period.1)

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
                 "Daily.Est","mean.beta",
                 "beta.sigma","beta","beta.sp","b.sp","sd.b.sp")

MCMC.params <- list(n.samples = 250000,
                    n.thin = 100,
                    n.burnin = 200000,
                    n.chains = 5)

# Needs initial values for N to avoid errors. N has to be significantly larger
# than observed n.
jags.inits <- function() list(N.sp = jags.data$n.sp[,1,] * 2 + 10)

# The first attempt used counts per day, which worked fine and estimates were
# pretty close to Laake's estimates. But, the model uses all data at observation
# period level
#out.file.name <- "RData/JAGS_pois_binom_results_Laake_Data.rds"

out.file.name <- "RData/JAGS_Spline_results_Laake_Data.rds"
jags.model <- paste0("models/Model_Nmix_Spline.txt")


if (!file.exists(out.file.name)){
  Start_Time<-Sys.time()
  
  jm <- jagsUI::jags(jags.data,
                     inits = jags.inits,
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
                 System = Sys.getenv())
  
  saveRDS(jm.out,
          file = out.file.name)
  
} else {
  jm.out <- readRDS(out.file.name)
}



