# LaakeData2Jags
# 
# Creates JAGS input data from Laake's ERAnalysis library.
# 
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

# the number of years in the dataset. A lot! 
all.years <- unique(Effort$Start.year)

# For jags and WinBugs code, what I need are
# 1. observed number of whales per day n[d, s, y], where d = # days since 12/1,
# s = station (1 = primary, 2 = secondary), y = year. For Laake's data, maximum 
# number of days per year was 94. 
# 2. Beaufort sea state bf[d,y]
# 3. Visibility code vs[d,y]
# 4. observer code obs[d,s,y]
# 5. the proportion of watch duration per day (out of 540 minutes or 9 hrs) watch.prop[d,y]
# 6. index of survey day, i.e., the number of days since 12/1 day[d,y]

# Need to count the number of days since 12-01 for each year.
# Then... 
# Count the number of whales per day and daily effort
# In early years, surveys were conducted 10 hrs. So, the watch proportion
# can be > 1.0, because we have used 9 hrs as maximum. 
Effort %>% 
  mutate(Day1 = as.Date(paste0(Start.year, "-12-01")),
         dt = as.numeric(as.Date(Date) - Day1)) %>%
  select(Start.year, nwhales, effort, vis, beaufort, Observer, dt) %>%
  group_by(Start.year, dt) %>%
  summarise(Start.year = first(Start.year),
            dt = first(dt),
            vs = max(vis),
            bf = max(beaufort),
            obs = first(Observer),
            effort = sum(effort),
            n = sum(nwhales)) %>%
  mutate(effort.min = effort * 24 * 60,
         watch.prop = effort.min/540) -> Effort.by.day

# Need to give numeric IDs to observers
Observer %>%
  mutate(ID.char = as.character(ID)) -> Observer

Observer.1 <- Observer[1:67,]

Effort.by.day %>%
  mutate(Initials = obs) %>%
  left_join(Observer, by = "Initials") %>%
  dplyr::select(-c(Initials, Observer, Name, Sex)) %>%
  rename(ID.1 = ID) %>%
  mutate(ID.char = obs) %>%
  left_join(Observer.1, by = "ID.char") %>%
  dplyr::select(-c(Initials, Observer, Name, Sex, ID.char)) -> Effort.by.day.1

Effort.by.day.1$ID[is.na(Effort.by.day.1$ID)] <- Effort.by.day.1$ID.1[is.na(Effort.by.day.1$ID)]
  
Effort.by.day %>% 
  select(Start.year) %>% 
  summarise(n = n()) -> n.year

# create day matrix - don't know how to do this in one line...
bf <- vs <- watch.prop <- day <- matrix(nrow = max(n.year$n), ncol = length(all.years))
n <- obs <- array(dim = c(max(n.year$n), 2, length(all.years)))
periods <- vector(mode = "numeric", length = length(all.years))
k <- 1
for (k in 1:length(all.years)){
  Effort.by.day.1 %>% 
    filter(Start.year == all.years[k]) -> tmp
    
  n[1:nrow(tmp), 1, k] <- tmp$n
  day[1:nrow(tmp), k] <- tmp$dt
  bf[1:nrow(tmp), k] <- tmp$bf
  vs[1:nrow(tmp), k] <- tmp$vs
  watch.prop[1:nrow(tmp), k] <- tmp$watch.prop
  
  obs[1:nrow(tmp), 1, k] <- tmp$ID
  
  periods[k] <- nrow(tmp)
}

jags.data <- list(n = n, 
                  n.station = rep(1, length(all.years)),
                  n.year = length(all.years),
                  n.obs = length(unique(Observer.1$ID)),
                  periods = periods,
                  obs = obs,
                  vs = scale(vs),
                  bf = scale(bf),
                  watch.prop = watch.prop,
                  day = day,
                  n.days = max(Effort.by.day.1$dt))

MCMC.params <- list(n.samples = 250000,
                    n.thin = 100,
                    n.burnin = 200000,
                    n.chains = 5)

out.file.name <- "RData/JAGS_pois_binom_results_Laake_Data.rds"
jags.model <- paste0("models/model_Richards_pois_bino.txt")


Start_Time<-Sys.time()

jm <- jagsUI::jags(jags.data.real,
                   inits = NULL,
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



# Using the filtered data, compute simple minded abundance estimates by treating the 
# sampled periods (watches) as a random sample of a period of a specified number
# of days (e.g., 4 days).  It uses raw counts with no correction for pod size or missed pods.
period=4
# compute the fraction of each period that was sampled (eg. 12 hours / (4 days *24 hrs/day))
sampled.fraction=with(Effort,
                      {
                        Day=as.numeric(as.Date(Date)-as.Date(paste(Start.year,"-12-01",sep="")))
                        tapply(effort,list(Start.year,cut(Day,seq(0,100,period))),sum)/period
                      })
# compute the number of whales counted in each period
whales.counted=with(Sightings,
                    {
                      Day=as.numeric(as.Date(Date)-as.Date(paste(Start.year,"-12-01",sep="")))
                      tapply(podsize,list(Start.year,cut(Day,seq(0,100,period))),sum)
                    })



