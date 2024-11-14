# WinBUGS_vs_Jags

# Compares results from WinBUGS Ver2.Rmd and Jags_Richards_AllData_V4.R
# The former analyzed data from 2015 to 2024, which were used for the 2024
# report of the abundance estimate. The latter is a new approach to the data,
# which I think is a better approach without using WinBUGS or its "cut" function.
#

rm(list = ls())
library(R2jags)
library(R2WinBUGS)
library(tidybayes)
library(tidyverse)

# Bring in the output from WinBUGS and Jags:
# 10-yr dataset with v2 extraction and minimum of 85 minutes observation
WinBUGS.out.10.years <- readRDS("RData/WinBUGS_10yr_v2_min85.rds")
Jags.out.10.years <- readRDS("RData/JAGS_Richards_pois_bino_v3_10yr_v2_2024-07-09.rds")
Jags.out.all.years <- readRDS("RData/JAGS_Richards_pois_bino_v4_AllYears_2024-10-21.rds")

# first compare input data to the two approaches
WinBUGS.data.10yr.n <- WinBUGS.out.10.years$BUGS.data$n 
Jags.data.10yr.n <- Jags.out.10.years$jags.data$n

# As it should be... no difference in observed counts
dif.n.10yr <- WinBUGS.data.10yr.n - Jags.data.10yr.n  

# It seems that Jags input for watch.prop, which is a measure of effort is a lot
# bigger values than those in WinBUGS input. This may explain the difference in
# estimated abundance... 
dif.watch.length <- WinBUGS.out.10.years$BUGS.data$Watch.Length - Jags.out.10.years$jags.data$watch.prop

# Bring in the Laake dataset as well as WinBUGS dataset to see what made this
# difference. The following code chunks are from LaakeData2WinBUGS.R and
# LaakeData2Jags.R

# WinBUGS:
data(Primary)      # on-effort sightings
data(ERSurveyData)
data("Observer")
data(PrimarySightings)
data(PrimaryEffort)
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
# 

# See LaakeAnalysis_NewData.R for the following code chunk
# Add new data
years <- c(2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024)
YEAR <- max(years)

# Observers need to be in integers, not their initials.
data("Observer")   # from ERAnalysis package.
#Observer$Set = "old"
# Need to give numeric IDs to observers
Observer %>%
  mutate(ID.char = as.character(ID)) -> Observer

# Lines 68 and 69 are duplicates. 
Observer.1 <- Observer[1:67,]

new.observers <- read.csv(file = "Data/ObserverList2023.csv") %>%
  transmute(ID = ID,
            Initials = obs)

Observer.1 %>%
  left_join(new.observers, by = "Initials") %>%
  filter(!is.na(ID.y)) %>%
  transmute(ID = ID.y,
            Initials = Initials) -> tmp

new.observers %>%
  anti_join(tmp, by = "ID") -> tmp.2

tmp.2$new.ID = 68:((68+nrow(tmp.2))-1)

tmp.2 %>%
  select(new.ID, Initials) %>%
  mutate(ID = new.ID,
         Initials = Initials) %>%
  select(-new.ID) -> tmp.3

all.observers <- rbind(Observer %>% select(ID, Initials), tmp.3) %>%
  mutate(Observer = Initials) %>%
  na.omit() %>%
  filter(ID != "21") %>%  # ID - 21 is ARV/AVS, which is not useful
  droplevels()

# add those initials back in
all.observers <- rbind(all.observers, data.frame(ID = c("21", "21"), 
                                                 Initials = c("ARV", "AVS"), 
                                                 Observer = c("ARV", "AVS")))

# sightings and effort
sightings.list.primary <- effort.list.primary  <- list()
sightings.list.secondary <- effort.list.secondary  <- list()
k <- 1
for (k in 1:length(years)){
  # These raw data files contain repeated observations of all groups.
  
  tmp.effort <- read.csv(paste0("Data/all_effort_", 
                                years[k], "_Tomo_v2.csv")) 
  
  tmp.effort %>%
    transmute(watch.key = watch.key,
              Start.year = years[k] - 1,
              key = key,
              begin = begin,
              end = end,
              npods = npods,
              nwhales = nwhales,
              effort = effort,
              vis = vis,
              beaufort = beaufort,
              Observer = toupper(observer),
              time = time,
              watch = shift,
              date.shift = date.shift,
              Use = T,
              Date = as.Date(Date, format = "%m/%d/%Y"),
              station = station) %>%
    filter(vis < 5, beaufort < 5) %>%
    na.omit() -> effort.all
  
  effort.all %>%
    filter(station == "P") -> effort.list.primary[[k]]
  
  effort.all %>%
    filter(station == "S") -> effort.list.secondary[[k]]
}

# Effort
effort.primary <- do.call("rbind", effort.list.primary)  %>%
  na.omit()

effort.primary %>% 
  left_join(all.observers, by = "Observer") %>%
  select(-c(Observer, Initials)) %>%
  rename(Observer = ID) -> effort.primary

effort.secondary <- do.call("rbind", effort.list.secondary)  %>%
  na.omit()

effort.secondary %>% 
  left_join(all.observers, by = "Observer") %>%
  select(-c(Observer, Initials)) %>%
  rename(Observer = ID) -> effort.secondary

##### end of code chunk from LaakeAnalysis_NewData.R
# Need to convert begin date to begin time, then to shift ID so they can be matched
# between primary and secondary sightings
frac.day2time <- function(x){
  decimal.day <- x - floor(x)
  hrs <- decimal.day * 24
  return(hrs)
}

# Convert fractional hours to shift ID
shift.definition.frac.hr <- function(time.dec.hrs){
  
  shift.id <- ifelse(time.dec.hrs <= 7.5, 0,
                     ifelse(time.dec.hrs <= 9, 1,
                            ifelse(time.dec.hrs <= 10.5, 2,
                                   ifelse(time.dec.hrs <= 12, 3,
                                          ifelse(time.dec.hrs <= 13.5, 4,
                                                 ifelse(time.dec.hrs <= 15, 5,
                                                        ifelse(time.dec.hrs <= 16.5, 6, 7)))))))
  return(shift.id)
}


Effort.1 %>% 
  mutate(Day1 = as.Date(paste0(Start.year, "-12-01")),
         dt = as.numeric(as.Date(Date) - Day1) + 1,
         obs = Observer,
         start.hr = frac.day2time(begin),
         shift = shift.definition.frac.hr(start.hr),
         station = 1) %>%
  select(Start.year, nwhales, effort, vis, beaufort, 
         obs, dt, begin, start.hr, shift, station) %>%
  filter(shift > 0) %>%
  group_by(Start.year) %>%
  mutate(effort.min = effort * 24 * 60,
         watch.prop = effort.min/540) -> Effort.1.by.period

Effort.2 %>% 
  mutate(Day1 = as.Date(paste0(Start.year, "-12-01")),
         dt = as.numeric(as.Date(Date) - Day1) + 1,
         obs = Observer,
         start.hr = frac.day2time(begin),
         shift = shift.definition.frac.hr(start.hr),
         station = 2) %>%
  select(Start.year, nwhales, effort, vis, beaufort, 
         obs, dt, begin, start.hr, shift, station) %>%
  filter(shift > 0) %>%
  group_by(Start.year) %>%
  mutate(effort.min = effort * 24 * 60,
         watch.prop = effort.min/540) -> Effort.2.by.period


Effort.1.by.period %>%
  mutate(Initials = obs) %>%
  left_join(Observer, by = "Initials") %>%
  dplyr::select(-c(Initials, Observer, Name, Sex)) %>%
  rename(ID.1 = ID) %>%
  mutate(ID.char = obs) %>%
  left_join(Observer.1, by = "ID.char") %>%
  dplyr::select(-c(Initials, Observer, Name, Sex, ID.char)) -> Effort.1.by.period.1

#Effort.by.day.1$ID[is.na(Effort.by.day.1$ID)] <- Effort.by.day.1$ID.1[is.na(Effort.by.day.1$ID)]
Effort.1.by.period.1$ID[is.na(Effort.1.by.period.1$ID)] <- Effort.1.by.period.1$ID.1[is.na(Effort.1.by.period.1$ID)]


create.WinBUGS.data <- function(in.data){
  # the number of years in the dataset. A lot! 
  all.years <- unique(in.data$Start.year)
  
  in.data %>% 
    group_by(Start.year) %>% 
    summarise(n = n()) -> n.year
  
  in.data.1 <- in.data
  # re-index observers
  obs.df <- data.frame(ID = unique(in.data.1$ID %>%  sort),
                       seq.ID = seq(1, length(unique(in.data.1$ID))))
  
  in.data.1 %>% 
    left_join(obs.df, by = "ID") -> in.data.1
  
  # create matrices - don't know how to do this in one line...  
  bf <- vs <- watch.prop <- day <- matrix(nrow = max(n.year$n), ncol = length(all.years))
  BUGS.day <- effort <- matrix(nrow = (max(n.year$n) + 2), 
                               ncol = length(all.years))
  
  BUGS.n <- matrix(data = 0, nrow= max(n.year$n), ncol= length(all.years))
  BUGS.obs <- matrix(data = nrow(obs.df)+1, nrow = max(n.year$n), ncol= length(all.years))
  
  periods <- vector(mode = "numeric", length = length(all.years))
  k <- 1
  for (k in 1:length(all.years)){
    in.data.1 %>% 
      filter(Start.year == all.years[k]) -> tmp
    
    BUGS.n[1:nrow(tmp), k] <- tmp$nwhales + 1
    BUGS.day[1:nrow(tmp), k] <- tmp$dt
    BUGS.day[(nrow(tmp)+1):(nrow(tmp)+2), k] <- c(1,90)
    bf[1:nrow(tmp), k] <- tmp$beaufort
    vs[1:nrow(tmp), k] <- tmp$vis
    effort[1:nrow(tmp), k] <- tmp$effort
    effort[(nrow(tmp)+1):(nrow(tmp)+2), k] <- c(1,1)
    BUGS.obs[1:nrow(tmp),  k] <- tmp$seq.ID
    
    periods[k] <- nrow(tmp)
  }
  
  BUGS.data <- list(n = BUGS.n,
                    n.com = BUGS.n,
                    n.sp = BUGS.n,
                    n.station = 2,
                    n.year = as.integer(length(all.years)),
                    n.obs = as.integer(length(unique(in.data.1$seq.ID))+1),
                    periods = as.integer(periods),
                    obs = BUGS.obs,
                    vs = vs,
                    bf = bf,
                    Watch.Length = effort,    
                    day = BUGS.day)
  
  return(BUGS.data)
}

BUGS.data.1 <- create.WinBUGS.data(Effort.1.by.period.1)

