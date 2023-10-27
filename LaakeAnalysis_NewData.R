#LaakeAnalysis_NewData
# Applies Laake's analysis on data since 2015

rm(list = ls())

library(ERAnalysis)
library(tidyverse)
library(ggplot2)
library(lubridate)

# For Laake's approach, I need primary effort, primary sightings, secondary effort,
# secondary sightings, distance, podsize, and associated visibility and Beaufort 
# sea state

# These are from Laake's code - good reference to compare to:
# sightings
Laake_PrimarySightings <- read.csv(file = "Data/Laake_PrimarySightings.csv")
Laake_SecondarySightings <- read.csv(file = "Data/Laake_SecondarySightings.csv")

# effort
Laake_PrimaryEffort <- read.csv(file = "Data/Laake_PrimaryEffort.csv") %>%
  filter(Use)
Laake_SecondaryEffort <- read.csv(file = "Data/Laake_SecondaryEffort.csv") %>%
  filter(Use)

# Necessary information (all sightings) was extracted from raw data files using
# Extract_Data_All_v2.Rmd (version of 2023-10-19)
# For these years, there were no secondary observers. I don't have raw data for
# previous years. 
# 
# In the Laake data format, the "key" is Date_groupID. The key variable is essential
# for the analysis. From the help file:
# The sightings and effort dataframes have some fairly specific requirements. 
# First and foremost, they must each have a field named key which is used to link 
# the sightings contained in a specific period of effort. 
# 
# Prior to 1987, effort was only recorded as an entire watch period but starting 
# in 1987 the watch period was broken into segments with constant environmental 
# data like visibility and wind force (beaufort). These smaller segments are called 
# effort periods. So for surveys before 1987 there is a single effort period for 
# each watch and for 1987 and after there can be several effort periods per watch. 
# 
# These small periods are identified by "key."
# 
# Effort must have fields named Start.year which is the year the survey began, 
# begin and end which are the decimal day values for the begin and end times for 
# the effort period, time which is the mid-point of the interval and effort which 
# is the length of the interval. It can also have optional fields vis and beaufort 
# which are the values during the effort period.
# 
# The sightings dataframe must have a field podsize, the recorded observed pod size, 
# and optionally corrected.podsize. All other fields are optional but there must be 
# a valid value (not NA) for each field used in dformula for computing detection probability.
# t241 is the computed time (decimal days since 1 Dec) at which whale crosses line 
# extending from shed at angle of 241; Used south or north time prior 1979; 
# 1979+ calculated value using 3.24 nm whale speed
# 
years <- c(2015, 2016, 2020, 2022, 2023)

#years <- 2015
# sightings and efort
sightings.list <- effort.list <- list()
k <- 1
for (k in 1:length(years)){
  # These raw data files contain repeated observations of all groups.
  # 
   
  tmp.sightings <- read.csv(paste0("data/all_sightings_", 
                         years[k], "_Tomo_v2.csv")) 
  
  tmp.sightings %>%
    mutate(Date1 = as.Date(Date, format = "%m/%d/%Y")) %>%
    #group_by(group) %>%
    transmute(Date = Date1,
              Time = Time_PST, 
              day = day(Date),
              month = month(Date),
              year = year(Date),
              watch = Shift,
              t241 = difftime((paste(Date, Time)),
                              (paste0((years[k] - 1), 
                                             "-12-01 00:00:00"))) %>%
                as.numeric(),
              Group_ID = Group_ID,
              distance = Distance,
              podsize = nwhales,
              vis = vis,
              beaufort = beaufort,
              Start.year = years[k]-1,
              Observer = observer,
              key = key) %>%
    arrange(Date, Group_ID) %>%
    filter(vis < 5, beaufort < 5) %>%
    na.omit() -> sightings.list[[k]]
  
  tmp.effort <- read.csv(paste0("data/all_effort_", 
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
              time = time,
              watch = shift,
              Use = T,
              Date = as.Date(Date, format = "%m/%d/%Y")) %>%
    filter(vis < 5, beaufort < 5) %>%
    na.omit() -> effort.list[[k]]
}

# sightings
sightings <- do.call("rbind", sightings.list) %>%
  na.omit()

# Need to replace observer code to integers.
all.observers <- unique(sightings$Observer)
observer.df <- data.frame(Observer = all.observers,
                          code = c(1:length(all.observers)))

sightings %>% left_join(observer.df, by = "Observer") %>%
  select(-Observer) %>%
  rename(Observer = code) -> sightings

# Effort
effort <- do.call("rbind", effort.list)  %>%
  na.omit()

# gsS: nmax x nmax pod size calibration matrix; each row is a true pod size 
# from 1 to nmax and the value for each column is the probability that a pod of 
# a true size S is recorded as a size s (1..nmax columns)
# 
naive.abundance.models <- list()
for (k in 1:length(years)){
  naive.abundance.models[[k]] <- estimate.abundance(spar = NULL,
                                                    dpar = NULL,
                                                    gsS = gsS,  # Does not matter when spar = NULL and dpar = NULL
                                                    effort =effort.list[[k]], 
                                                    sightings =sightings.list[[k]], 
                                                    final.time = 90,
                                                    lower.time = 0,
                                                    gformula = ~s(time),
                                                    dformula = NULL)
  
}

# How do I adjust these estimates for other factors...? Perhaps... I can add these
# new data to Laake's primary dataframes (sightings and effort)

sightings %>%
  transmute(X = seq(from = 35608, to = 35608 + nrow(sightings) - 1),
            Date = Date, day = day, month = month,
            year = year, 
            watch = watch,
            t241 = t241, distance = distance,
            podsize = podsize,
            vis = vis,
            beaufort = beaufort,
            wind.direction = NA,
            key = key,
            pphr = NA,
            Start.year = Start.year,
            original.watch = watch,
            only = TRUE,
            hours = NA,
            Sex = NA,
            Observer = Observer) -> sightings.Laake.format

sightings.all <- rbind(Laake_PrimarySightings, sightings.Laake.format)

# Need to add watch.key, which is yyyy_seq, where seq is the sequential number 
# of watch. 
effort %>%
  transmute(X = seq(from = 7293, to = 7292 + nrow(effort)),
            watch.key = )