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
    arrange(Date, Group_ID) -> sightings.list[[k]]
  
  tmp.effort <- read.csv(paste0("data/all_effort_", 
                                years[k], "_Tomo_v2.csv")) 
  
  tmp.effort %>%
    transmute(Start.year = years[k] - 1,
              key = key,
              begin = begin,
              end = end,
              npods = npods,
              nwhales = nwhales,
              time = time,
              effort = effort,
              vis = vis,
              beaufort = beaufort,
              Date = as.Date(Date, format = "%m/%d/%Y")) -> effort.list[[k]]
}

sightings <- do.call("rbind", sightings.list) %>%
  na.omit()
  #select(Date, Shift, Distance, n, Visibility, Beaufort, Observer1) %>%
  # transmute(Date = as.Date(Date, format = "%m/%d/%Y"),
  #           day = day(Date),
  #           month = month(Date),
  #           year = year(Date),
  #           watch = Shift,
  #           t241 = NA,
  #           distance = Distance,
  #           podsize = n,
  #           vis = Visibility,
  #           beaufort = Beaufort,
  #           Start.year = years[k]-1,
  #           Observer = Observer1) 

# Effort
effort <- do.call("rbind", effort.list)  %>%
  na.omit()
  # select(-ff) %>%
  # transmute(Start.year = Start.year,
  #           begin = begin,
  #           end = end,
  #           nwhales = n,
  #           effort = dur,
  #           vis = vs,
  #           beaufort = bf,
  #           Observer = obs,
  #           watch = i)


