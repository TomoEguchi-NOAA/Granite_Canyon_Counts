#LaakeAnalysis_NewData
# Applies Laake's analysis on data since 2015

# This script is used to run the method of Laake et al. 2012. Their approach 
# requires matched sightings between two stations (primary and secondary). Such
# information is lacking in the data from the 2009/2010 season to the most recent
# data (2024/2025). So, I used the ratio between "naive" and "corrected" estiamtes
# for earlier results (up to 2006/2007) to adjust for new estimates. 
# 
# Because I cannot build the ERAnalysis package in newer R versions, I created
# a script that contains all functions that were used in Laake's analysis (Laake_functions.R).
# Raw data from recent years (i.e., 2009/2010 to current) need to be extracted
# using Extract_Data_All_v2.Rmd. Output of this script can be saved into a CVS
# file by making save.file <- TRUE in the beginning of this script. 

rm(list = ls())

#library(ERAnalysis)
library(tidyverse)
library(ggplot2)
library(lubridate)

source("Granite_Canyon_Counts_fcns.R")
source("Laake_functions.R")

save.file <- F
# These are the second year of each season.
years <- c(2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024, 2025)
YEAR <- max(years)
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

# Observers need to be in integers, not their initials.
#data("Observer")   # from ERAnalysis package.
# The Observer rda file was converted into a csv file.
Obserer <- read.csv("Data/Observer_Laake.csv")

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

# Only 2010 and 2011 had two observation stations in the recent years, which are 
# needed for applying Laake's method beyond naive estimates. Even for those years,
# there is no information about which pod was sighted by one or both stations.
# So, the method cannot be applied directly to all data. 

# sightings and effort
sightings.list.primary <- effort.list.primary  <- list()
sightings.list.secondary <- effort.list.secondary  <- list()

# Bring in processed data and create dataframes that match the format of Laake's
# input data
for (k in 1:length(years)){
  # These raw data files contain repeated observations of all groups.
  # They are created from Extract_Data_All_v2.Rmd 
  tmp.sightings <- read.csv(paste0("Data/all_sightings_", 
                                   years[k], "_Tomo_v2.csv")) 
  
  tmp.sightings %>%
    mutate(Date1 = as.Date(Date, format = "%m/%d/%Y")) %>%
    #group_by(group) %>%
    transmute(key = as.factor(key),
              Date = Date1,
              Time = Time_PST, 
              day = day(Date),
              month = month(Date),
              year = year(Date),
              watch = shift,
              t241 = difftime((paste(Date, Time)),
                              (paste0((years[k] - 1), 
                                             "-11-30 00:00:00"))) %>%
                as.numeric(),
              distance = Distance,
              podsize = nwhales,
              vis = vis,
              beaufort = beaufort,
              wind.direction = NA,
              pphr = pphr,
              Start.year = years[k]-1,
              original.watch = NA,
              only = TRUE,
              hours = floor(time - min(time)),
              Sex = NA,   # observer's gender
              Observer = toupper(observer),
              seq = seq(1, nrow(tmp.sightings)),
              Group_ID = Group_ID,
              station = station) %>%
    arrange(Date, Group_ID) %>%
    filter(vis < 5, beaufort < 5) %>%
    select(-c(Group_ID, Time)) -> sightings.all
    #na.omit() -> sightings.all

  sightings.all %>% 
    filter(station == "P") -> sightings.list.primary[[k]]
  
  sightings.all %>% 
    filter(station == "S") -> sightings.list.secondary[[k]]
  
  tmp.effort <- read.csv(paste0("Data/all_effort_", 
                                years[k], "_Tomo_v2.csv")) 
  
  tmp.effort %>%
    transmute(watch.key = as.factor(watch.key),
              Start.year = years[k] - 1,
              key = as.factor(key),
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
              #date.shift = date.shift,
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

# sightings for all years combined
sightings.primary <- do.call("rbind", sightings.list.primary) 

# if (!file.exists(paste0("Data/all_observers_", max(years), ".csv"))){
  # # Find Observers in Laake's observer list who are also in the new observer list
  # # There are multiple initials per person in some cases. ID numbers should be the
  # # same for these initials
  # Observer %>%
    # filter(is.na(Observer)) %>%
    # droplevels() %>%
    # mutate(Identifier = Initials) %>%
    # select(ID, Identifier, Name) -> Observer.NA
  
  # Observer %>%
    # filter(is.na(Initials)) %>%
    # droplevels() %>%
    # mutate(Identifier = as.factor(Observer)) %>%
    # select(ID, Identifier, Name) -> Observer.Initial.NA
  
  # Observer.1 <- rbind(Observer.NA, Observer.Initial.NA)
  
  # uniq.Observer.1 <- data.frame(Name = unique(Observer.1$Name),
                                # ID.new = seq(1:length(unique(Observer.1$Name))))
  
  # Observer.1 %>%
    # left_join(uniq.Observer.1, by = "Name") %>%
    # select(-c(ID, Name)) %>%
    # rename(ID = ID.new) -> Observer.2
  
  # # Find observers from more recent data (2009/2010 to present)
  # Observers.new <- unique(sightings.primary$Observer)
  # new.observers <- data.frame(Identifier = Observers.new,
                              # ID = seq((max(Observer.2$ID)+1), 
                                       # (max(Observer.2$ID) + length(Observers.new))))
  
  # # Combine "old" and "new" observer lists
  # Observer.2 %>%
    # left_join(new.observers, by = "Identifier") %>% #
    # filter(!is.na(ID.y)) %>%
    # transmute(ID = ID.y,
              # Identifier = Identifier) -> tmp
  
  # # Remove those (tmp) from the new observer list
  # new.observers %>%
    # anti_join(tmp, by = "ID") -> tmp.2
  
  # tmp.4 <- rbind(Observer.2, tmp.2) %>%
    # mutate(Observer = Identifier) %>%
    # select(-Identifier) %>%
    # na.omit() 
  
  # # Fix one observer because "ARV/AVS" is not usable 
  # ARV.ID <- tmp.4[tmp.4$Observer == "ARV/AVS", "ID"]
  # tmp.4 %>%
    # filter(Observer != "ARV/AVS") %>%  
    # droplevels() -> tmp.5
  
  # # add those initials back in with the same ID
  # all.observers <- rbind(tmp.5, 
                         # data.frame(ID = c(ARV.ID, ARV.ID), 
                                    # Observer = c("ARV", "AVS")))
  
  # # Write the new list to a file
  # write.csv(all.observers,
            # file = paste0("Data/all_observers_", max(years), ".csv"))
# } else {
  # all.observers <- read.csv(paste0("Data/all_observers_", max(years), ".csv"))
# }

# 2025-06-18 Changed the previous block to a function call. 
all.observers <- create.observer.list(sightings.primary)

# Rearrange sightings dataframe columns
sightings.primary %>% 
  left_join(all.observers, by = "Observer") %>%
  select(-ID) %>%
  #rename(Observer = ID) %>%
  mutate(Observer = as.factor(Observer)) %>%
  relocate(Observer, .after = Sex) %>%
  select(-station) -> sightings.primary

sightings.secondary <- do.call("rbind", sightings.list.secondary) 

sightings.secondary %>% 
  left_join(all.observers, by = "Observer") %>%
  select(-ID) %>%
  #rename(Observer = ID) %>%
  mutate(Observer = as.factor(Observer)) %>%
  relocate(Observer, .after = Sex) %>%
  select(-station) -> sightings.secondary

# Effort
effort.primary <- do.call("rbind", effort.list.primary)  

effort.primary %>% 
  left_join(all.observers, by = "Observer") %>%
  select(-ID) %>%
  #rename(Observer = ID) %>% 
  mutate(Observer = as.factor(Observer)) %>%
  select(-station) %>%
  relocate(Observer, .after = beaufort) -> effort.primary

effort.secondary <- do.call("rbind", effort.list.secondary)  

effort.secondary %>% 
  left_join(all.observers, by = "Observer") %>%
  select(-ID) %>%
  #rename(Observer = ID) %>%
  mutate(Observer = as.factor(Observer)) %>%
  select(-station) %>%
  relocate(Observer, .after = beaufort) -> effort.secondary

# gsS: nmax x nmax pod size calibration matrix; each row is a true pod size 
# from 1 to nmax and the value for each column is the probability that a pod of 
# a true size S is recorded as a size s (1..nmax columns)
naive.abundance.models.new <- list()
for (k in 1:length(years)){
  sightings.primary %>%
    filter(Start.year == (years[k]-1)) -> sightings
  
  effort.primary %>%
    filter(Start.year == (years[k]-1)) -> effort
  
  naive.abundance.models.new[[k]] <- estimate.abundance(spar = NULL,
                                                        dpar = NULL,
                                                        gsS = gsS,  # Does not matter when spar = NULL and dpar = NULL
                                                        effort =effort, 
                                                        sightings =sightings, 
                                                        final.time = 90,
                                                        lower.time = 0,
                                                        gformula = ~s(time),
                                                        dformula = NULL)
  
}

# Without sightings from the secondary team, it seems that I can't go beyond
# computing naive abundance estimates. 2023-12-01
# 
# In Laake's analysis, they computed the multiplication factor for naive-to-true conversion from 
# years with two stations. Then the average of those conversion factors was used for years without
# two stations. 
# 
# # Define set of models to be evaluated for detection
models=c("podsize+Dist+Observer",
         "podsize+Dist+Observer+beaufort",
         "podsize+Dist+Observer+vis")

# # compute.series (or compute.series.new) can take a Match dataframe, which
# # contains information about sighting and non-sighting of each whale pod by
# # primary and secondary observers. The match dataframe from Laake's analysis
# # was obtained from running their code and saving the output Match to a .csv
# # file. A similar dataframe needs to be created for 2010 and 2011.
# 
# Match.Laake <- read.csv("Data/match_1987_2006.csv")
# 
# # for new datasets, only 2010 and 2011 have the secondary observers:
# # need to create the "seen" variable, which is 0/1
# sightings.primary %>%
#   filter(Start.year == 2009 | Start.year == 2010)  %>%
#   mutate(seen = 1,
#          new.key = ifelse(Group_ID > 9, 
#                           paste0(Date, "-", Group_ID),
#                           paste0(Date, "-0", Group_ID))) -> sightings.primary.1
# 
# sightings.secondary %>%
#   filter(Start.year == 2009 | Start.year == 2010)  %>%
#   mutate(seen = 1,
#          new.key = ifelse(Group_ID > 9, 
#                           paste0(Date, "-", Group_ID),
#                           paste0(Date, "-0", Group_ID))) -> sightings.secondary.1
# 
# # Need to reduce the primary sightings to match secondary survey dates:
# unique.secondary.dates <- unique(sightings.secondary.1$Date)
# filtered.primary.sightings.1 <- sightings.primary.1[sightings.primary.1$Date %in% unique.secondary.dates, ]
# filtered.primary.effort <- effort.primary[effort.primary$Date %in% unique.secondary.dates, ]
# 
# # Combine with the secondary sightings
# rbind(filtered.primary.sightings.1, sightings.secondary.1) %>%
#   arrange(new.key) -> sightings.all.1
# 
# # Seen by both primary and secondary
# filtered.primary.sightings.1 %>%
#   semi_join(sightings.secondary.1 %>%
#               select(new.key), by = "new.key") -> seen.by.both
#                            
# sightings.all.1 %>%
#   filter(new.key %in% seen.by.both$new.key) -> sightings.seen.by.both
# 
# # Seen by primary, but not seen by secondary
# filtered.primary.sightings.1 %>%
#   anti_join(sightings.secondary.1 %>%
#               select(new.key), by = "new.key") -> seen.by.primary.only
# 
# seen.by.primary <- seen.by.primary.only
# seen.by.primary$seen <- 0
# seen.by.primary$station <- "S"
# sightings.by.primary <- rbind(seen.by.primary.only, seen.by.primary) 
# 
# # Seen by secondary, but not seen by primary
# sightings.secondary.1 %>%
#   anti_join(sightings.primary.1 %>%
#               select(new.key), by = "new.key") -> seen.by.secondary.only
# 
# seen.by.secondary <- seen.by.secondary.only
# seen.by.secondary$seen <- 0
# seen.by.secondary$station <- "P"
# sightings.by.secondary <- rbind(seen.by.secondary.only, seen.by.secondary) 
# 
# sightings.all.2 <- rbind(sightings.seen.by.both,
#                          sightings.by.primary,
#                          sightings.by.secondary) %>%
#   arrange(new.key) 
#   
# sightings.all.2$seq <- seq(1, nrow(sightings.all.2))
# 
# # Need to find right observers for not-seen pods.
# # date.shift is the index for date and observer combinations
# not.seen <- filter(sightings.all.2, seen == 0)
# k <- 1
# for (k in 1:nrow(not.seen)){
#   if (not.seen$station[k] == "S"){
#     effort.secondary %>%
#       filter(date.shift == not.seen$date.shift[k]) -> effort.k
#     if (nrow(effort.k) > 0){
#       not.seen$Observer[k] <- effort.k$Observer[1]
#     } else {
#       not.seen$Observer[k] <- NA
#     }
#   } else {
#     filtered.primary.effort %>%
#       filter(date.shift == not.seen$date.shift[k]) -> effort.k
#     if (nrow(effort.k) > 0){
#       not.seen$Observer[k] <- effort.k$Observer[1]
#     } else {
#       not.seen$Observer[k] <- NA
#     }
#   }
# }
# 
# sightings.all.3 <- rbind(sightings.all.2 %>% filter(seen == 1),
#                          not.seen) %>% 
#   arrange(new.key) 
# 
# sightings.all.3 %>%
#   transmute(X = seq(1, nrow(sightings.all.2)),
#             key = key,
#             Date = Date,
#             seen = seen, 
#             station = station,
#             day = day,
#             watch = watch,
#             t241 = t241,
#             distance = distance,
#             podsize = podsize,
#             vis = vis,
#             beaufort = beaufort,
#             wind.direction = NA,
#             Start.year = Start.year,
#             seq = seq,
#             hours = NA,
#             pphr = pphr,
#             Sex = NA,
#             Observer = Observer,
#             Use = TRUE,
#             new.key = new.key) -> Match.new
# 
# # Somehow we need to remove NAs in Observers... but can't remove just those lines
# # because Primary and Secondary observations need to match.
# # When there was no secondary or primary observers, the other has to be removed.
# Match.new %>% 
#   filter(is.na(Observer)) %>%
#   select(new.key) %>%
#   unique() -> new.key.obs.NA
# 
# Match.new %>% 
#   anti_join(new.key.obs.NA, by = "new.key") -> Match.new.1
# 
# Match.all <- rbind(Match.Laake, Match.new.1 %>% select(-new.key)) 
# 
# # In Match.new, it appears that there are 2544 Primary and 2534 Secondary...
# # There should be the same number of primary and secondary...  
# # Debugging here...
# # Match.new %>% group_by(new.key) %>% summarize(n.new.key = n()) -> new.key.summary
# # new.key.summary %>% filter(n.new.key > 2) -> problem.new.keys
# # Match.new %>%
# #   filter(new.key == problem.new.keys$new.key[1])
# 
# # 2010-01-25-38 in 2010-01-25_4 pphr is 14.446 AND 0.963 in two lines at 12:00:30
# # The data file has been fixed for 2010-01-25. The problem was when the observers
# # didn't change until a few minutes after 1200. Sightings occurred and observers
# # changed without specifically having a line to define it (I might have deleted
# # the line...). I created a short "shift" that consisted of the observers from 
# # the previous shift and a new shift starting right after with new observers.
# 
# 
# # Next compute the series of abundance estimates for 8 years plus 
# # 2 years for which secondary observers exist (2010, 2011) by
# # fitting and selecting the best detection model but not applying the pod size correction.
# # From those 10 estimates and the naive estimates, compute an average ratio and 
# # apply it to generate the estimates for the first 15 surveys prior to 1987 and the most
# # recent 5 years (2015, 2016, 2020, 2022, 2023).
# # years <- c(2010, 2011, 2015, 2016, 2020, 2022, 2023)
# # # I created a new function compute.series.new by adding 2010 and 2011.
# # DistBreaks=c(0,1,2,3,4,20)
# # cutoff=4
# # 
# 
# # Sightings <- merge(sightings.primary,
# #                    subset(effort.primary,
# #                           select=c("key","Use")), by="key")
# 
# # 2023-12-15 Need to find correct inputs here; sightings, effort, and match
# # Have to combine Laake's inputs and new ones here
# sightings.primary %>%
#   transmute(X = seq((max(Laake_PrimarySightings$X)+1):(max(Laake_PrimarySightings$X) + nrow(sightings.primary))),
#             Date = Date,
#             day = day,
#             month = month,
#             year = year,
#             watch = watch,
#             t241 = t241,
#             distance = distance,
#             podsize = podsize,
#             vis = vis,
#             beaufort = beaufort,
#             wind.direction = NA,
#             key = key,
#             pphr = pphr,
#             Start.year = Start.year,
#             original.watch = watch,
#             only = TRUE,
#             hours = NA,
#             Sex = NA,
#             Observer = Observer) -> sightings.
# 
# Sightings <- rbind(Laake_PrimarySightings, sightings.)
# 
# # Filter effort and sightings and store in dataframes Effort and Sightings
# effort.primary %>%
#   transmute(X = seq((max(Laake_PrimaryEffort$X)+1), (max(Laake_PrimaryEffort$X) + nrow(effort.primary))),
#             watch.key = watch.key,
#             Start.year = Start.year,
#             key = key,
#             begin = begin,
#             end = end,
#             npods = npods,
#             nwhales = nwhales,
#             effort = effort,
#             vis = vis,
#             beaufort = beaufort,
#             Observer = Observer,
#             time = time,
#             watch = watch,
#             Use = Use,
#             Date = Date) -> effort.
#  
# Effort <- rbind(Laake_PrimaryEffort, effort.) 
# 
# # The following line is needed even if there was no correction made. 
# Sightings$corrected.podsize = Sightings$podsize

# There was a podsize estimate of 50 on 2011-01-31 11:43:21. This value is a bit
# dubious given others are almost all less than 10, and the maximum is usually
# 20 or less. So, I changed it to 5. 

# I give up... 2023-12-20 I was unable to use 2010 and 2011 data to update the 
# ratio between naive and "true" estimates. So, I will use the analysis up to 
# 2006 by Laake et al. and apply the ratio to the new estimates from 2010 to 2023. 
# 

#library(ERAnalysis)  

#Run example code from the help file
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
# of the off-effort sightings.  The code below shows that the counts of primary
# sighting records matches the total counts of the on and off effort sightings split into
# the 2 dataframes.
#

data(PrimaryOff)
data(Primary)
data(ERSurveyData)
NorthYears=c(2000,2001)

data(PrimarySightings)
data(PrimaryEffort)

PrimarySightings=PrimarySightings[!(PrimarySightings$Observer=="MAS"&PrimarySightings$Start.year==1995),]
PrimaryEffort=PrimaryEffort[!(PrimaryEffort$Observer=="MAS"&PrimaryEffort$Start.year==1995),]

# Define arguments used in the analysis
all.years=unique(PrimaryEffort$Start.year)
recent.years=all.years[all.years>=1987]
early.years=all.years[all.years<1987]
final.time = sapply(tapply(floor(PrimaryEffort$time),
                           PrimaryEffort$Start.year,max),
                    function(x) ifelse(x>90,100,90))
lower.time=rep(0,length(final.time))
fn=1.0817
se.fn=0.0338

# These are the number of sightings that were included/excluded based on Use
Sightings=PrimarySightings
Sightings=merge(Sightings,subset(PrimaryEffort,select=c("key","Use")))

# Filter effort and sightings and store in dataframes Effort and Sightings
Effort=PrimaryEffort[PrimaryEffort$Use,]  
Sightings=PrimarySightings
Sightings$seq=1:nrow(Sightings)
Sightings=merge(Sightings,subset(Effort,select=c("key")))
Sightings=Sightings[order(Sightings$seq),]

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

# Compute simple minded population estimate and plot it
period.estimate=apply(whales.counted/sampled.fraction,1,sum,na.rm=TRUE)

# Compute naive estimates of abundance for the 23 surveys.  These use the uncorrected
# counts of whales from the primary observer during watches in which neither Beaufort nor
# vis exceeded 4.  For each year a gam with a smooth over time is fitted and this is
# used to predict total abundance throughout the migration from the counts of whales
# during the sampled periods.  There is no correction for missed pods or for measurement
# error in podsize. Each fitted migration gam is plotted with the observed values and
# saved in the file NaiveMigration.pdf.
#pdf("NaiveMigration.pdf")
#naive.abundance.models=vector("list",23)

# Run naive abundance estimates for Laake's data and new data.

all.sightings <- rbind(Sightings,  sightings.primary )

all.effort <- rbind(Effort, effort.primary)

final.time <- c(final.time, rep(90, times = length(effort.list.primary)))
lower.time <- c(lower.time, rep(0, times = length(effort.list.primary)))

naive.abundance.models=vector("list",length(c(all.years, (years))))
i=0
for (year in c(all.years, (years-1))){
  i=i+1
  primary=all.sightings[all.sightings$Start.year==year,]
  primary$Start.year=factor(primary$Start.year)
  ern=subset(all.effort,
             subset=as.character(Start.year)==year,
             select=c("Start.year","key","begin","end","effort","time","vis","beaufort"))
  
  ern$Start.year=factor(ern$Start.year)
  naive.abundance.models[[i]]=estimate.abundance(spar=NULL,
                                                 dpar=NULL,
                                                 gsS=gsS,
                                                 effort=ern, 
                                                 sightings=primary, 
                                                 final.time=final.time[i],
                                                 lower.time=lower.time[i],
                                                 gformula=~s(time),
                                                 dformula=NULL)
}
#dev.off()
#
#
# Nhat.naive=sapply(naive.abundance.models,function(x) x$Total)
# 
# # Define set of models to be evaluated for detection
# models=c("podsize+Dist+Observer",
#          "podsize+Dist+Observer+beaufort",
#          "podsize+Dist+Observer+vis",
#          "podsize+Dist+Observer+Vis")
# 
# # Create time series of estimates based on previous approach using Reilly pod size
# # correction method but using 1978 data for surveys <=1987 and 1992-1994 aerial data
# # for surveys >=1992
# data(add.cf.reilly)
# data(add.cf.laake)
# Sightings$corrected.podsize[Sightings$Start.year<=1987]=reilly.cf(Sightings$podsize[Sightings$Start.year<=1987],
#                                                                   add.cf.reilly)
# Sightings$corrected.podsize[Sightings$Start.year>1987]=reilly.cf(Sightings$podsize[Sightings$Start.year>1987],
#                                                                  add.cf.laake)
# #pdf("ReillyApproach.pdf")
# reilly.estimates=compute.series(models,
#                                 naive.abundance.models,
#                                 sightings=Sightings,
#                                 effort=Effort,TruePS=FALSE)
# 
# Nhat.reilly=reilly.estimates$Nhat
# Nhat.reilly[1:15]=Nhat.naive[1:15]*(Nhat.reilly[16]/Nhat.naive[16])
# avg.Reilly.podsize=tapply(Sightings$corrected.podsize,Sightings$Start.year,mean)
# 
# # Next compute the series of abundance estimates for the most recent 8 years by
# # fitting and selecting the best detection model but not applying the pod size correction.
# # From those 8 estimates and the naive estimates, compute an average ratio and 
# # apply it to generate the estimates for the first 15 surveys prior to 1987.
# Sightings$corrected.podsize = Sightings$podsize
# abundance.estimates.nops.correction = compute.series(models, 
#                                                      naive.abundance.models,
#                                                      sightings=Sightings,
#                                                      effort=Effort,
#                                                      TruePS=FALSE)
# Sightings$corrected.podsize=NULL

# Next compute the series of abundance estimates for the most recent 8 years by
# fitting and selecting the best detection model.  From those 8 estimates and the
# naive estimates, compute an average ratio and apply it to generate the estimates
# for the first 15 surveys prior to 1987. Note with hessian=TRUE, the analysis can
# take about 30-60 minutes to complete. (TE: Takes about 6 minutes now. But added
# the if-else. 2023-08-31)
# 

# Because these require some time to run, I saved the output in an rds file. To
# run the code, uncomment the following block:

####################### Block starts here:
# naive.abundance.models.Laake=vector("list",length(all.years))
# i=0
# for (year in all.years){
#   i=i+1
#   primary=Sightings[Sightings$Start.year==year,]
#   primary$Start.year=factor(primary$Start.year)
#   ern=subset(Effort,
#              subset=as.character(Start.year)==year,
#              select=c("Start.year","key","begin","end","effort","time","vis","beaufort"))
#   
#   ern$Start.year=factor(ern$Start.year)
#   naive.abundance.models.Laake[[i]]=estimate.abundance(spar=NULL,
#                                                        dpar=NULL,
#                                                        gsS=gsS,
#                                                        effort=ern, 
#                                                        sightings=primary, 
#                                                        final.time=final.time[i],
#                                                        lower.time=lower.time[i],
#                                                        gformula=~s(time),
#                                                        dformula=NULL)
# }
# 
# saveRDS(naive.abundance.models.Laake,
#         file = "RData/naive_abundance_Laake.rds")
####################### Block ends here

naive.abundance.models.Laake <- readRDS(file = "RData/naive_abundance_Laake.rds")

# To modify naive estimates to better estimates, I also ran the code and saved
# the results in an rds file. To run the code, uncomment the following block:

####################### Block starts here:
##### From Laake_example_code.R

#Next compute the series of abundance estimates for the most recent 8 years by
# fitting and selecting the best detection model but not applying the pod size correction.
# From those 8 estimates and the naive estimates, compute an average ratio and 
# apply it to generate the estimates for the first 15 surveys prior to 1987.
# Sightings$corrected.podsize=Sightings$podsize
# abundance.estimates.nops.correction=compute.series(models, 
#                                                    naive.abundance.models.Laake,
#                                                    sightings=Sightings,
#                                                    effort=Effort,
#                                                    TruePS=FALSE)
# saveRDS(abundance.estimates.nops.correction,
#         file = "RData/Nhats_nops_correction.rds")
# 
# Sightings$corrected.podsize=NULL
# abundance.estimates.1967to2006 = compute.series(models,
#                                                 naive.abundance.models.Laake,
#                                                 sightings=Sightings,
#                                                 effort=Effort,
#                                                 hessian=TRUE)

# saveRDS(abundance.estimates.1967to2006,
#         file = "RData/Nhats_1967to2006.rds")

####################### Block ends here

abundance.estimates.1967to2006 <- readRDS(file = "RData/Nhats_1967to2006.rds")
##### End of Laake_example_code.R 

# The ratio between "naive" and non-naive Nhats when there were no two stations
# was ~1.57 (1.567350 to be exact).
# lapply(naive.abundance.models.Laake, FUN = function(x) x$Total) %>% unlist -> naive.Nhat
# abundance.estimates.1967to2006$Nhat/naive.Nhat

naive.Nhat.2010to2024 <- lapply(naive.abundance.models.new, 
                                FUN = function(x) x$Total) %>% unlist
Nhat.2010to2024 <- naive.Nhat.2010to2024 * 1.567

ratio <- abundance.estimates.1967to2006$ratio
ratio.SE <- 0.03  # from Laake et al 2012, Table 8

# Compute series of estimates for before 1987 without nighttime correction factor (eqn 24)  
W.hat.1 <- c(sapply(naive.abundance.models[1:15], 
                    function(x)x$Total)*ratio)

# Apply nighttime correction factor (eqn 29)
fn = 1.0817
SE.fn <- 0.0338

# Need to add CI or SE and var-cov matrix. 
# Bring in the results from Laake_example_code.R
abundance.vc <- read_rds(file = "RData/abundance.vc.rds")

# Var(W.hat) for year < 1987
W.tilde.1 <- c(sapply(naive.abundance.models[1:15], function(x) x$Total))  # naive abundance
var.W.tilde.1 <- c(sapply(naive.abundance.models[1:15], function(x) x$var.Total))
var.W.hat.1 <- W.tilde.1^2 * ratio.SE^2 * 9 + ratio^2 * var.W.tilde.1   # eqn 27

N.hat.1 <- W.hat.1 * fn

# the following four lines are not quite right - see Table 9 to compare CV values 2024-01-31
# var.Nhat.1 <- (fn * W.hat.1)^2 * ((SE.fn/fn)^2 + (var.W.hat.1/((W.hat.1)^2)))  # eqn 30
# SE.Nhat.1 <- sqrt(var.Nhat.1)
# CV.Nhat.1 <- SE.Nhat.1/N.hat.1


# SE values are a little different from what I calculated above (SE.Nhat.1) but not by much
SE.Nhat.1 <- abundance.vc$se[1:length(N.hat.1)]

# var(W.hat) for year > 1985 eqn. 25
# From Table 8 in Laake et al. 2012
W.hat.2 <- setNames(c(24883, 14571, 18585, 19362, 19539, 15133, 14822, 17682),
                    c("1987", "1992", "1993", "1995", "1997", "2000", "2001", "2006"))

#W.hat <- c(W.hat.1, W.hat.2)
N.hat.2 <- abundance.estimates.1967to2006$summary.df$Nhat
SE.Nhat.2 <- abundance.vc$se[(length(N.hat.1)+1):length(abundance.vc$se)]

# The same approach for year < 1987 will be used for years 2009 - most recent
# Although the values didn't match exactly, they were quite close. So, 
# I'm not going to worry too much about it.
W.tilde.3 <- c(sapply(naive.abundance.models.new, function(x) x$Total))
var.W.tilde.3 <- c(sapply(naive.abundance.models.new, function(x) x$var.Total))
W.hat.3 <- c(sapply(naive.abundance.models.new, 
                    function(x)x$Total)*ratio)

var.W.hat.3 <- W.tilde.3^2 * ratio.SE^2 * 9 + ratio^2 * var.W.tilde.3   # eqn 27
N.hat.3 <- W.hat.3 * fn

# Fix the following three lines according to what I find on lines 721-724
var.Nhat.3 <- (fn * W.hat.3)^2 * ((SE.fn/fn)^2 + (var.W.hat.3/((W.hat.3)^2)))  # eqn 30
SE.Nhat.3 <- sqrt(var.Nhat.3)

# Function from Laake's code
conf.int=function(abundance, CV, alpha=0.05, digits=2, prt=FALSE){
  # Computes confidence intervals based on lognormal distr.
  # JMB / NMML / 11 Sep 2008
  
  if (alpha <0 || alpha > .999) stop("alpha must be in (0,1)")
  z = round(abs(qnorm(alpha/2)),2)
  if (prt) cat("N:",abundance,"  cv:",CV,"  alpha:",alpha,"  z:",z,"\n")
  C <- exp(z * sqrt(log(1 + CV^2)))
  SL <- round(abundance/C,digits)
  SU <- round(abundance * C,digits)
  data.frame(SL,SU)
}

all.estimates <- data.frame(Start.Year = c(all.years, 
                                           lapply(naive.abundance.models.new,
                                                  function(x) names(x$Total)) %>% 
                                             unlist() %>% as.numeric), 
                            Nhat = c(N.hat.1, N.hat.2, N.hat.3),
                            SE = c(SE.Nhat.1, SE.Nhat.2, SE.Nhat.3)) %>%
  mutate(CV = SE/Nhat,
         Season = paste0(Start.Year, "/", (Start.Year + 1)))

CI <- conf.int(all.estimates$Nhat, all.estimates$CV)

all.estimates$CL.low <- CI$SL
all.estimates$CL.high <- CI$SU

seq.years <- data.frame(Start.Year = seq(min(all.estimates$Start.Year), 
                                         max(all.estimates$Start.Year)))

seq.years %>% left_join(all.estimates, by = "Start.Year") %>%
  relocate(Season, .after = Start.Year) -> all.estimates

if (!file.exists(paste0("Data/all_estimates_Laake_",
                        max(years), ".csv")))
  write.csv(all.estimates, paste0("Data/all_estimates_Laake_",
                                  max(years), ".csv"), 
            quote = FALSE, row.names = FALSE)

ggplot(all.estimates) +
  geom_point(aes(x = Start.Year, y = Nhat)) +
  geom_errorbar(aes(x = Start.Year, ymin = CL.low, ymax = CL.high))

# Debugging purposes  #######################################
# models <- models 
# naive.abundance <- naive.abundance.models 
# sightings <- Sightings
# effort <- Effort
# gsS=NULL
# best=TRUE
# Match=Match.all
# cutoff=4
# final.time=c(rep(90,20),100,100,90)
# lower.time=rep(0,23) 
# lcut= -0.05 # 0.2 #  without linking sightings, the function does not run
# mcut=1.0
# twt= 0.01 #0.18 in minutes - make is short to keep all?
# dwt= 3.95
# pwt=0.05
# crittype="ellipse"
# DistBreaks=c(0,1,2,3,4,20)
# Use=TRUE
# hessian=FALSE
# debug=FALSE
# TruePS=FALSE
# recent.years=c(1987,1992,1993,1995,1997,2000,2001,2006,2009, 2010)
# fn=1.0817

# Then go into compute.series.new.r and run each line separately
#############################################################


# 
# data.fit <- sightings.all %>%
#   select(Start.year, podsize, distance, vis, beaufort, Observer) %>%
#   mutate(seen = ifelse(podsize > ))
# 
# initial.models=vector("list",length(recent.years))
# i=0
# for (year in recent.years){
#   i=i+1
#   zz=data.fit[data.fit$Start.year==year, ]
#   zz$Start.year=factor(zz$Start.year)
#   zz$Dist=cut(zz$distance,DistBreaks)
#   zz$Observer=factor(zz$Observer)
#   zz$Vis=cut(zz$vis,c(0,3,6))
#   initial.models[[i]]=select.detection.models(zz,models,cutoff)
#   print(summary(initial.models[[i]][[1]]))
#   cat("\n",length(initial.models[[i]]))
# }
# 
# 
# 
# 
# sightings.all$corrected.podsize = sightings.all$podsize
# abundance.estimates.nops.correction=compute.series.new(models, 
#                                                        naive.abundance.models,
#                                                        sightings=sightings.all,
#                                                        effort=effort.all,
#                                                        TruePS=FALSE)

