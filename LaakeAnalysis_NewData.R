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

# Observers need to be in integers, not their initials.
data("Observer")   # from ERAnalysis package.
#Observer$Set = "old"

new.observers <- read.csv(file = "Data/ObserverList2023.csv") %>%
  transmute(ID = ID,
            Initials = obs)

Observer %>%
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

# Only 2010 and 2011 had two observation stations, which are needed for applying
# Laake's method beyond naive estimates

years <- c(2010, 2011, 2015, 2016, 2020, 2022, 2023)

#years <- 2015
# sightings and efort
sightings.list.primary <- effort.list.primary  <- list()
sightings.list.secondary <- effort.list.secondary  <- list()
k <- 1
for (k in 1:length(years)){
  # These raw data files contain repeated observations of all groups.
  # 
   
  tmp.sightings <- read.csv(paste0("Data/all_sightings_", 
                                   years[k], "_Tomo_v2.csv")) 
  
  tmp.sightings %>%
    mutate(Date1 = as.Date(Date, format = "%m/%d/%Y")) %>%
    #group_by(group) %>%
    transmute(Date = Date1,
              Time = Time_PST, 
              day = day(Date),
              month = month(Date),
              year = year(Date),
              watch = shift,
              t241 = difftime((paste(Date, Time)),
                              (paste0((years[k] - 1), 
                                             "-11-30 00:00:00"))) %>%
                as.numeric(),
              Group_ID = Group_ID,
              distance = Distance,
              podsize = nwhales,
              vis = vis,
              beaufort = beaufort,
              Start.year = years[k]-1,
              Observer = toupper(observer),
              key = key,
              date.shift = date.shift,
              station = station) %>%
    arrange(Date, Group_ID) %>%
    filter(vis < 5, beaufort < 5) %>%
    na.omit() -> sightings.all

  sightings.all %>% 
    filter(station == "P") -> sightings.list.primary[[k]]
  
  sightings.all %>% 
    filter(station == "S") -> sightings.list.secondary[[k]]
  
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

# sightings
sightings.primary <- do.call("rbind", sightings.list.primary) %>%
  na.omit()

sightings.primary %>% 
  left_join(all.observers, by = "Observer") %>%
  select(-c(Observer, Initials)) %>%
  rename(Observer = ID) -> sightings.primary

sightings.secondary <- do.call("rbind", sightings.list.secondary) %>%
  na.omit()

sightings.secondary %>% 
  left_join(all.observers, by = "Observer") %>%
  select(-c(Observer, Initials)) %>%
  rename(Observer = ID) -> sightings.secondary

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

# gsS: nmax x nmax pod size calibration matrix; each row is a true pod size 
# from 1 to nmax and the value for each column is the probability that a pod of 
# a true size S is recorded as a size s (1..nmax columns)
# 
naive.abundance.models <- list()
for (k in 1:length(years)){
  sightings.primary %>%
    filter(Start.year == (years[k]-1)) -> sightings
  
  effort.primary %>%
    filter(Start.year == (years[k]-1)) -> effort
  
  naive.abundance.models[[k]] <- estimate.abundance(spar = NULL,
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
# two stations. Here, I compute those factors for the 2009/2010 and 2010/2011 seasons, which were
# not part of Laake's analysis. 

sightings.primary %>%
  transmute(X = (max(Laake_PrimarySightings$X) + 1) : (max(Laake_PrimarySightings$X) + 1 + nrow(sightings.primary) - 1),
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

effort.primary %>%
  mutate(X = (max(Laake_PrimaryEffort$X) + 1) : (max(Laake_PrimaryEffort$X) + 1 + nrow(effort.primary) - 1)) %>%
  select(-c(station)) %>%
  relocate(Observer, .after = beaufort) -> effort.Laake.format
  
effort.all <- rbind(Laake_PrimaryEffort, effort.Laake.format)

final.time = sapply(tapply(floor(effort.all$time), effort.all$Start.year, max), function(x) ifelse(x>90,100,90))
lower.time = rep(0,length(final.time))
 
all.years <- unique(sightings.all$Start.year)
naive.abundance.models=vector("list", length(all.years))
i <- 0
for (year in all.years){
  i <- i+1
  primary <- sightings.all[sightings.all$Start.year==year,]
  primary$Start.year=factor(primary$Start.year)
  ern=subset(effort.all,
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

# Define set of models to be evaluated for detection
models=c("podsize+Dist+Observer",
         "podsize+Dist+Observer+beaufort",
         "podsize+Dist+Observer+vis",
         "podsize+Dist+Observer+Vis")

# Needs functions in the following scripts
source("~/R/ERAnalysis/R/MatchingLinking.r")
source("compute.series.new.r")
source("~/R/ERAnalysis/R/io.glm.R")

# compute.series (or compute.series.new) can take a Match dataframe, which
# contains information about sighting and non-sighting of each whale pod by
# primary and secondary observers. The match dataframe from Laake's analysis
# was obtained from running their code and saving the output Match to a .csv
# file. A similar dataframe needs to be created for 2010 and 2011.

Match.Laake <- read.csv("Data/match_1987_2006.csv")

# for new datasets, only 2010 and 2011 have the secondary observers:
# need to create the "seen" variable, which is 0/1
sightings.primary %>%
  filter(Start.year == 2009 | Start.year == 2010)  %>%
  mutate(seen = 1,
         new.key = paste0(Date, "-", Group_ID)) -> sightings.primary.1

sightings.secondary %>%
  filter(Start.year == 2009 | Start.year == 2010)  %>%
  mutate(seen = 1,
         new.key = paste0(Date, "-", Group_ID)) -> sightings.secondary.1

# Need to reduce the primary sightings to match secondary survey dates:
unique.secondary.dates <- unique(sightings.secondary.1$Date)
filtered.primary.sightings.1 <- sightings.primary.1[sightings.primary.1$Date %in% unique.secondary.dates, ]
filtered.primary.effort <- effort.primary[effort.primary$Date %in% unique.secondary.dates, ]

#START FROM HERE! 2023-12-12

# Combine with the secondary sightings
rbind(filtered.primary.sightings.1, sightings.secondary.1) %>%
  arrange(by = new.key) -> sightings.all.1

new.key <- unique(c(sightings.primary$new.key, sightings.secondary$new.key))

new.key.df <- data.frame(new.key = new.key)
         
new.key.df %>%
  left_join(sightings.primary, by = "new.key") -> sightings.primary.new.key

sightings.primary.new.key$seen[is.na(sightings.primary.new.key$seen)] <- 0

# Fill in data when "seen" is 0
sightings.primary.new.key %>%
  filter(seen == 0) %>%
  mutate(Date = new.key %>% str_sub(start = 1, end = 10),
         Time = NA,
         day = new.key %>% str_sub(start = 9, end = 10) %>% as.numeric(),
         month = new.key %>% str_sub(start = 6, end = 7)%>% as.numeric(),
         year = new.key %>% str_sub(start = 1, end = 4)%>% as.numeric(),
         watch = NA,
         t241 = NA,
         Group_ID = new.key %>% str_sub(start = 12, end = 13)%>% as.numeric(),
         distance = NA, podsize = NA, vis = NA, beaufort = NA,
         Start.year = ifelse(month > 10, year, year-1),
         key = NA,
         station = "P",
         Observer = NA,
         seen = seen) -> tmp

sightings.primary.all <- rbind(sightings.primary.new.key %>% filter(seen == 1),
                               tmp) %>%
  arrange(by = new.key)

# Because there were less effort for the secondary team than the primary team,
# I match primary to secondary.
new.key.df %>%
  left_join(sightings.secondary, by = "new.key") -> sightings.secondary.new.key

sightings.secondary.new.key$seen[is.na(sightings.secondary.new.key$seen)] <- 0

# Fill in data when "seen" is 0
sightings.secondary.new.key %>%
  filter(seen == 0) %>%
  mutate(Date = new.key %>% str_sub(start = 1, end = 10),
         Time = NA,
         day = new.key %>% str_sub(start = 9, end = 10)%>% as.numeric(),
         month = new.key %>% str_sub(start = 6, end = 7)%>% as.numeric(),
         year = new.key %>% str_sub(start = 1, end = 4)%>% as.numeric(),
         watch = NA,
         t241 = NA,
         Group_ID = new.key %>% str_sub(start = 12, end = 13)%>% as.numeric(),
         distance = NA, podsize = NA, vis = NA, beaufort = NA,
         Start.year = ifelse(month > 10, year, year-1),
         key = NA,
         station = "S",Observer = NA,
         seen = seen) -> tmp

sightings.secondary.all <- rbind(sightings.secondary.new.key %>% filter(seen == 1),
                               tmp) %>%
  arrange(by = new.key)

sightings.all <- rbind(sightings.primary.all, sightings.secondary.all) %>%
  arrange(by = new.key)

# fill in watch and observer for not-seen pods.
sightings.primary.new.key %>% 
  filter(seen == 1) %>%
  select(Date, Observer, watch) -> primary.observers

sightings.secondary.new.key %>%
  filter(seen == 1) %>%
  select(Date, Observer, watch) -> secondary.observers


seen.0 <- which(sightings.all$seen == 0)
k <- 1
for (k in 1:length(seen.0)){
  tmp.0 <- sightings.all[seen.0[k],]
  tmp <- sightings.all %>% filter(new.key == tmp.0$new.key)
  sightings.all[seen.0[k], "watch"] <- na.omit(tmp$watch)
  sightings.all[seen.0[k], "key"] <- na.omit(tmp$key)
  
  if (tmp.0$station == "S"){   # if secondary is missing
    secondary.observers %>%
      filter(Date == tmp.0$Date,
             watch == sightings.all[seen.0[k], "watch"]) %>%
      select(Observer) %>%
      pull() -> obs
    
    sightings.all[seen.0[k], "watch"] <- ifelse(length(obs) == 0,
                                                NA, obs)
  } else {
    primary.observers %>%
      filter(Date == tmp.0$Date,
             watch == sightings.all[seen.0[k], "watch"]) %>%
      select(Observer) %>%
      pull() -> obs
    
    sightings.all[seen.0[k], "watch"] <- ifelse(length(obs) == 0,
                                                NA, obs)
    
  }
  
}

         mutate(wind.direction = NA,
         seq = seq(1:nrow(sightings.primary)),
         hours = NA,
         pphr = NA,
         Sex = NA,
         Use = T) %>%
  
  select(c(key, Date, seen, station, day, watch, t241, distance, podsize, vis, beaufort, wind.direction,
           Start.year, seq, hours, pphr, Sex, Observer, Use)) -> Match.primary


         wind.direction = NA,
         seq = seq(1:nrow(sightings.primary)),
         hours = NA,
         pphr = NA,
         Sex = NA,
         Use = T) %>%
  select(c(key, Date, seen, station, day, watch, t241, distance, podsize, vis, beaufort, wind.direction,
           Start.year, seq, hours, pphr, Sex, Observer, Use)) -> Match.secondary




# This function is in compute.series.r (or compute.series.new.r)
# select.detection.models=function(x, models, cutoff=4){
#   mod=vector("list",length(models))
#   for(i in 1:length(models))
#     mod[[i]]=io.glm(x,as.formula(paste("seen~",models[i],sep="")))
# 
#   AICValues=sapply(mod,function(x) AIC(x))
#   DeltaAIC=AICValues-min(AICValues)
#   mod.numbers=(1:length(models))[order(DeltaAIC)]
#   return(mod[mod.numbers[sort(DeltaAIC)<cutoff]])
# }

# Next compute the series of abundance estimates for 8 years plus 
# 2 years for which secondary observers exist (2010, 2011) by
# fitting and selecting the best detection model but not applying the pod size correction.
# From those 10 estimates and the naive estimates, compute an average ratio and 
# apply it to generate the estimates for the first 15 surveys prior to 1987 and the most
# recent 5 years (2015, 2016, 2020, 2022, 2023).
# years <- c(2010, 2011, 2015, 2016, 2020, 2022, 2023)
# # I created a new function compute.series.new by adding 2010 and 2011.
# DistBreaks=c(0,1,2,3,4,20)
# cutoff=4
# 

Sightings <- merge(sightings.primary, subset(effort.primary,
                                             select=c("key","Use")), by="key")

# Filter effort and sightings and store in dataframes Effort and Sightings
Effort <- effort.primary 

# Sightings$seq=1:nrow(Sightings)
# Sightings=merge(Sightings,subset(Effort,select=c("key")))
# Sightings=Sightings[order(Sightings$seq),]

Sightings$corrected.podsize=Sightings$podsize
abundance.estimates.nops.correction=compute.series.new(models, 
                                                       naive.abundance.models,
                                                       sightings=Sightings,
                                                       effort=Effort,
                                                       TruePS=FALSE)

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

