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
  na.omit()

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
              Observer = observer,
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

sightings %>% left_join(all.observers, by = "Observer") %>%
  select(-c(Observer, Initials)) %>%
  rename(Observer = ID) -> sightings

# Effort
effort <- do.call("rbind", effort.list)  %>%
  na.omit()

effort %>% left_join(all.observers, by = "Observer") %>%
  select(-c(Observer, Initials)) %>%
  rename(Observer = ID) -> effort

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
  transmute(X = (max(Laake_PrimarySightings$X) + 1) : (max(Laake_PrimarySightings$X) + 1 + nrow(sightings) - 1),
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

effort.Laake.format <- data.frame(X = (max(Laake_PrimaryEffort$X) + 1) : (max(Laake_PrimaryEffort$X) + 1 + nrow(effort) - 1),
                                  watch.key = effort$watch.key,
                                  Start.year = effort$Start.year,
                                  key = effort$key,
                                  begin = effort$begin,
                                  end = effort$end,
                                  npods = effort$npods,
                                  nwhales = effort$nwhales,
                                  effort = effort$effort,
                                  vis = effort$vis,
                                  beaufort = effort$beaufort,
                                  Observer = effort$Observer,
                                  time = effort$time, 
                                  watch = effort$watch,
                                  Use = effort$Use,
                                  Date = effort$Date)

effort.all <- rbind(Laake_PrimaryEffort, effort.Laake.format)


final.time = sapply(tapply(floor(effort.all$time), effort.all$Start.year, max), function(x) ifelse(x>90,100,90))
lower.time = rep(0,length(final.time))

all.years <- unique(sightings.all$Start.year)
naive.abundance.models=vector("list", length(all.years))
i <- 0
for (year in all.years)
{
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

# Needs functions in the following script
source("~/R/ERAnalysis/R/MatchingLinking.r")
source("compute.series.new.r")
source("~/R/ERAnalysis/R/io.glm.R")

select.detection.models=function(x,models,cutoff=4){
  mod=vector("list",length(models))      
  for(i in 1:length(models))
    mod[[i]]=io.glm(x,as.formula(paste("seen~",models[i],sep="")))
  AICValues=sapply(mod,function(x) AIC(x))
  DeltaAIC=AICValues-min(AICValues)
  mod.numbers=(1:length(models))[order(DeltaAIC)]
  return(mod[mod.numbers[sort(DeltaAIC)<cutoff]])
}

#Next compute the series of abundance estimates for the most recent 8 years plus 
# most recent 5 years (2015, 2016, 2020, 2022, 2023) by
# fitting and selecting the best detection model but not applying the pod size correction.
# From those 8 estimates and the naive estimates, compute an average ratio and 
# apply it to generate the estimates for the first 15 surveys prior to 1987.
# I created a new function compute.series.new by adding the recent 5 years.
recent.years <- c(1987,1992,1993,1995,1997,2000,2001,2006,2014,2015,2019,2021,2022)
DistBreaks=c(0,1,2,3,4,20)
cutoff=4

# To fit binomial models, I need to have the response variable, which is 0 or 1,
# I don't see such line in Laake's code. Sightings doesn't contain non-sighting
# data. So, effort data frame needs to be joined with sightings, which Laake
# does it in Line 40 of compute.series.r. But... that doesn't work... 
# 2023-11-30 

tmp <- merge(sightings.all, subset(effort.all,
                                   select=c("key","Use")), by="key")

data.fit <- sightings.all %>%
  select(Start.year, podsize, distance, vis, beaufort, Observer) %>%
  mutate(seen = ifelse(podsize > ))

initial.models=vector("list",length(recent.years))
i=0
for (year in recent.years){
  i=i+1
  zz=data.fit[data.fit$Start.year==year, ]
  zz$Start.year=factor(zz$Start.year)
  zz$Dist=cut(zz$distance,DistBreaks)
  zz$Observer=factor(zz$Observer)
  zz$Vis=cut(zz$vis,c(0,3,6))
  initial.models[[i]]=select.detection.models(zz,models,cutoff)
  print(summary(initial.models[[i]][[1]]))
  cat("\n",length(initial.models[[i]]))
}




sightings.all$corrected.podsize = sightings.all$podsize
abundance.estimates.nops.correction=compute.series.new(models, 
                                                       naive.abundance.models,
                                                       sightings=sightings.all,
                                                       effort=effort.all,
                                                       TruePS=FALSE)

