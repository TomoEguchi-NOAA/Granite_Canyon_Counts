# LaakeData2WinBUGS
# 
# Creates WinBUGS input data from Laake's ERAnalysis library and runs WinBUGS.
# 
# Code chunks with var1 = expressions are by Laake. I use var1 <- expressions

# I can't seem to get this one to work. The error is usually value of bernoulli z[1,1] must be an integer
# 
# I even rewrote the code to make the number of stations to be 1, explicitly. Don't know what to do here. 2023-09-12

# I have not run this script on the entire dataset. This should be done to compare the two approaches. 2024-03-14

rm(list = ls())

library(ERAnalysis)
library(tidyverse)
library(ggplot2)
library(R2WinBUGS)

# This brings in WinBUGS inputs that actually worked. I can compare all the objects
# to see what's different in Laake's data.
#run.date <- "2022-10-21"   # the date that WinBUGS Ver2.Rmd was run and the output saved.
#data.worked <- readRDS(paste0("RData/BUGS_data_runs_",
#                              run.date, ".rds"))

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

# Estimates from Laake et al. are here:
# col.defs <- cols(Year = col_character(),
#                  Nhat = col_double(),
#                  CV = col_double())
# 
# Laake.estimates <- read_csv(file = "Data/Laake et al 2012 Table 9 Nhats.csv",
#                             col_types = col.defs) %>% 
#   mutate(SE = CV * Nhat,
#          LCL = Nhat - 1.96 * SE,
#          UCL = Nhat + 1.96 * SE,
#          Season = lapply(strsplit(Year, "_"), 
#                          FUN = function(x) paste0(x[1], "/", x[2])) %>% 
#            unlist) %>%
#   dplyr::select(Season, Nhat, SE, LCL, UCL)  %>%
#   mutate(Year = lapply(str_split(Season, "/"), 
#                        FUN = function(x) x[2]) %>% 
#            unlist() %>% 
#            as.numeric())

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
# There is no need for sightings in Durban's approach
Effort.1 = PrimaryEffort[PrimaryEffort$Use,]  

# Sightings = PrimarySightings
# Sightings$seq = 1:nrow(Sightings)
# Sightings = merge(Sightings, subset(Effort, select=c("key")))
# Sightings = Sightings[order(Sightings$seq),]

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

# sightings and efort
sightings.list.primary <- effort.list.primary  <- list()
sightings.list.secondary <- effort.list.secondary  <- list()
k <- 1
for (k in 1:length(years)){
  # These raw data files contain repeated observations of all groups.
  # 
  
  # tmp.sightings <- read.csv(paste0("Data/all_sightings_", 
  #                                  years[k], "_Tomo_v2.csv")) 
  # 
  # tmp.sightings %>%
  #   mutate(Date1 = as.Date(Date, format = "%m/%d/%Y")) %>%
  #   #group_by(group) %>%
  #   transmute(Date = Date1,
  #             Time = Time_PST, 
  #             day = day(Date),
  #             month = month(Date),
  #             year = year(Date),
  #             watch = shift,
  #             t241 = difftime((paste(Date, Time)),
  #                             (paste0((years[k] - 1), 
  #                                     "-11-30 00:00:00"))) %>%
  #               as.numeric(),
  #             Group_ID = Group_ID,
  #             distance = Distance,
  #             podsize = nwhales,
  #             vis = vis,
  #             beaufort = beaufort,
  #             Start.year = years[k]-1,
  #             Observer = toupper(observer),
  #             key = key,
  #             pphr = pphr,
  #             date.shift = date.shift,
  #             station = station) %>%
  #   arrange(Date, Group_ID) %>%
  #   filter(vis < 5, beaufort < 5) %>%
  #   na.omit() -> sightings.all
  # 
  # sightings.all %>% 
  #   filter(station == "P") -> sightings.list.primary[[k]]
  # 
  # sightings.all %>% 
  #   filter(station == "S") -> sightings.list.secondary[[k]]
  
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
# sightings.primary <- do.call("rbind", sightings.list.primary) %>%
#   na.omit()
# 
# sightings.primary %>% 
#   left_join(all.observers, by = "Observer") %>%
#   select(-c(Observer, Initials)) %>%
#   rename(Observer = ID) -> sightings.primary
# 
# sightings.secondary <- do.call("rbind", sightings.list.secondary) %>%
#   na.omit()
# 
# sightings.secondary %>% 
#   left_join(all.observers, by = "Observer") %>%
#   select(-c(Observer, Initials)) %>%
#   rename(Observer = ID) -> sightings.secondary

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
         shift = shift.definition.frac.hr(start.hr)) %>%
  select(Start.year, nwhales, effort, vis, beaufort, obs, dt, begin, start.hr, shift) %>%
  group_by(Start.year) %>%
  mutate(effort.min = effort * 24 * 60,
         watch.prop = effort.min/540) -> Effort.1.by.period

Effort.2 %>% 
  mutate(Day1 = as.Date(paste0(Start.year, "-12-01")),
         dt = as.numeric(as.Date(Date) - Day1) + 1,
         obs = Observer,
         start.hr = frac.day2time(begin),
         shift = shift.definition.frac.hr(start.hr)) %>%
  select(Start.year, nwhales, effort, vis, beaufort, obs, dt, begin, start.hr, shift) %>%
  group_by(Start.year) %>%
  mutate(effort.min = effort * 24 * 60,
         watch.prop = effort.min/540) -> Effort.2.by.period

# Need to combine old and new effort dataframes, for primary and secondary
 


#Effort.by.day %>%
Effort.1.by.period %>%
  mutate(Initials = obs) %>%
  left_join(Observer, by = "Initials") %>%
  dplyr::select(-c(Initials, Observer, Name, Sex)) %>%
  rename(ID.1 = ID) %>%
  mutate(ID.char = obs) %>%
  left_join(Observer.1, by = "ID.char") %>%
  dplyr::select(-c(Initials, Observer, Name, Sex, ID.char)) -> Effort.1.by.period.1

#Effort.by.day.1$ID[is.na(Effort.by.day.1$ID)] <- Effort.by.day.1$ID.1[is.na(Effort.by.day.1$ID)]
Effort.by.period.1$ID[is.na(Effort.by.period.1$ID)] <- Effort.by.period.1$ID.1[is.na(Effort.by.period.1$ID)]



#############################################
# Also run WinBugs and see how that compares.

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
    
    BUGS.n[1:nrow(tmp), k] <- tmp$n + 1
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

BUGS.data.1 <- create.WinBUGS.data(Effort.by.period.1)



MCMC.params <- list(n.iter = 125000,
                    n.thin = 50,
                    n.burnin = 25000,
                    n.chains = 5)

# Bring in the most recent output from WiBUGS Ver2.Rmd
# this file contains all necessary inputs for 2006 - 2019:
data.0 <- readRDS("RData/2006-2019_GC_Formatted_Data.RDS")

all.years <- c(2007, 2008, 2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024)

# output from Ver2.0 extraction
years <- c("2015", "2016", "2020", "2022", "2023", "2024")

#out.v2 <- lapply(years, FUN = function(x) readRDS(paste0("RData/V2.1_Aug2022/out_", x, "_Tomo_v2.rds")))

# 85 or 30
min.duration <- 85 #30

out.v2 <- lapply(years, 
                 FUN = function(x) readRDS(paste0("RData/V2.1_Feb2024/out_", x,
                                                  "_min", min.duration, "_Tomo_v2.rds")))

# just for 8 weeks in 2023
#out.v2 <- lapply(years, FUN = function(x) readRDS(paste0("RData/V2.1_Mar2023_8weeks/out_", x, 
#                                                         "_min", min.duration, "_Tomo_v2.rds")))

begin. <- lapply(out.v2, FUN = function(x) x$Final_Data$begin)
end. <- lapply(out.v2, FUN = function(x) x$Final_Data$end)

# Number of watch periods in each year's survey - before the 2019/2020 season
# plus the new ones
# This has been changed for 2024. I now have edited data for 2010 - 2024 seasons.
# So, I can just use the first two (2006/2007 and 2007/2008). The two numbers for 
# 2009/2010 and 2010/2011 don't match. In Durban's analysis, they were 164 and 178,
# respectively. For the new extraction, they are 180 and 145. I will stick with
# Durban's data for now. 2024-03-14
periods <-c(136, 135, 164, 178,
            lapply(begin., FUN = function(x) length(x)) %>% unlist)

# Shorten data to first x years only to replicate 
# the analysis in Durban et al 2016. 
# Or use it for other purposes.
x <- length(periods)

out.file.name <- paste0("RData/WinBUGS_", x, "yr_v2_min", min.duration, ".rds")

# Bring in the most recent WinBUGS run results, which contain input data, initial values
# and other essentials
Ver2.results <- readRDS(out.file.name)




# In WinBUGS code, the number of days per season is fixed at 90. So, without
# adjusting the code, which I don't want to do for comparing results, the data
# needs to be adjusted so that there are no days > 90.
WinBUGS.dir <- paste0(Sys.getenv("HOME"), "/WinBUGS14")

# For Durban's WinBUGS model, each sampling period is treated as is, rather than 
# grouping them by day.

Effort %>%
  mutate(Day1 = as.Date(paste0(Start.year, "-12-01")),
         dt = as.numeric(as.Date(Date) - Day1) + 1) %>%
  select(Start.year, nwhales, effort, vis, beaufort, Observer, dt) -> Durban.data

Durban.data %>%
  mutate(Initials = Observer) %>%
  left_join(Observer, by = "Initials") %>%
  mutate(obs = Observer.x) %>%
  dplyr::select(-c(Initials, Observer.x, Observer.y, Name, Sex)) %>%
  rename(ID.1 = ID) %>%
  mutate(ID.char = obs) %>%
  left_join(Observer.1, by = "ID.char") %>% #-> tmp
  dplyr::select(-c(Initials, Observer, Name, Sex, ID.char)) -> Durban.data.1

Durban.data.1$ID[is.na(Durban.data.1$ID)] <- Durban.data.1$ID.1[is.na(Durban.data.1$ID)]

# Also shrunk the dataset to troubleshoot... no avail. 
Durban.data.1 %>% 
  filter(dt < 90) -> Durban.data.2
  #filter(dt < 90, Start.year > 1985) -> Durban.data.2

BUGS.data <- create.WinBUGS.data(Durban.data.2)

x <- length(BUGS.data$periods)

#we're going to make N a partially observed data object with anchor points at day 1 and 90
# TE: I don't know how these numbers were created... they are generally 2x n (not all)
# N_inits <- as.matrix(read.table("Data/Initial Values/N_inits.txt",
#                                 header=T))

# For this run, there is no 2nd station. 
n <- BUGS.data$n
for (k in 1:ncol(n)){
  if (BUGS.data$periods[k] < nrow(n))
    n[(BUGS.data$periods[k]+1):nrow(n),k] <- NA
}
N_inits1 <- n * 2 + 2
#N_inits2 <- BUGS.data$n[, 2,] * 2 + 2 

N_inits <- N_inits1
#N_inits[N_inits1 < N_inits2] <- N_inits2[N_inits1 < N_inits2]

N_inits <- rbind(N_inits,
                 matrix(data = NA, 
                        nrow = 2, 
                        ncol = length(BUGS.data$periods)))

# for (k in 1:length(BUGS.data$periods)){
#   N_inits[(BUGS.data$periods[k]+1):nrow(N_inits), k] <- NA  
# }

# NAs for all Ns that are estimated, and 0s for the days 1 and 90
N <- matrix(data = NA, 
            nrow=max(BUGS.data$periods)+2, 
            ncol=length(BUGS.data$periods)) 

for(i in 1:length(BUGS.data$periods)){
  #True number of whales passing fixed at 0 for day 1 and 90
  N[(BUGS.data$periods[i]+1):(BUGS.data$periods[i]+2), i] <- 0 
}

BUGS.data$N <- N
BUGS.data$N.com <- N
BUGS.data$N.sp <- N

BUGS.data$knot <-  c(-1.46,-1.26,-1.02,-0.78,
                     -0.58,-0.34,-0.10,0.10,
                     0.34,0.57,0.78,1.02,1.26,1.46)

BUGS.data$n.knots <- 14

# the u data is whether there were observers on watch. 
# 0 counts are often associated with years/shifts with 
# no second observer. So if u=0, it will fix observation probability at 0
# the second column for each year is for the second station - not the second
# observer.
# 
# For this comparison, I only use the primary observer, so I can just have 1s
# in the primary observer place, and zeros elsewhere
u <- array(data = 0, dim = dim(BUGS.data$n))

for (y in 1:BUGS.data$n.year){
  for (p in 1:BUGS.data$periods[y]){
    #for (s in 1:BUGS.data$n.station){
      u[p,y] <- 1
    #}
  }
}

BUGS.data$u <- u

BUGS.inits <- function() list(mean.prob = 0.5,
                               BF.Fixed = 0,
                               VS.Fixed = 0,
                               mean.prob.sp = 0.5,
                               BF.Fixed.sp = 0,
                               VS.Fixed.sp = 0,
                               mean.prob.com = 0.5,
                               BF.Fixed.com = 0,
                               VS.Fixed.com = 0,
                               mean.beta = c(0,0,0), #mean.beta = c(5,0.14,-3.5),
                               beta.sigma = c(1,1,1),#beta.sigma = c(7,7,7),
                               BF.Switch = 1,
                               VS.Switch = 1,
                               OBS.Switch = 1,
                               sigma.Obs = 1,
                               BF.Switch.sp = 1,
                               VS.Switch.sp = 1,
                               OBS.Switch.sp = 1,
                               sigma.Obs.sp = 1,
                               BF.Switch.com = 1,
                               VS.Switch.com = 1,
                               OBS.Switch.com = 1,
                               sigma.Obs.com = 1,
                               N = N_inits,
                               N.com = N_inits,
                               N.sp = N_inits,
                               beta.sp = array(data=0, dim=c(2,x)),
                               sd.b.sp = rep(1, times = x),
                               z = matrix(1, nrow=90, ncol= x)) #)

## Compare between the data that worked and the new one that does not work
#BUGS.data.worked <- data.worked$BUGS.data
# 
# BUGS.inits1 <- BUGS.inits()
#BUGS.inits.worked <- data.worked$BUGS.inits
# 
# 
# BUGS.data$n[,1,2]
# BUGS.data.worked$n[,1,2]

## End of comparison


BUGS.parameters <- c("lambda","OBS.RF","OBS.Switch",
                     "BF.Switch","BF.Fixed","VS.Switch",
                     "VS.Fixed","mean.prob","mean.prob.com",
                     "mean.prob.sp","BF.Fixed.com",
                     "BF.Fixed.sp","VS.Fixed.com",
                     "VS.Fixed.sp",
                     "Corrected.Est","Raw.Est","z",
                     "com","sp","Daily.Est","mean.beta",
                     "beta.sigma","beta","beta.sp","b.sp","sd.b.sp")

out.file.name <- paste0("RData/WinBUGS_Laake_Data.rds")

if (!file.exists(out.file.name)){
  
  #Run time: 
  Start_Time<-Sys.time()
  
  BUGS_out <- bugs(data = BUGS.data,
                   inits = BUGS.inits,
                   parameters = BUGS.parameters,
                   model.file="GW_Nmix_Orig.bugs",
                   n.chains = MCMC.params$n.chains,
                   n.iter = MCMC.params$n.iter, #MCMC.params$n.samples, 
                   n.burnin = MCMC.params$n.burnin, 
                   n.thin = MCMC.params$n.thin,
                   debug=F,
                   bugs.directory = WinBUGS.dir,
                   DIC = FALSE)
  
  Run_Time <- Sys.time() - Start_Time
  BUGS.out <-list(BUGS.data = BUGS.data,
                  N_inits = N_inits,
                  BUGS_out = BUGS_out,
                  Run_Time = Run_Time,
                  Sys.info = Sys.info())
  
  saveRDS(BUGS.out, file = out.file.name)
  
} else {
  BUGS.out <- readRDS(out.file.name)
}

# 2024-02-15: Compare BUGS.data and code to find discrepancies in input array or vector
# sizes vs. what's in the code.

# 2023-03-03 (problem solved but leave the comments below for future reference)
# ERROR: NIL dereference (read). According to the user manual, this error may
# happen "at compilation in some circumstances when an inappropriate
# transformation is made, for example an array into a scalar." 
# (https://www.mrc-bsu.cam.ac.uk/wp-content/uploads/manual14.pdf)
# https://stackoverflow.com/questions/21969600/bugs-error-messages

# DIC problems: Surprisingly, sometimes when getting a trap (including one with the very
# informative title “NIL dereference (read)”), setting the argument DIC = FALSE in the bugs()
# function has helped. (https://www.mbr-pwrc.usgs.gov/software/kerybook/AppendixA_list_of_WinBUGS_tricks.pdf)

# I only changed the effort (Watch.Length) for this analysis using 30 minutes as the 
# minimum observation duration. That changed the size of data arrays, which shouldn't be an issue. 

# It turned out a new observer (JSD) was introduced when a shorter minimum was used to filter
# observation period. "JSD" was only found in 2020 (200109_080043). A strange thing happened for that day.
# During the 0800 shift, JSD/JWG changed to SJC/JDS on 0900, then the shift continued until 0930. The new
# cutoff time (30 min) picked up the first (1 hr) and the second (30 min) as separate observation periods.


BUGS.out$BUGS_out$summary %>% 
  as.data.frame() %>%
  rownames_to_column(var = "Parameter") -> Summary.stats

Corrected.Est <- Summary.stats[grep("Corrected", Summary.stats$Parameter),]

Nhat.df <- data.frame(median = Corrected.Est$`50%`,
                      LCL = Corrected.Est$`2.5%`,
                      UCL = Corrected.Est$`97.5%`,
                      mean = Corrected.Est$mean,
                      Season =  Laake.estimates$Season)

Nhat.df %>% 
  left_join(Laake.estimates, by = "Season") %>%
  mutate(delta_Nhat = median - Nhat) -> all.estimates

ggplot(all.estimates) +
  geom_point(aes(x = Season, y = median ),
             color = "blue") +
  geom_errorbar(aes(x = Season, ymin = LCL.x, ymax = UCL.x),
                color = "blue") +
  geom_point(aes(x = Season, y = Nhat ),
             color = "green") +
  geom_errorbar(aes(x = Season, ymin = LCL.y, ymax = UCL.y),
                color = "green")


                    