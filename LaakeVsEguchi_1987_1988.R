# LaakeVsEguchi_1987-1988.R
# 

# Figure out why the 1987/1988 season has very different estimates between Laake's and my estimates
# 

rm(list = ls())
source("Granite_Canyon_Counts_fcns.R")
source("Laake_functions.R")

library(tidyverse)

# Data from ERAnalysis is in the following:
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
# I use the example code from the ERAnalysis package. 
data(PrimaryOff)
data(Primary)
data(ERSurveyData)
NorthYears=c(2000,2001)

# Likewise, the secondary sightings are those with EXPERIMENT==2 but the LOCATION that
# is not designated as primary.  The following show that the counts match for ERSurveyData
# and the dataframe SecondarySightings that was created by ERAbund
data(SecondarySightings)

# The data in PrimarySightings are all southbound sightings for all years in which visibility and beaufort
# are less than or equal to 4. Below the counts are shown for the 2 dataframes for
# recent surveys since 1987/88.
data(PrimarySightings)
table(PrimarySightings$Start.year[PrimarySightings$Start.year>=1987])

# The following code produces the series of abundance estimates used in the paper and
# some comparative naive estimates.  It will also optionally (commented out) compute
# the variance-covariance matrix for the estimates.
#
data(PrimaryEffort)

# The data from observer MAS in 1995 is excluded because this observer
# did not participate in the double count in that year. 
PrimarySightings = PrimarySightings[!(PrimarySightings$Observer=="MAS"&PrimarySightings$Start.year==1995),]
PrimaryEffort=PrimaryEffort[!(PrimaryEffort$Observer=="MAS" & PrimaryEffort$Start.year==1995),]

#table(PrimarySightings$Start.year)

# Define arguments used in the analysis
all.years = unique(PrimaryEffort$Start.year)
recent.years = all.years[all.years>=1987]
early.years = all.years[all.years<1987]
final.time = sapply(tapply(floor(PrimaryEffort$time), 
                           PrimaryEffort$Start.year,max), 
                    function(x) ifelse(x>90,100,90))

lower.time=rep(0,length(final.time))
fn=1.0817
se.fn=0.0338

# Effort and sightings prior to 1987 were filtered for an entire watch if vis or beaufort 
# exceeded 4 at any time during the watch.  This is done for surveys starting in 1987 with the
# Use variable which is set to FALSE for all effort records in a watch if at any time the vis or
# beaufort exceeded 4 during the watch.
# Here are the hours of effort that are excluded (FALSE) and included (TRUE) by each year
# Note that for most years <1987 there are no records with Use==FALSE because the filtered records
# were excluded at the time the dataframe was constructed. The only exception is for 1978 in which  
# one watch (5 hours) was missing a beaufort value so it was excluded.
tapply(PrimaryEffort$effort,
       list(PrimaryEffort$Use, PrimaryEffort$Start.year), sum)*24

# These are the number of sightings that were included/excluded based on Use
Sightings = PrimarySightings
Sightings = merge(Sightings, subset(PrimaryEffort,select=c("key","Use")))
tapply(Sightings$Start.year, list(Sightings$Use,Sightings$Start.year), length)

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
sampled.fraction = with(Effort,
                      {
                        Day = as.numeric(as.Date(Date)-as.Date(paste(Start.year,"-12-01", sep="")))
                        tapply(effort,list(Start.year,cut(Day,seq(0,100,period))), sum)/period
                      })
# compute the number of whales counted in each period
whales.counted = with(Sightings,
                    {
                      Day = as.numeric(as.Date(Date)-as.Date(paste(Start.year,"-12-01", sep="")))
                      tapply(podsize, list(Start.year, cut(Day,seq(0, 100, period))), sum)
                    })

# Compute simple minded population estimate and plot it
period.estimate = apply(whales.counted/sampled.fraction, 1, sum, na.rm=TRUE)

as.data.frame(period.estimate) %>% 
  rownames_to_column(var = "Start.year") %>%
  transmute(Start.year = Start.year,
            Simple.estimate = period.estimate) %>%
  mutate(Season = paste0(Start.year, "/", 
                         (as.numeric(Start.year) + 1))) -> Simple.estimates

# Compute naive estimates of abundance for the 23 surveys. These use the uncorrected
# counts of whales from the primary observer during watches in which neither Beaufort nor
# vis exceeded 4.  For each year a gam with a smooth over time is fitted and this is
# used to predict total abundance throughout the migration from the counts of whales
# during the sampled periods.  There is no correction for missed pods or for measurement
# error in podsize. Each fitted migration gam is plotted with the observed values and
# saved in the file NaiveMigration.pdf.
#pdf("NaiveMigration.pdf")
naive.abundance.models=vector("list",23)
i=0
for (year in all.years){
  i=i+1
  primary=Sightings[Sightings$Start.year==year,]
  primary$Start.year=factor(primary$Start.year)
  ern=subset(Effort,
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
                                                 dformula=NULL,
                                                 plotit = FALSE)
}

Nhat.naive = sapply(naive.abundance.models,function(x) x$Total)

as.data.frame(Nhat.naive) %>% 
  rownames_to_column(var = "Start.year") %>%
  transmute(Start.year = Start.year, 
            Naive.estimate = Nhat.naive) %>%
  mutate(Season = paste0(Start.year, "/", (as.numeric(Start.year) + 1))) %>%
  select(Season, Naive.estimate) -> Naive.estimates

Reported.estimates <- read.csv(file = "Data/Nhats_2025.csv") %>%
  filter(Method == "Laake") %>%
  filter(Year0 < 2007) %>%
  transmute(Start.year = Year0,
            Season = Season,
            Reported.estimate = Nhat) %>%
  select(Season, Reported.estimate)

# Pod size corrections"
# Define set of models to be evaluated for detection
models=c("podsize+Dist+Observer",
         "podsize+Dist+Observer+beaufort",
         "podsize+Dist+Observer+vis",
         "podsize+Dist+Observer+Vis")

# Create time series of estimates based on previous approach using Reilly pod size
# correction method but using 1978 data for surveys <=1987 and 1992-1994 aerial data
# for surveys >=1992
data(add.cf.reilly)
data(add.cf.laake)
Sightings$corrected.podsize[Sightings$Start.year <= 1987] = reilly.cf(Sightings$podsize[Sightings$Start.year<=1987],
                                                                  add.cf.reilly)
Sightings$corrected.podsize[Sightings$Start.year>1987]=reilly.cf(Sightings$podsize[Sightings$Start.year>1987],
                                                                 add.cf.laake)
#pdf("ReillyApproach.pdf")
if (!file.exists("RData/Reilly_estimates.rds")){
  reilly.estimates=compute.series(models,
                                  naive.abundance.models,
                                  sightings=Sightings,
                                  effort=Effort,TruePS=FALSE)
  saveRDS(reilly.estimates, 
          file = "RData/Reilly_estimates.rds")
} else {
  reilly.estimates <- readRDS("RData/Reilly_estimates.rds")
}

Nhat.reilly=reilly.estimates$Nhat
Nhat.reilly[1:15] = Nhat.naive[1:15]*(Nhat.reilly[16]/Nhat.naive[16])
avg.Reilly.podsize=tapply(Sightings$corrected.podsize,Sightings$Start.year,mean)

as.data.frame(Nhat.reilly) %>% 
  rownames_to_column(var = "Start.year") %>%
  transmute(Start.year = Start.year,
            Reilly.estimate = Nhat.reilly) %>%
  mutate(Season = paste0(Start.year, "/", 
                         (as.numeric(Start.year) + 1))) -> Reilly.estimates

# Next compute the series of abundance estimates for the most recent 8 years by
# fitting and selecting the best detection model but not applying the pod size correction.
# From those 8 estimates and the naive estimates, compute an average ratio and 
# apply it to generate the estimates for the first 15 surveys prior to 1987.
Sightings$corrected.podsize = Sightings$podsize

if (!file.exists("RData/Nhats_nops_correct.rds")){
  abundance.estimates.nops.correction = compute.series(models, 
                                                       naive.abundance.models,
                                                       sightings=Sightings,
                                                       effort=Effort,
                                                       TruePS=FALSE)
  saveRDS(abundance.estimates.nops.correction,
          file = "RData/Nhats_nops_correct.rds")  
} else {
  abundance.estimates.nops.correction <- readRDS("RData/Nhats_nops_correct.rds")
}

ratio = abundance.estimates.nops.correction$Nhat/Nhat.naive/fn

# the following creates the same estimates as the ones in the report (Reported.estimate):
# 
if (!file.exists("RData/Nhats_Laake.rds")){
  abundance.estimates = compute.series(models,naive.abundance.models,
                                       sightings=Sightings,
                                       effort=Effort,
                                       hessian=TRUE)
  
  saveRDS(abundance.estimates,
          file = "RData/Nhats_Laake.rds")  
} else {
  abundance.estimates <- readRDS("RData/Nhats_Laake.rds")
}

Nhat.final <- abundance.estimates$Nhat

as.data.frame(Nhat.final) %>%
  rownames_to_column(var = "Start.year") %>%
  transmute(Start.year = Start.year,
            Final.estimate = Nhat.final) %>%
  mutate(Season = paste0(Start.year, "/", (as.numeric(Start.year) + 1))) -> Final.estimates

Simple.estimates %>%
  left_join(Naive.estimates, by = "Season") %>%
  left_join(Reported.estimates, by = "Season") %>%
  left_join(Reilly.estimates, by = "Season") -> all.estimates.1

ggplot(all.estimates.1) +
  geom_point(aes(x = Season, y = Simple.estimate), color = "darkgreen") +
  geom_point(aes(x = Season, y = Naive.estimate), color = "orange") +
  geom_point(aes(x = Season, y = Reported.estimate), color = "blue") +
  geom_point(aes(x = Season, y = Reilly.estimate), color = "red")

# Somehow, effort in 1987/1988 contains a lot of less than 3 hrs... All thsoe
# sightings add up to more than 1000 whales, which seems to be the cause of 
# the large differences. 

Effort %>%
  mutate(Start.year.f = as.factor(Start.year)) %>%
  group_by(Start.year.f) %>%
  summarize(Start.year = first(Start.year),
            mean.effort = mean(effort),
            min.effort = min(effort),
            max.effort = max(effort),
            total.whales = sum(nwhales)) -> effort.summary

ggplot(effort.summary) +
  geom_point(aes(x = Start.year, y = mean.effort, size = total.whales))

ggplot(effort.summary) +
  geom_point(aes(x = Start.year, y = min.effort, size = total.whales))

ggplot(Effort %>% 
         mutate(Start.year.f = as.factor(Start.year)) %>%
         group_by(Start.year.f)) +
  geom_boxplot(aes(x = Start.year.f, y = effort))

Effort %>%
  filter(effort > 1/24) %>%
  mutate(Start.year.f = as.factor(Start.year)) %>%
  group_by(Start.year.f) %>%
  summarize(Start.year = first(Start.year),
            mean.effort = mean(effort),
            min.effort = min(effort),
            max.effort = max(effort),
            total.whales = sum(nwhales)) -> effort.summary.1hr

ggplot(effort.summary.1hr) +
  geom_point(aes(x = Start.year, y = mean.effort, size = total.whales))

Effort %>%
  filter(effort < 1/24) %>%
  mutate(Start.year.f = as.factor(Start.year)) %>%
  group_by(Start.year.f) %>%
  summarize(Start.year = first(Start.year),
            mean.effort = mean(effort),
            min.effort = min(effort),
            max.effort = max(effort),
            total.whales = sum(nwhales)) -> effort.summary.less.1hr

ggplot(effort.summary.less.1hr) +
  geom_point(aes(x = Start.year, 
                 y = mean.effort, 
                 size = total.whales))


Sightings %>%
  filter(Start.year == 1987) -> Sightings.1988

Effort %>%
  filter(Start.year == 1987) -> Effort.1988

Sightings %>%
  filter(Start.year == 1985) -> Sightings.1986

Effort %>%
  filter(Start.year == 1985) -> Effort.1986

###########################################################################
###########################################################################

# My data extraction function for Laake data
# In the function below, add shift # to each effort line and create a cumulative
# effort for each shift. Then filter shifts that meet the minimum shift duration
# criteria. Currently, each line is compared to the threshold, which is not the
# right way to do it. 2025-10-01 This was completed.

min.dur <- 60
Laake.data.jags <- LaakeData2JagsInput(min.dur = min.dur, max.day = 100)

# Data in the JAGS analysis is found in the saved output file. 
ver <- 5 #out.table$model[1] #5
jm.out <- readRDS(paste0("RData/JAGS_Richards_Nmixture_v",
                         ver,
                         "a_1968to2025_min", min.dur,
                         "_NoBUGS.rds"))
# Create an observed counts dataframe
obs.day <- jm.out$jags.input$jags.data$day[,1,]
obs.n <- jm.out$jags.input$jags.data$n[,1,]
obs.prop <- jm.out$jags.input$jags.data$watch.length[,1,]
obs.prop <- rbind(c(rep(0, ncol(obs.n))),
                  obs.prop,
                  c(rep(0, ncol(obs.n))))

all.start.year <- c(jm.out$jags.input$jags.input.Laake$all.start.year,
                    jm.out$jags.input$jags.input.new$start.years)

# Create a dataframe with all years, including unsampled years.
all.years <- data.frame(start.year = seq(min(all.start.year), max(all.start.year))) %>%
  mutate(Season = paste0(start.year, "/", start.year + 1))

obs.n.df <- data.frame(start.year = rep(all.start.year, each = 100),
                       Season = rep(paste0(all.start.year, "/",
                                           all.start.year+1),
                                    each = 100),
                       Day = rep(1:100, times = length(all.start.year)),
                       n = NA,
                       prop = NA)

k1 <- k2 <- 1
for (k1 in 1:length(all.start.year)){
  obs.day.1 <- obs.day[, k1]
  obs.n.1 <- obs.n[, k1]
  obs.day.1.uniq <- na.omit(unique(obs.day.1))
  for (k2 in 1:length(obs.day.1.uniq)){
    obs.n.df[obs.n.df$start.year == all.start.year[k1] &
               obs.n.df$Day == obs.day.1.uniq[k2], "n"] <- sum(obs.n[obs.day[,k1] == obs.day.1.uniq[k2], k1], na.rm = T)
    obs.n.df[obs.n.df$start.year == all.start.year[k1] &
               obs.n.df$Day == obs.day.1.uniq[k2], "prop"] <- sum(obs.prop[obs.day[,k1] == obs.day.1.uniq[k2], k1], na.rm = T)
  }
  
}

# This is for daily estimates
N.hats.day <- data.frame(Season = rep(paste0(all.start.year, 
                                             "/", all.start.year+1),
                                      each = nrow(jm.out$jm$mean$N)), 
                         start.year = rep(all.start.year,
                                          each = nrow(jm.out$jm$mean$N)),
                         Day = rep(1:nrow(jm.out$jm$mean$N), 
                                   times = length(all.start.year)),
                         Mean = as.vector(jm.out$jm$mean$N),
                         LCL = as.vector(jm.out$jm$q2.5$N),
                         UCL = as.vector(jm.out$jm$q97.5$N),
                         mean.N = as.vector(rbind(jm.out$jm$mean$mean.N, 
                                                  rep(0, length(all.start.year)))),
                         LCL.mean.N = as.vector(rbind(jm.out$jm$q2.5$mean.N, 
                                                      rep(0, length(all.start.year)))),
                         UCL.mean.N = as.vector(rbind(jm.out$jm$q97.5$mean.N, 
                                                      rep(0, length(all.start.year)))),
                         obs.n = obs.n.df$n,
                         prop = obs.n.df$prop,
                         Method = paste0("Eguchi M", ver)) %>%
  mutate(Nhat = obs.n/prop)


obs.n.df %>%
  filter(Season == "1987/1988") -> obsd.n.1988

N.hats.day %>%
  filter(Season == "1987/1988") -> N.hats.day.1988

ggplot(N.hats.day.1988) +
  geom_line(aes(x = Day, y = mean.N)) +
  geom_ribbon(aes(x = Day, ymin = LCL.mean.N, ymax = UCL.mean.N),
              fill = "gold", alpha = 0.5) +
  geom_point(aes(x = Day, y = Nhat)) 

effort.jags.1988 <- Laake.data.jags$jags.data$watch.prop[,1,16]
