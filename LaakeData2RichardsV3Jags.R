# LaakeData2RichardsV3Jags
# 
# Creates JAGS input data from Laake's ERAnalysis library and runs JAGS on 
# model_Richards_pois_bino_v3_Laake.txt. 
# 
# This v2 version has year spcific K parameters. It has conversion issues but so was
# v1 with a common K parameter.

# This Laake version allows different numbers of sampling periods between primary
# and secondary sampling stations.
# 

# Code chunks with var1 = expressions are by Laake. I use var1 <- expressions

rm(list = ls())

library(ERAnalysis)
library(tidyverse)
library(ggplot2)
library(loo)
library(bayesplot)

source("Granite_Canyon_Counts_fcns.R")

jags.model <- "models/model_Richards_pois_bino_v3_Laake.txt"

# Estimates from Laake et al. are here:
col.defs <- cols(Year = col_character(),
                 Nhat = col_double(),
                 CV = col_double())

Laake.estimates <- read_csv(file = "Data/Laake et al 2012 Table 9 Nhats.csv",
                            col_types = col.defs) %>% 
  mutate(SE = CV * Nhat,
         LCL = Nhat - 1.96 * SE,
         UCL = Nhat + 1.96 * SE,
         Season = lapply(strsplit(Year, "_"), 
                         FUN = function(x) paste0(x[1], "/", x[2])) %>% 
           unlist) %>%
  dplyr::select(Season, Nhat, SE, LCL, UCL)  %>%
  mutate(Year = lapply(str_split(Season, "/"), 
                       FUN = function(x) x[2]) %>% 
           unlist() %>% 
           as.numeric())

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
#data(Primary)      # on-effort sightings
data(ERSurveyData)
data("Observer")

# The data in PrimarySightings are all southbound sightings for all years in which visibility and beaufort
# are less than or equal to 4. Below the counts are shown for the 2 dataframes for
# recent surveys since 1987/88.
#table(Primary$Start.year[Primary$vis<=4 & Primary$beaufort<=4])
#data(PrimarySightings)
#data(PrimaryEffort)

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
# Effort = PrimaryEffort[PrimaryEffort$Use,]  
# 
# Sightings = PrimarySightings
# Sightings$seq = 1:nrow(Sightings)
# Sightings = merge(Sightings, subset(Effort, select=c("key")))
# Sightings = Sightings[order(Sightings$seq),]

# sightings
Laake_PrimarySightings <- read.csv(file = "Data/Laake_PrimarySightings.csv")
Laake_SecondarySightings <- read.csv(file = "Data/Laake_SecondarySightings.csv")

# effort
Laake_PrimaryEffort <- read.csv(file = "Data/Laake_PrimaryEffort.csv") %>%
  filter(Use)
Laake_SecondaryEffort <- read.csv(file = "Data/Laake_SecondaryEffort.csv") %>%
  filter(Use)

# observers
obs <- c(unique(Laake_PrimaryEffort$Observer),
         unique(Laake_SecondaryEffort$Observer)) %>% unique()
Laake.obs <- data.frame(obs = obs,
                        obs.ID = seq(1:length(obs)))

Laake_PrimaryEffort %>% 
  group_by(Date) %>% 
  summarize(sum.effort = sum(effort)) -> tmp
# maximum daily effort was 0.447917, or 10.75 hrs.

max.effort <- max(tmp$sum.effort)
Laake_PrimaryEffort %>% 
  group_by(Start.year) %>%
  reframe(Date = Date,
            Year = first(Start.year),
            n = nwhales,
            Day = as.Date(Date) - as.Date(paste0(Year, "-12-01")),
            Season = paste0(Year, "/", Year + 1),
            obs = Observer,
            vs = vis,
            bf = beaufort,
            watch.prop = effort/max.effort) -> Laake.primary.counts

# Although the maximum effort for the secondary effort was 10 hrs, I use
# the same max effort as the primary
# Laake_SecondaryEffort %>% 
#   group_by(Date) %>% 
#   summarize(sum.effort = sum(effort)) -> tmp

Laake_SecondaryEffort %>% 
  group_by(Start.year) %>%
  reframe(Date = Date,
            Year = first(Start.year),
            n = nwhales,
            Day = as.Date(Date) - as.Date(paste0(Year, "-12-01")),
            Season = paste0(Year, "/", Year + 1),
            obs = Observer,
            vs = vis,
            bf = beaufort,
            watch.prop = effort/max.effort) -> Laake.secondary.counts

Laake.primary.counts %>% 
  group_by(Date) %>%
  reframe(Year = first(Start.year),
            Date = Date,
            Day = as.Date(Date) - as.Date(paste0(Year, "-12-01")),
            Season = paste0(Year, "/", Year + 1),
            Daily.n = sum(n)) -> Laake.primary.daily.counts


# Richards function fit starts here:
run.date.Laake <- Sys.Date() #"2023-09-29" # "2023-08-11"
out.file.name <- paste0("RData/JAGS_Richards_v3_Laake_", 
                        run.date.Laake, ".rds")

# Find double observer years
double.obs.year <- unique(Laake_SecondaryEffort$Start.year) 
all.year <- unique(Laake_PrimaryEffort$Start.year)

n.station <- rep(1, length(all.year))
n.station[all.year %in% double.obs.year] <- 2

Laake.primary.counts %>%
  group_by(Start.year) %>%
  summarize(Season = first(Season),
            periods = first(n())) -> Laake.primary.periods

Laake.secondary.counts %>%
  group_by(Start.year) %>%
  summarize(Season = first(Season),
            periods = first(n())) -> Laake.secondary.periods

Laake.primary.periods %>% 
  left_join(Laake.secondary.periods, by = "Season") %>%
  select(Start.year.x, Season, periods.x, periods.y) %>%
  transmute(Start.year = Start.year.x,
            Season = Season,
            periods.1 = periods.x,
            periods.2 = periods.y) -> Laake.periods

day <- bf <- vs <- watch.prop <- obs.input <- n.Laake <- array(dim = c(max(Laake.primary.periods$periods),
                                                                       2, length(all.year)))

y <- 3
y2 <- 1
for (y in 1:length(all.year)){
  temp.data <- Laake.primary.counts %>%
    filter(Start.year == all.year[y]) %>%
    arrange(Day) %>%
    left_join(Laake.obs, by = "obs")
  
  n.Laake[1:Laake.primary.periods$periods[y], 1, y] <- temp.data$n
  obs.input[1:Laake.primary.periods$periods[y], 1, y] <- temp.data$obs.ID
  watch.prop[1:Laake.primary.periods$periods[y], 1, y] <- temp.data$watch.prop
  bf[1:Laake.primary.periods$periods[y], 1, y] <- scale(temp.data$bf)
  vs[1:Laake.primary.periods$periods[y], 1, y] <- scale(temp.data$vs)
  day[1:Laake.primary.periods$periods[y], 1, y] <- as.numeric(temp.data$Day)
  
  # fill in the secondary observations
  if (isTRUE(all.year[y] %in% double.obs.year)){
    temp.data <- Laake.secondary.counts %>%
      filter(Start.year == all.year[y])%>%
      left_join(Laake.obs, by = "obs")
    
    n.Laake[1:Laake.secondary.periods$periods[y2], 2, y] <- temp.data$n
    obs.input[1:Laake.secondary.periods$periods[y2], 2, y] <- temp.data$obs.ID
    watch.prop[1:Laake.secondary.periods$periods[y2], 2, y] <- temp.data$watch.prop
    bf[1:Laake.secondary.periods$periods[y2], 2, y] <- scale(temp.data$bf)
    vs[1:Laake.secondary.periods$periods[y2], 2, y] <- scale(temp.data$vs)
    day[1:Laake.secondary.periods$periods[y2], 2, y] <- as.numeric(temp.data$Day)
    
    y2 <- y2 + 1
    
  }
}

jags.data.Laake <- list(  n = n.Laake,
                          n.station = n.station,
                          n.year = length(all.year),
                          n.obs = nrow(Laake.obs),
                          periods = Laake.periods %>% 
                            select(periods.1, periods.2) %>% simplify2array(),
                          obs = obs.input,
                          vs = vs,
                          bf = bf,
                          watch.prop = watch.prop,
                          day = day)

jags.params <- c("OBS.RF", "OBS.Switch",
                 "BF.Switch", "BF.Fixed",
                 "VS.Switch", "VS.Fixed",
                 "mean.prob", "mean.N", "Max",
                 "Corrected.Est", "Raw.Est", "N",
                 "K", "S1", "S2", "P",
                 "Max.alpha", "Max.beta",
                 "S1.alpha", "S2.alpha",
                 "S1.beta", "S2.beta",
                 "P.alpha", "P.beta",
                 "log.lkhd")

MCMC.params <- list(n.samples = 250000,
                    n.thin = 100,
                    n.burnin = 200000,
                    n.chains = 5)

if (!file.exists(out.file.name)){
  
  Start_Time<-Sys.time()
  
  jm <- jagsUI::jags(jags.data.Laake,
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
  jm.out.Laake <- list(jm = jm,
                       jags.data = jags.data.Laake,
                       jags.params = jags.params,
                       jags.model = jags.model,
                       MCMC.params = MCMC.params,
                       Run_Time = Run_Time,
                       Sys.env = Sys.getenv())
  
  saveRDS(jm.out.Laake,
          file = out.file.name)
  
} else {
  
  jm.out.Laake <- readRDS(out.file.name)
}

# need to turn zeros into NAs when there were no second station:
data.array <- jm.out.Laake$jags.data$n
data.array[,2,which(jm.out.Laake$jags.data$n.station == 1)] <- NA

LOOIC.n.Laake <- compute.LOOIC(loglik.array = jm.out.Laake$jm$sims.list$log.lkhd,
                               data.array = data.array,
                               MCMC.params = jm.out.Laake$MCMC.params)

# Look at Rhat statistics
max.Rhat.Laake <- lapply(jm.out.Laake$jm$Rhat, FUN = max, na.rm = T) %>%
  unlist()
max.Rhat.Laake.big <- max.Rhat.Laake[which(max.Rhat.Laake > 1.1)]



Nhat.Laake.df <- data.frame(median = jm.out.Laake$jm$q50$Corrected.Est,
                      LCL = jm.out.Laake$jm$q2.5$Corrected.Est,
                      UCL = jm.out.Laake$jm$q97.5$Corrected.Est,
                      mean = jm.out.Laake$jm$mean$Corrected.Est,
                      Season =  Laake.estimates$Season)

Nhat.Laake.df %>% 
  left_join(Laake.estimates, by = "Season") %>%
  mutate(delta_Nhat = median - Nhat) -> all.estimates

ggplot(all.estimates) +
  geom_point(aes(x = Season, y = median ),
             color = "blue") +
  #geom_errorbar(aes(x = Season, ymin = LCL.x, ymax = UCL.x),
  #              color = "blue") +
  geom_point(aes(x = Season, y = Nhat ),
             color = "green") +
  geom_errorbar(aes(x = Season, ymin = LCL.y, ymax = UCL.y),
                color = "green")

# mcmc_dens(jm.out.Laake$jm$samples, c("Max.alpha", "Max.beta", 
#                                      "S1.alpha", "S1.beta", 
#                                      "S2.alpha", "S2.beta",
#                                      "P.alpha", "P.beta"))

P1 <- c("P[1]", "P[2]", "P[3]", "P[4]", "P[5]", "P[6]", "P[7]", "P[8]", "P[9]")
P2 <- c("P[10]", "P[11]", "P[12]", "P[13]", "P[14]", "P[15]", "P[16]", "P[17]", "P[18]")
P3 <- c("P[19]", "P[20]", "P[21]", "P[22]", "P[23]")

K3 <- c("K[19]", "K[20]", "K[21]", "K[22]", "K[23]")
K2 <- c("K[10]", "K[11]", "K[12]", "K[13]", "K[14]", "K[15]", "K[16]", "K[17]", "K[18]")
K1 <- c("K[1]", "K[2]", "K[3]", "K[4]", "K[5]", "K[6]", "K[7]", "K[8]", "K[9]")

S1.1 <- c("S1[1]", "S1[2]", "S1[3]", "S1[4]", "S1[5]", "S1[6]", "S1[7]", "S1[8]", "S1[9]")
S2.1 <- c("S2[1]", "S2[2]", "S2[3]", "S2[4]", "S2[5]", "S2[6]", "S2[7]", "S2[8]", "S2[9]")

S1.2 <- c("S1[10]", "S1[11]", "S1[12]", "S1[13]", "S1[14]", "S1[15]", "S1[16]", "S1[17]", "S1[18]")
S2.2 <- c("S2[10]", "S2[11]", "S2[12]", "S2[13]", "S2[14]", "S2[15]", "S2[16]", "S2[17]", "S2[18]")

S1.3 <- c("S1[19]", "S1[20]", "S1[21]", "S1[22]", "S1[23]")
S2.3 <- c("S2[19]", "S2[20]", "S2[21]", "S2[22]", "S2[23]")

