# Checking input data for 1990 to 2006 (start years) because of
# wildly different Nhats between Laake's and mine. Other years
# seem to be in the same ballpark. 
# 

rm(list=ls())
library(tidyverse)

source("Granite_Canyon_Counts_fcns.R")

min.dur <- 30

## get Laake's data:
load("Data/Primary.rda")      # on-effort sightings
load("Data/ERSurveyData.rda")
load("Data/Observer.rda")

# The data in PrimarySightings are all southbound sightings for all years in which visibility and beaufort
# are less than or equal to 4. Below the counts are shown for the 2 dataframes for
# recent surveys since 1987/88.
#table(Primary$Start.year[Primary$vis<=4 & Primary$beaufort<=4])
#load("Data/PrimarySightings.rda")
load("Data/PrimaryEffort.rda")
load("Data/SecondaryEffort.rda")
#load("Data/SecondarySightings.rda")

# Filter effort and sightings and store in dataframes Effort and Sightings
Laake_PrimaryEffort <- PrimaryEffort[PrimaryEffort$Use,]  

Laake_SecondaryEffort <- SecondaryEffort[SecondaryEffort$Use,]

# Filter to minimum duration:
Laake_PrimaryEffort %>%
  filter(effort >= min.dur / (60 * 24) ) -> Laake_PrimaryEffort

Laake_SecondaryEffort %>%
  filter(effort >= min.dur / (60 * 24) ) -> Laake_SecondaryEffort

Laake_PrimaryEffort %>%
  filter(Start.year > 1989, 
         Start.year < 2007) -> Primary.data.1990to2006

Laake_SecondaryEffort %>%
  filter(Start.year > 1989, 
         Start.year < 2007) -> Secondary.data.1990to2006

all.years <- unique(Laake_PrimaryEffort$Start.year)
years.1990to2006 <- unique(Primary.data.1990to2006$Start.year)

Primary.data.1990to2006 %>%
  select(Start.year, time, effort, nwhales) %>%
  rename(n = nwhales) %>%
  mutate(method = "Laake") -> primary.Laake.1990to2006

Secondary.data.1990to2006 %>%
  select(Start.year, time, effort, nwhales) %>%
  rename(n = nwhales) %>%
  mutate(method = "Laake") -> secondary.Laake.1990to2006

## Get Jags input using the function:
jags.input <- LaakeData2JagsInput(min.dur)

# Extract jags$n for the same years (1990 to 2006)
jags.n.1990to2006 <- jags.input$n[,,all.years %in% years.1990to2006]
jags.effort.1990to2006 <- jags.input$watch.prop[,,all.years %in% years.1990to2006]
jags.day.1990to2006 <- jags.input$day[,,all.years %in% years.1990to2006]

primary.jags.1990to2006 <- data.frame(Start.year = rep(years.1990to2006, 
                                                       each = dim(jags.n.1990to2006)[1]),
                                      time = as.vector(jags.day.1990to2006[,1,]),
                                      effort = as.vector(jags.effort.1990to2006[,1,]),
                                      n = as.vector(jags.n.1990to2006[,1,]),
                                      method = "jags") 



secondary.jags.1990to2006 <- data.frame(Start.year = rep(years.1990to2006, 
                                                         each = dim(jags.n.1990to2006)[1]),
                                        time = as.vector(jags.day.1990to2006[,2,]),
                                        effort = as.vector(jags.effort.1990to2006[,2,]),
                                        n = as.vector(jags.n.1990to2006[,2,]),
                                        method = "jags") 

all.primary.1990to2006 <- rbind(primary.Laake.1990to2006, 
                                primary.jags.1990to2006)

all.secondary.1990to2006 <- rbind(secondary.Laake.1990to2006, 
                                  secondary.jags.1990to2006)

ggplot(all.primary.1990to2006) +
  geom_point(aes(x = time, y = n, color = method)) +
  facet_wrap(~ Start.year)

ggplot(all.secondary.1990to2006) +
  geom_point(aes(x = time, y = n, color = method)) +
  facet_wrap(~ Start.year)

# They look fine... 

ggplot(all.primary.1990to2006) +
  geom_point(aes(x = time, y = effort, color = method)) +
  facet_wrap(~ Start.year)

ggplot(all.secondary.1990to2006) +
  geom_point(aes(x = time, y = effort, color = method)) +
  facet_wrap(~ Start.year)

