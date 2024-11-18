# LaakeData2Jags
# 
# Creates JAGS input data from Laake's ERAnalysis library and runs JAGS.
# 
# Code chunks with var1 = expressions are by Laake. I use var1 <- expressions

#rm(list = ls())

#library(ERAnalysis)
# because I can't build new ERAnalysis library, I moved all the ERAnalysis data 
# files to the Granite_Canyon_Counts/Data folder. This way, I can still use the
# data files using the data() function - I can't use it on my linux laptop... So,
# I changed it to the load function. 2024-11-14

# min.dur is the minimum effort duration in minutes to be included in the analysis
LaakeData2JagsInput <- function(min.dur){
  
  library(tidyverse)
  # library(ggplot2)
  # library(R2WinBUGS)
  
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
  load("Data/Primary.rda")      # on-effort sightings
  load("Data/ERSurveyData.rda")
  load("Data/Observer.rda")
  
  # The data in PrimarySightings are all southbound sightings for all years in which visibility and beaufort
  # are less than or equal to 4. Below the counts are shown for the 2 dataframes for
  # recent surveys since 1987/88.
  #table(Primary$Start.year[Primary$vis<=4 & Primary$beaufort<=4])
  load("Data/PrimarySightings.rda")
  load("Data/PrimaryEffort.rda")
  load("Data/SecondaryEffort.rda")
  
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
  Laake_PrimaryEffort = PrimaryEffort[PrimaryEffort$Use,]  
  
  Laake_SecondaryEffort <- SecondaryEffort[SecondaryEffort$Use,]

    # Filter to minimum duration:
  Laake_PrimaryEffort %>%
    filter(effort >= min.dur / (60 * 24) ) -> Laake_PrimaryEffort
  
  Laake_SecondaryEffort %>%
    filter(effort >= min.dur / (60 * 24) ) -> Laake_SecondaryEffort
  
  Sightings = PrimarySightings
  Sightings$seq = 1:nrow(Sightings)
  Sightings = merge(Sightings, subset(Laake_PrimaryEffort, select=c("key")))
  Sightings = Sightings[order(Sightings$seq),]
  
  # filter off-effort sightings and high Beaufort/vis lines from secondary sightings
  # but... there is no effort data for the secondary sightings... so, can't use it
  # for BUGS/jags - ignore it for now.
  # SecondarySightings %>% 
  #   filter(vis < 5, beaufort < 5, is.na(off)) -> secondary.sightings
  
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
  
  # Summarizing by day worked fine but the model requires counts per observation
  # period. Needs to be redone. 2023-09-15 DONE.
  
  Laake_PrimaryEffort %>% 
    mutate(Day1 = as.Date(paste0(Start.year, "-12-01")),
           dt = as.numeric(as.Date(Date) - Day1) + 1,
           obs = Observer) %>%
    select(Start.year, nwhales, effort, vis, beaufort, obs, dt) %>%
    group_by(Start.year) %>%
    mutate(effort.min = effort * 24 * 60,
           watch.prop = effort.min/540) -> Effort.by.period
  
  # observers
  obs <- c(unique(Laake_PrimaryEffort$Observer),
           unique(Laake_SecondaryEffort$Observer)) %>% unique()
  Laake.obs <- data.frame(obs = obs,
                          obs.ID = seq(1:length(obs)))
  
  # Laake_PrimaryEffort %>% 
  #   group_by(Date) %>% 
  #   summarize(sum.effort = sum(effort)) -> tmp
  # maximum daily effort was 0.447917 days, or 10.75 hrs.
  # Effort is measured in decimal days
  
  #max.effort <- max(tmp$sum.effort)
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
            watch.prop = effort) -> Laake.primary.counts
  
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
            watch.prop = effort) -> Laake.secondary.counts
  
  Laake.primary.counts %>% 
    group_by(Date) %>%
    reframe(Year = first(Start.year),
            Date = Date,
            Day = as.Date(Date) - as.Date(paste0(Year, "-12-01")),
            Season = paste0(Year, "/", Year + 1),
            Daily.n = sum(n)) -> Laake.primary.daily.counts
  
  
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
    bf[1:Laake.primary.periods$periods[y], 1, y] <- temp.data$bf #scale(temp.data$bf)
    vs[1:Laake.primary.periods$periods[y], 1, y] <- temp.data$vs #scale(temp.data$vs)
    day[1:Laake.primary.periods$periods[y], 1, y] <- as.numeric(temp.data$Day)
    
    # fill in the secondary observations
    if (isTRUE(all.year[y] %in% double.obs.year)){
      temp.data <- Laake.secondary.counts %>%
        filter(Start.year == all.year[y])%>%
        left_join(Laake.obs, by = "obs")
      
      n.Laake[1:Laake.secondary.periods$periods[y2], 2, y] <- temp.data$n
      obs.input[1:Laake.secondary.periods$periods[y2], 2, y] <- temp.data$obs.ID
      watch.prop[1:Laake.secondary.periods$periods[y2], 2, y] <- temp.data$watch.prop
      bf[1:Laake.secondary.periods$periods[y2], 2, y] <- temp.data$bf #scale(temp.data$bf)
      vs[1:Laake.secondary.periods$periods[y2], 2, y] <- temp.data$vs #scale(temp.data$vs)
      day[1:Laake.secondary.periods$periods[y2], 2, y] <- as.numeric(temp.data$Day)
      
      y2 <- y2 + 1
      
    }
  }
  
  jags.data <- list(  n = n.Laake,
                      n.station = n.station,
                      n.year = length(all.year),
                      n.obs = nrow(Laake.obs),
                      periods = Laake.periods %>% 
                        select(periods.1, periods.2) %>% simplify2array(),
                      obs = obs.input,
                      vs = vs,
                      bf = bf,
                      watch.prop = watch.prop,
                      day = day,
                      n.days = 94)
  
  return(jags.data)  
}



