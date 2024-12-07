# plot_Jags_results.R
# Creates plots from results of Jags analysis
# It uses output from either Jags_Richards_AllData.R or Jags_Richards_NoLaakeData.R.
# 

rm(list=ls())

source("Granite_Canyon_Counts_fcns.R")
.data <- "all" #"no Laake" #  or 
min.dur <- 30 #85 or 30 #
WinBUGS.years <- c(2007, 2008, 2010, 2011, 2015, 2016, 
                   2020, 2022, 2023, 2024)

jags.only.jm.out <- readRDS(paste0("RData/JAGS_Richards_pois_bino_v5_min",
                                   min.dur, "_since2010_2024-12-05.rds"))

jags.NoBUGS.jm.out <- readRDS(paste0("RData/JAGS_Richards_pois_bino_v5_min",
                                     min.dur, "_NoBUGS_2024-12-06.rds"))

jags.NoBUGS.all.years <- c(jags.NoBUGS.jm.out$jags.input$jags.input.Laake$all.start.year,
                           jags.NoBUGS.jm.out$jags.input$jags.input.new$start.years)

jags.NoBUGS.all.seasons <- paste0(jags.NoBUGS.all.years, "/",
                                  jags.NoBUGS.all.years + 1)

if (min.dur == 85){
  jags.run.date <- "2024-12-04"  # for min.dur == 85 and "no Laake"
  
} else {
  jags.run.date <- "2024-12-03"  # for min.dur == 30
  
}

BUGS.run.date <- "2024-11-23"
if (.data == "all"){

  jm.out <- readRDS(paste0("RData/JAGS_Richards_pois_bino_v5_min", min.dur, "_AllYears_",
                      jags.run.date, ".rds"))
} else {
  jags.run.date <- "2024-12-04"
  jm.out <- readRDS(paste0("RData/JAGS_Richards_pois_bino_v5_min", 
                    min.dur, "_Since2006_",
                    jags.run.date, ".rds"))
}

jags.input <- jm.out$jags.input

# make all start years in numeric
start.years <- as.numeric(jags.input$start.years)

# Create a dataframe with all years, including unsampled years.
all.years <- data.frame(start.year = seq(min(start.years), 
                                   max(start.years))) %>%
  mutate(Season = paste0(start.year, "/", start.year + 1))

# Look at the annual abundance estimates:
Nhat. <- data.frame(Season = paste0(start.years, "/", 
                                    start.years+1),
                    Nhat = jm.out$jm$q50$Corrected.Est,
                    LCL = jm.out$jm$q2.5$Corrected.Est,
                    UCL = jm.out$jm$q97.5$Corrected.Est) %>%
  right_join(all.years, by = "Season") %>%
  arrange(start.year) %>%
  mutate(Method = "Richards")

Nhat.NoBUGS <- data.frame(Season = jags.NoBUGS.all.seasons,
                          Nhat = jags.NoBUGS.jm.out$jm$q50$Corrected.Est,
                          LCL = jags.NoBUGS.jm.out$jm$q2.5$Corrected.Est,
                          UCL = jags.NoBUGS.jm.out$jm$q97.5$Corrected.Est) %>%
  right_join(all.years, by = "Season") %>%
  arrange(start.year) %>%
  mutate(Method = "Richards-NoBUGS")


Nhat.jags.only <- data.frame(Season = jags.only.jm.out$jags.input$seasons,
                             Nhat = jags.only.jm.out$jm$q50$Corrected.Est,
                             LCL = jags.only.jm.out$jm$q2.5$Corrected.Est,
                             UCL = jags.only.jm.out$jm$q97.5$Corrected.Est) %>%
  right_join(all.years, by = "Season") %>%
  arrange(start.year) %>%
  mutate(Method = "Richards-JagsOnly")

#N.hats.day.jags.only <- 
# This is for daily estimates
N.hats.day <- data.frame(Season = rep(paste0(start.years, "/", start.years+1), 
                                      each = nrow(jm.out$jm$mean$N)), #rep(Nhat.$Season, each = nrow(jm.out$jm$mean$N)),
                         Day = rep(1:nrow(jm.out$jm$mean$N), times = length(start.years)),
                         Mean = as.vector(jm.out$jm$mean$N),
                         LCL = as.vector(jm.out$jm$q2.5$N),
                         UCL = as.vector(jm.out$jm$q97.5$N)) 

# Daily estimates plots
p.daily.Richards <- ggplot(N.hats.day %>% group_by(Season)) + 
  geom_ribbon(aes(x = Day, ymin = LCL, ymax = UCL),
              fill = "blue", alpha = 0.5) +
  geom_point(aes(x = Day, y = Mean), 
             shape = 4,
             color =  "green") +
  geom_path(aes(x = Day, y = Mean)) + 
  facet_wrap(~ Season)

# These are not the best estimates
Reported.estimates <- read.csv(file = "Data/all_estimates_2024.csv") %>%
  transmute(Season = Season,
            Nhat = Nhat,
            LCL = LCL,
            UCL = UCL,
            Method = paste0(Method, "-Reported")) %>%
  right_join(all.years, by = "Season") %>%
  arrange(start.year) %>%
  relocate(Method, .after = start.year)

WinBUGS.out <- readRDS(file = paste0("RData/WinBUGS_", 
                                     min(WinBUGS.years), "to", 
                                     max(WinBUGS.years), "_v2_min", 
                                     min.dur, "_",
                                     BUGS.run.date, ".rds"))

Corrected.Est <- WinBUGS.out$BUGS.out$sims.list$Corrected.Est
WinBUGS.seasons <- WinBUGS.out$BUGS.input$seasons 

WinBUGS.abundance.df <- data.frame(Season = WinBUGS.seasons,
                                  Nhat = apply(Corrected.Est,
                                               FUN = mean,
                                               MARGIN = 2),
                                  LCL = apply(Corrected.Est, 
                                              MARGIN = 2, 
                                              FUN = quantile, 0.025),
                                  UCL = apply(Corrected.Est, 
                                              MARGIN = 2, 
                                              FUN = quantile, 0.975)) %>%
  right_join(all.years, by = "Season") %>%
  arrange(start.year) %>%
  mutate(Method = "Durban")

# Create a dataframe for daily estimates:
daily.estim <- WinBUGS.out$BUGS.out$sims.list$Daily.Est

# get stats:
mean.mat <- LCL.mat <- UCL.mat <- matrix(data = NA, 
                                         nrow = dim(daily.estim)[2], 
                                         ncol = dim(daily.estim)[3])

for (k1 in 1:dim(daily.estim)[2]){
  for (k2 in 1:dim(daily.estim)[3]){
    mean.mat[k1, k2] <- mean(daily.estim[,k1,k2])
    LCL.mat[k1, k2] <- quantile(daily.estim[,k1,k2], 0.025)
    UCL.mat[k1, k2] <- quantile(daily.estim[,k1,k2], 0.975)
  }
  
}

N.hats.day.Durban <- data.frame(Season = rep(WinBUGS.seasons,
                                             each = dim(daily.estim)[2]),
                                Day = rep(1:dim(daily.estim)[2],
                                          length(WinBUGS.seasons)),
                                Mean = as.vector(mean.mat),
                                LCL = as.vector(LCL.mat),
                                UCL = as.vector(UCL.mat))

# Daily estimates plots
p.daily.Durban <- ggplot(N.hats.day.Durban %>% group_by(Season)) + 
  geom_ribbon(aes(x = Day, ymin = LCL, ymax = UCL),
              fill = "blue", alpha = 0.5) +
  geom_path(aes(x = Day, y = Mean)) + 
  #geom_point(data = Nhats.HT.all[Nhats.HT.all$Season %in% seasons,],
  #           aes(x = Day, y = Nhat),
  #           alpha = 0.3) + 
  facet_wrap(~ Season)

# Include non-survey years - no Laake estimates for 2007/2008 because I don't have
# raw data for that year. Only the WinBUGS inputs. 
Laake.abundance.new <- read.csv(file = "Data/all_estimates_Laake_2024.csv") %>%
  mutate(LCL = CL.low,
         UCL = CL.high) %>%
  select(c(Season, Nhat, LCL, UCL)) %>%
  right_join(all.years, by = "Season") %>%
  arrange(start.year) %>%
  mutate(Method = "Laake")

# In reported estimates, there are two 2006/2007.
Reported.estimates %>%
  na.omit() %>%
  select(Season) %>% 
  unique() -> sampled.seasons 

Laake.abundance.new %>%
  rbind(WinBUGS.abundance.df) %>%
  rbind(Nhat.) %>%
  rbind(Nhat.jags.only) %>%
  rbind(Nhat.NoBUGS) -> all.estimates

p.Nhats <- ggplot(all.estimates) +
  geom_point(aes(x = start.year, y = Nhat,
                 color = Method),
             alpha = 0.5) +
  geom_errorbar(aes(x = start.year, 
                    ymin = LCL, ymax = UCL,
                    color = Method)) +
  ylim(0, max(all.estimates$UCL, na.rm = T)+1000) + 
  labs(title = paste0("Min. dur. = ", min.dur, " minutes"))

#### Some observations in results over time. ###########################
# Precision of the estimates using Richards function is a lot better than other 
# methods. They are a bit lower (negatively biased) than other two. I think the
# problem is detection probabilities. mean.prob gets very small (mean = 0.04) 
# even though the prior is beta(0.95, 0.05). 

# Changed the prior to unif(0.9, 1.0) - 2024-10-21

# This problem has been fixed as of November 1, 2024. The problem was how I set
# up the proportion of watch period. It was set at hours observed over maximum
# possible. But, this should have been over 24 hrs (or one day), which was 
# already computed when data were extracted. By dividing that (times 60 min times
# 24 hrs) by 540 minutes made those numbers bigger by 2.667 times (24*60/540), 
# resulting in a smaller estimated abundance. 

# Although the estimates are a lot better now they are a bit off, especially in 
# some years. These need to be looked at carefully. 2024-11-13

####    ###########################

