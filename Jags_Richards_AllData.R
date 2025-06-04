# Jags_Richards_AllData.R
# 
# Combines Laake data and more recent data and runs
# model_Richards_pois_bino_vX, where X is version number. See below.  
# Some diagnostics are conducted.

# M1 <- (1 + (2 * exp(K) - 1) * exp((1/S1) * (P - d))) ^ (-1/exp(K))
# M2 <- (1 + (2 * exp(K) - 1) * exp((1/S2) * (P - d))) ^ (-1/exp(K))
# N <- min + (max - min) * (M1 * M2)
#
# d is the number of days from the beginning of nesting season
# S1 < 0 and S2 > 0 define the "fatness" of the function
# K > 0 defines the "flatness" at the peak of the function
# P defines where the peak is relative to the range of d; min(d) < P < max(d)
# min is "the basal level of the number of whales outside the migration season"
# max > min. min == 0.

# Model version - depending on which parameters are constant/year-specific
# v3: Max, P, S1, S2, and K are time specific
# v4: Max, S1, and S2 are time specific. K and P are constant.
# v5: Max, P, S1 and S2 are time specific. K is constant.

rm(list = ls())

#library(ERAnalysis)
library(tidyverse)
library(ggplot2)
library(loo)
library(bayesplot)

source("Granite_Canyon_Counts_fcns.R")

# "2025-04-21" for v3/60 
# "2025-04-17" for v4/60
# "2025-04-18" for v5/60

Run.date <- Sys.Date() #"2025-04-21" #"2025-04-17" #

# Minimum length of observation periods in minutes
# In order to run a new minimum duration, WinBUGS needs to be run first.
min.dur <- 60 #10 #85 #

ver <- "v3"

# These are the ending year of each season - different from the standard, which
# uses the beginning year. So, for example, 2022 in the following vector indicates
# for the 2021/2022 season
years <- c(2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024, 2025)

# MCMC.params <- list(n.samples = 250000,
#                     n.thin = 100,
#                     n.burnin = 200000,
#                     n.chains = 5)

# v3 does not converge well with the above MCMC setting so increasing samples
MCMC.params <- list(n.samples = 750000,
                    n.thin = 500,
                    n.burnin = 500000,
                    n.chains = 5)

# MCMC.params <- list(n.samples = 2500,
#                     n.thin = 10,
#                     n.burnin = 200,
#                     n.chains = 5)

jags.params <- c("OBS.RF", "BF.Fixed",
                 "VS.Fixed",
                 "mean.prob", "mean.N", "Max",
                 "Corrected.Est", "Raw.Est", "N",
                 "K", "S1", "S2", "P",
                 "Max.alpha", "Max.beta",
                 "S1.alpha", "S2.alpha",
                 "S1.beta", "S2.beta",
                 "P.alpha", "P.beta",
                 "K.alpha", "K.beta",
                 "N.alpha",
                 "log.lkhd")

###############################
# The following function uses "new" data since 2010 as well as those from Laake's 
# analysis to compute abundance since the 1967/1968 season. 
jm.out <- NoBUGS_Richards_fcn(min.dur = min.dur, 
                              ver = ver, 
                              years = years, 
                              data.dir = "RData/V2.1_Feb2025", 
                              jags.params = jags.params, 
                              MCMC.params = MCMC.params,
                              Run.date = Run.date,
                              max.day = 100)

# need to turn zeros into NAs when there were no second station:
data.array <- jm.out$jags.input$jags.data$n
data.array[,2,which(jm.out$jags.input$jags.data$n.station == 1)] <- NA

LOOIC.n <- compute.LOOIC(loglik.array = jm.out$jm$sims.list$log.lkhd,
                         data.array = data.array,
                         MCMC.params = jm.out$MCMC.params)

# There are some (< 0.5%) bad ones. I should look at which ones are not fitting well.

# Look at Rhat statistics
max.Rhat <- lapply(jm.out$jm$Rhat, FUN = max, na.rm = T) %>%
  unlist()
max.Rhat.big <- max.Rhat[which(max.Rhat > 1.1)]

mcmc_dens(jm.out$jm$samples, c("S1.alpha", "S1.beta",
                               "S2.alpha", "S2.beta",
                               "P.alpha", "P.beta",
                               "K.alpha", "K.beta"))
# P.alpha and P.beta seem to be not behaving well - the right tails are not 
# captured. 
mcmc_trace(jm.out$jm$samples, c("S1.alpha", "S1.beta",
                                "S2.alpha", "S2.beta",
                                "P.alpha", "P.beta",
                                "K.alpha", "K.beta"))

mcmc_dens(jm.out$jm$samples, c("BF.Fixed", "VS.Fixed"))

# v4 has one P and one K.
par.idx <- c(1:nrow(jm.out$jm$mean$S1))

if (ver == "v4"){
  mcmc_trace(jm.out$jm$samples, c("P", "K"))
  mcmc_dens(jm.out$jm$samples, c("P", "K"))
} else if (ver == "v3"){
  mcmc_trace(jm.out$jm$samples, paste0("P[", par.idx, "]"))
  mcmc_trace(jm.out$jm$samples, paste0("K[", par.idx, "]"))
  
} else {
  mcmc_trace(jm.out$jm$samples, paste0("P[", par.idx, "]"))
  mcmc_trace(jm.out$jm$samples, "K")
}

mcmc_trace(jm.out$jm$samples, paste0("S1[", par.idx, "]"))
mcmc_trace(jm.out$jm$samples, paste0("S2[", par.idx, "]"))
mcmc_trace(jm.out$jm$samples, paste0("Max[", par.idx, "]"))


# plot.trace.dens function is in Granite_Canyon_Counts_fcns.R
# ps.K <- plot.trace.dens(param = "K", 
#                         jags.out = "jm.out$jm")
# 
# ps.P <- plot.trace.dens(param = "P", 
#                         jags.out = "jm.out$jm")
# 
# 
# ps.S1 <- plot.trace.dens(param = "S1", 
#                          jags.out = "jm.out$jm")
# 
# ps.S2 <- plot.trace.dens(param = "S2", 
#                          jags.out = "jm.out$jm")
# 
# ps.Max <- plot.trace.dens(param = "Max", 
#                           jags.out = "jm.out$jm")

# These three year-specific parameters seemed to behave fine.

all.start.year <- c(jm.out$jags.input$jags.input.Laake$all.start.year,
                    jm.out$jags.input$jags.input.new$start.years)

# Create a dataframe with all years, including unsampled years.
all.years <- data.frame(start.year = seq(min(all.start.year), max(all.start.year))) %>%
  mutate(Season = paste0(start.year, "/", start.year + 1))

# Look at the annual abundance estimates:
Nhat. <- data.frame(Season = paste0(all.start.year, "/", all.start.year+1),
                    Nhat = jm.out$jm$q50$Corrected.Est,
                    LCL = jm.out$jm$q2.5$Corrected.Est,
                    UCL = jm.out$jm$q97.5$Corrected.Est) %>%
  right_join(all.years, by = "Season") %>%
  arrange(start.year) %>%
  mutate(Method = "Eguchi")

# This is for daily estimates
N.hats.day <- data.frame(Season = rep(paste0(all.start.year, "/", all.start.year+1), 
                                      each = nrow(jm.out$jm$mean$N)), #rep(Nhat.$Season, each = nrow(jm.out$jm$mean$N)),
                         Day = rep(1:nrow(jm.out$jm$mean$N), times = length(all.start.year)),
                         Mean = as.vector(jm.out$jm$mean$N),
                         LCL = as.vector(jm.out$jm$q2.5$N),
                         UCL = as.vector(jm.out$jm$q97.5$N)) 

# Daily estimates plots
p.daily.Richards <- ggplot(N.hats.day %>% group_by(Season)) + 
  geom_ribbon(aes(x = Day, ymin = LCL, ymax = UCL),
              fill = "blue", alpha = 0.5) +
  geom_path(aes(x = Day, y = Mean)) + 
  geom_point(aes(x = Day, y = Mean)) +
  facet_wrap(~ Season)

# These are not the best estimates because they were not updated as more data
# were collected. I should use the output from the most recent WinBUGS run for 
# the last x years.
Reported.estimates <- read.csv(file = "Data/all_estimates_2024.csv") %>%
  transmute(Season = Season,
            Nhat = Nhat,
            LCL = LCL,
            UCL = UCL,
            Method = paste0(Method, "-Reported")) %>%
  right_join(all.years, by = "Season") %>%
  arrange(start.year) %>%
  relocate(Method, .after = start.year)

WinBugs.run.date <- "2025-04-11"
WinBugs.out <- readRDS(file = paste0("RData/WinBUGS_2007to2025_v2_min", min.dur, 
                                     "_100000_",
                                     WinBugs.run.date, ".rds"))

Corrected.Est <- WinBugs.out$BUGS.out$sims.list$Corrected.Est

# We don't have raw data for 2006/2007 and 2007/2008 seasons
seasons <- c("2006/2007", "2007/2008", jm.out$jags.input$jags.input.new$seasons)

Durban.abundance.df <- data.frame(Season = seasons,
                                  Nhat = apply(Corrected.Est,
                                               FUN = mean,
                                               MARGIN = 2),
                                  # CV = apply(Corrected.Est,
                                  #            FUN = function(x) 100*sqrt(var(x))/mean(x),
                                  #            MARGIN = 2),
                                  # median = apply(Corrected.Est, 
                                  #                FUN = median, 
                                  #                MARGIN = 2),
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
daily.estim <- WinBugs.out$BUGS.out$sims.list$Daily.Est

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

N.hats.day.Durban <- data.frame(Season = rep(seasons, each = dim(daily.estim)[2]),
                                Day = rep(1:dim(daily.estim)[2], length(seasons)),
                                Mean = as.vector(mean.mat),
                                LCL = as.vector(LCL.mat),
                                UCL = as.vector(UCL.mat))

# Daily estimates plots
p.daily.Durban <- ggplot(N.hats.day.Durban %>% group_by(Season)) + 
  geom_ribbon(aes(x = Day, ymin = LCL, ymax = UCL),
              fill = "blue", alpha = 0.5) +
  geom_path(aes(x = Day, y = Mean)) + 
  facet_wrap(~ Season)

# Include non-survey years - no estimates for 2007/2008 because I don't have
# raw data for that year. Only the WinBUGS inputs. 
Laake.abundance.new <- read.csv(file = "Data/all_estimates_Laake_2025.csv") %>%
  mutate(LCL = CL.low,
         UCL = CL.high) %>%
  select(c(Season, Nhat, LCL, UCL)) %>%
  right_join(all.years, by = "Season") %>%
  arrange(start.year) %>%
  mutate(Method = "Laake")

#Laake.output <- read_rds(file = "RData/Laake_abundance_estimates_2024.rds")

# In reported estimates, there are two 2006/2007.
Reported.estimates %>%
  na.omit() %>%
  select(Season) %>% 
  unique() -> sampled.seasons 

Laake.abundance.new %>%
  rbind(Durban.abundance.df) %>%
  rbind(Nhat.) -> all.estimates
#  rbind(spline.Nhat) 
#  rbind(Reported.estimates %>% na.omit()) -> all.estimates

p.Nhats <- ggplot(all.estimates) +
  geom_point(aes(x = start.year, y = Nhat,
                 color = Method),
             alpha = 0.5) +
  geom_errorbar(aes(x = start.year, ymin = LCL, ymax = UCL,
                    color = Method)) +
  ylim(5000, 35000)

# ggsave(plot = p.Nhats,
#        filename = paste0("figures/Nhats_", ver, "_", min.dur, "min.png"),
#        device = "png",
#        dpi = 600)

Nhat. %>% 
  select(Season, start.year, Nhat, LCL, UCL) %>%
  rename(Nhat.Eguchi = Nhat,
         LCL.Eguchi = LCL,
         UCL.Eguchi = UCL) %>%
  cbind(Laake.abundance.new %>%
          select(Nhat, LCL, UCL) %>%
          rename(Nhat.Laake = Nhat,
                 LCL.Laake = LCL,
                 UCL.Laake = UCL)) %>%
  cbind(Durban.abundance.df %>%
          select(Nhat, LCL, UCL) %>%
          rename(Nhat.Durban = Nhat,
                 LCL.Durban = LCL,
                 UCL.Durban = UCL)) %>%
  mutate(d.Laake.Eguchi = Nhat.Laake - Nhat.Eguchi,
         d.Durban.Eguchi = Nhat.Durban - Nhat.Eguchi) -> Nhat.all.wide

# Check convergence
high.Rhat <- function(x){
  return(data.frame(idx = which(x > 1.1),
                    start.year = all.start.year[which(x > 1.1)],
                    Rhat = x[which(x > 1.1)]))
}

high.Rhat.Max <- high.Rhat(jm.out$jm$Rhat$Max)

high.Rhat.K <- high.Rhat(jm.out$jm$Rhat$K)
high.Rhat.S1 <- high.Rhat(jm.out$jm$Rhat$S1)
high.Rhat.S2 <- high.Rhat(jm.out$jm$Rhat$S2)
high.Rhat.P <- high.Rhat(jm.out$jm$Rhat$P)


