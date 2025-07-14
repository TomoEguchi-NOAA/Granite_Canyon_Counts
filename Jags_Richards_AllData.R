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

# Model versions - depending on which parameters are constant/year-specific
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
options(mc.cores = 5)

WinBUGS.Run.Date <- "2025-04-11"
#WinBUGS.Run.Date <- "2025-06-06"

Run.date <- Sys.Date() #"2025-04-21" #"2025-04-17" #

# Minimum length of observation periods in minutes
min.dur <- 60 #10 #85 #

# v3 has conversion issues. v4 and v5 seem to work fine. 2025-06-25
ver <- "v5" #  "v3" #"v5" # "v4" # 

# These are the ending year of each season - for example, 2022 in the following vector indicates
# for the 2021/2022 season. These data were extracted using Extract_Data_All_v2.Rmd
# Data prior to the 2009/2010 season are in Laake's ERAnalayis package. 
years <- c(2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024, 2025)
data.dir <- "RData/V2.1_Feb2025"
max.day <- 100

MCMC.params <- list(n.samples = 250000,
                    n.thin = 100,
                    n.burnin = 200000,
                    n.chains = 5)
# 
# # v3 does not converge well with the above MCMC setting so increasing samples
# MCMC.params <- list(n.samples = 1000000,
#                     n.thin = 500,
#                     n.burnin = 750000,
#                     n.chains = 5)

# MCMC.params <- list(n.samples = 50000,
#                     n.thin = 20,
#                     n.burnin = 10000,
#                     n.chains = 5)

# MCMC.params <- list(n.samples = 25000,
#                     n.thin = 10,
#                     n.burnin = 5000,
#                     n.chains = 5)
# 
# MCMC.params <- list(n.samples = 10000,
#                     n.thin = 10,
#                     n.burnin = 500,
#                     n.chains = 5)

jags.params <- c("VS.Fixed", "BF.Fixed",
                 "Max", "K", "S1", "S2", "P",
                 "mean.prob", "prob", "obs.prob",
                 "mean.N", "Corrected.Est", "N", "obs.N",
                 "OBS.RF", "sigma.Obs",
                 "Max.alpha", "Max.beta",
                 "S1.alpha", "S2.alpha",
                 "S1.beta", "S2.beta",
                 "P.alpha", "P.beta",
                 "K.alpha", "K.beta",
                 #"beta.1",
                 #"N.alpha", "N.obs",
                 "log.lkhd")

# jags.params <- c("VS.Fixed", "BF.Fixed",
#                  "mean.prob", "prob", 
#                  "Corrected.Est", "Obs.N", "N",
#                  "beta0", "mu", "sigma1", "sigma2",
#                  "OBS.RF", "log.lkhd")

###############################
### Testing model using just new data  ###
# model.name <- paste0("Richards_pois_bino_", ver)
# jags.model <- paste0("models/model_", model.name, ".txt")
#jags.model <- "models/model_Split_Gaussian_N_Mixture_v3.txt"

# jags.data.list <- data2Jags_input_NoBUGS(min.dur = 60,
#                                          years = years,
#                                          data.dir = "RData/V2.1_Feb2025",
#                                          max.day = 100)

# Remove day = 1 and day = 100, rearrange, put 1 on top, 100 at the bottom 
# jags.data <- jags.data.list$jags.data
# jags.data$periods <- jags.data$periods - 2
# #
# jags.data$day[jags.data$day == 1] <- NA
# jags.data$day[jags.data$day == max.day] <- NA
# 
# jags.data$n[jags.data$day == 1] <- NA
# jags.data$n[jags.data$day == max.day] <- NA
# 
# for (k in 1:jags.data$n.year){
#   jags.data$day[(jags.data$periods[k, 1]+1), 1, k] <- 100
#   jags.data$n[(jags.data$periods[k, 1]+1), 1, k] <- 0
#   if (jags.data$n.station[k] == 2){
#     jags.data$day[(jags.data$periods[k, 2]+1), 2, k] <- 100
#     jags.data$n[(jags.data$periods[k, 2]+1), 2, k] <- 0
#   }
# 
# }

# jags.data$day <- abind::abind(array(data = 1,
#                                     dim = c(1, 2, jags.data$n.year)),
#                               jags.data$day, along = 1)
# 
# jags.data$n <- abind::abind(array(data = 0,
#                                   dim = c(1, 2, jags.data$n.year)),
#                               jags.data$n, along = 1)
# 
# jags.data$scaled.day <- jags.data$day-(max.day/2)
# jags.data$x_day <- matrix(data = NA, nrow = max.day, ncol = dim(jags.data$n)[3])
# for (k in 1:dim(jags.data$n)[3])
#   jags.data$x_day[,k] <- (seq(1, max.day))/max.day

# pool observers with < 20 observation periods
# obs.vec <- as.vector(jags.data$obs)
# data.frame(obs = obs.vec) %>%
#   mutate(obs.f = as.factor(obs)) %>%
#   group_by(obs.f) %>%
#   summarize(n = n(),
#             obs = first(obs)) -> obs.summary
# 
# obs.n.min <- 10
# obs.too.few <- obs.summary %>% filter(n < obs.n.min)
# 
# obs.to.keep <- obs.summary %>% filter(n >= obs.n.min)
# obs.to.keep$new.ID <- seq(1, dim(obs.to.keep)[1])
# obs.others <- max(obs.to.keep$new.ID)
# obs <- jags.data$obs
# new.no.obs <- obs.others + 1
# old.no.obs <- max(obs.to.keep$obs)
# for (k in 1:nrow(obs.too.few)){
#   obs[obs == obs.too.few$obs[k]] <- NA
# }
# 
# for (k in 1:(nrow(obs.to.keep)-1)){
#   obs[obs == obs.to.keep$obs[k]]  <- obs.to.keep$new.ID[k]
# }
# 
# obs[is.na(obs)] <- obs.others
# obs[obs == old.no.obs] <- new.no.obs
# 
# jags.data$obs <- obs
# jags.data$n.obs <- max(obs) - 1

# min.N <- matrix(nrow = max.day, ncol = dim(jags.data$n)[3])
# 
# k <- 1
# d <- 44
# for (k in 1:dim(jags.data$n)[3]){
#   for (d in 1:max.day){
#     n.sum.1 <- sum(jags.data$n[jags.data$day[,1,k] == d,1,k], na.rm = T)    
#     if (jags.data$n.station[k] == 2){
#       n.sum.2 <- sum(jags.data$n[jags.data$day[,2,k] == d, 2, k], na.rm = T)
#     }
#     
#     min.N[d,k] <- max(n.sum.1, n.sum.2)
#   }
# }
# 
# jags.data$min.N <- min.N

#jags.data["N"] <- NULL
#jags.data$n.year <- 5

# jm <- jagsUI::jags(jags.data,
#                    inits = NULL,
#                    parameters.to.save= jags.params,
#                    model.file = jags.model,
#                    n.chains = MCMC.params$n.chains,
#                    n.burnin = MCMC.params$n.burnin,
#                    n.thin = MCMC.params$n.thin,
#                    n.iter = MCMC.params$n.samples,
#                    DIC = T,
#                    parallel=T)
# 
# jm.out <- list(jm = jm,
#                jags.input = jags.data.list,
#                #start.year = all.start.year,
#                jags.params = jags.params,
#                jags.model = jags.model,
#                MCMC.params = MCMC.params,
#                Sys.env = Sys.getenv())


###############################

# The following function uses "new" data since 2010 as well as those from Laake's 
# analysis to compute abundance since the 1967/1968 season. There were two seasons
# where the survey continued beyond the 90th day. So, max.day needs to be increased
# from 90. I used 100. 
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
data.array[,2,which(jm.out$jags.input$jags.data$n.station == 1)] <- NA

LOOIC.n <- compute.LOOIC(loglik.array = jm.out$jm$sims.list$log.lkhd,
                         data.array = data.array,
                         MCMC.params = MCMC.params)

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
  mutate(Method = paste0("Eguchi ", ver))

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
  #geom_point(aes(x = Day, y = Mean)) +
  facet_wrap(~ Season)

# These are not the best estimates because they were not updated as more data
# were collected. I should use the output from the most recent WinBUGS run for 
# the last x years.
#Reported.estimates <- read.csv(file = "Data/all_estimates_2024.csv") %>%
Reported.estimates <- read.csv(file = "Data/Nhats_2025.csv") %>%  
  transmute(Season = Season,
            Nhat = Nhat,
            LCL = LCL,
            UCL = UCL,
            Method = paste0(Method, "-Reported")) %>%
  right_join(all.years, by = "Season") %>%
  arrange(start.year) %>%
  relocate(Method, .after = start.year)

#WinBugs.run.date <- "2025-04-11"
WinBugs.out <- readRDS(file = paste0("RData/WinBUGS_2007to2025_v2_min", min.dur,
                                     "_100000_",
                                     WinBUGS.Run.Date, ".rds"))

# WinBugs.out <- readRDS(file = paste0("RData/WinBUGS_1968to2025_v2_min", min.dur, 
#                                      "_85000_",
#                                      WinBUGS.Run.Date, ".rds"))

Corrected.Est <- WinBugs.out$BUGS.out$sims.list$Corrected.Est

# We don't have raw data for 2006/2007 and 2007/2008 seasons
seasons <- c("2006/2007", "2007/2008", jm.out$jags.input$jags.input.new$seasons)

all.season <- paste0(all.start.year, "/", all.start.year+1)
Durban.abundance.df <- data.frame(Season = WinBugs.out$BUGS.input$seasons,
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

N.hats.day.Durban <- data.frame(Season = rep(WinBugs.out$BUGS.input$seasons, 
                                             each = dim(daily.estim)[2]),
                                Day = rep(1:dim(daily.estim)[2], 
                                          length(WinBugs.out$BUGS.input$seasons)),
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
  ylim(2000, 40000) +
  theme(legend.position = "top")

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

# Compare how daily sums among years
obsd.periods.primary <- jm.out$jags.input$jags.data$periods[,1]
watch.prop.primary <- jm.out$jags.input$jags.data$watch.prop[,1,]
obsd.effort.primary <- rbind(rep(0, times = dim(watch.prop.primary)[2]), 
                             watch.prop.primary, 
                             rep(0, times = dim(watch.prop.primary)[2]))

obsd.n.primary <- jm.out$jags.input$jags.data$n[,1,]
obsd.day.primary <- jm.out$jags.input$jags.data$day[,1,]
obsd.n.prop <- obsd.n.primary[,] * obsd.effort.primary

obsd.n.df <- data.frame(Season = rep(all.season, each = dim(obsd.n.prop)[1]),
                        obsd.n = as.vector(obsd.n.prop),
                        day = as.vector(obsd.day.primary),
                        effort = as.vector(obsd.effort.primary)) %>%
  na.omit()


ggplot(obsd.n.df) +
  geom_point(aes(x = day, y = obsd.n)) +
  facet_wrap(~ Season)

ggplot(obsd.n.df) +
  geom_point(aes(x = day, y = effort)) +
  facet_wrap(~ Season)

# Simple comparison between observed counts per hour vs. estimated abundance
# obsd.n.prop.sum <- data.frame(Season = all.season,
#                               n.sum = colSums(obsd.n.prop, na.rm = T))
# 
# Nhat.all.wide %>% 
#   left_join(obsd.n.prop.sum, by = "Season") -> Nhat.all.wide 
# 
# ggplot(Nhat.all.wide) +
#   geom_point(aes(y = Nhat.Laake, x = n.sum), color = "blue") +
#   geom_point(aes(y = Nhat.Eguchi, x = n.sum), color = "red")
# 
# Nhat.all.wide %>%
#   filter(n.sum < 1000) -> Nhat.all.wide.1000
# 
# ggplot(Nhat.all.wide.1000) +
#   geom_point(aes(y = Nhat.Laake, x = n.sum), color = "blue") +
#   geom_point(aes(y = Nhat.Eguchi, x = n.sum), color = "red")
