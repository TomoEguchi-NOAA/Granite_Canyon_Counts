#Jags_Richards_NoLaakeData.R
#
# Uses only data since 2006 without data from Laake. It is the
# same as Jags Richards Ver1.Rmd. Some diagnostics are conducted.
# 
# It runs model_Richards_pois_bino_vX, where X is version number. 
# See below.  

# M1 <- (1 + (2 * exp(K) - 1) * exp((1/S1) * (P - d))) ^ (-1/exp(K))
# M2 <- (1 + (2 * exp(K) - 1) * exp((1/S2) * (P - d))) ^ (-1/exp(K))
# N <- min + (max - min) * (M1 * M2)
#
# d is the number of days from the beginning of nesting season
# S1 < 0 and S2 > 0 define the "fatness" of the function
# K > 0 defines the "flatness" at the peak of the function
# P defines where the peak is relative to the range of d min(d) < P < max(d)
# min is "the basal level of nesting outside the nesting season"
# max > min

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

options(mc.cores = 5)

source("Granite_Canyon_Counts_fcns.R")
#source("DataSince2006_Jags_Input.R")

#Run.date <- Sys.Date()
#Run.date <- "2024-12-04"

# Minimum length of observation periods in minutes
min.dur <- 60# 85 #30 #
data.dir <- "RData/V2.1_Feb2025"
max.day <- 100

# This pulls out data from 2010 onward and create jags input
years <- c(2008, 2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024, 2025)
jags.input.new <- data2Jags_input_NoBUGS(min.dur = min.dur, 
                                         years = years,
                                         data.dir = data.dir,
                                         max.day = max.day)

#jags.input.new$max.n <- max(jags.input.new$jags.data$n, na.rm = T)

ver <- "v5a" # "v5a"
model.name <- paste0("Richards_Nmixture_", ver) 
jags.model <- paste0("models/model_", model.name, ".txt")

out.file.name <- paste0("RData/temp/JAGS_", model.name, "_",
                        min(jags.input.new$years), "to", 
                        max(jags.input.new$years), 
                        "_min", min.dur, "NoBUGS.rds")

jags.params <- c("VS.Fixed", "BF.Fixed",
                 "Max", "K", "K1", "K2", "S1", "S2", "P",
                 "P1", "P2",
                 "mean.prob", "prob", "obs.prob",
                 "mean.N", "Corrected.Est", "N", "obs.N",
                 "OBS.RF", "sigma.Obs",
                 "Max.alpha", "Max.beta",
                 "S1.alpha", "S2.alpha",
                 "S1.beta", "S2.beta",
                 #"P.alpha", "P.beta",
                 "K.alpha", "K.beta",
                 #"beta.1",
                 #"N.alpha", "N.obs",
                 "log.lkhd")

# MCMC.params <- list(n.samples = 250000,
#                     n.thin = 100,
#                     n.burnin = 200000,
#                     n.chains = 5)

MCMC.params <- list(n.samples = 550000,
                    n.thin = 100,
                    n.burnin = 500000,
                    n.chains = 5)

if (!file.exists(out.file.name)){
  
  Start_Time<-Sys.time()
  
  jm <- jagsUI::jags(jags.input.new$jags.data,
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
  jm.out <- list(jm = jm,
                 jags.input = jags.input.new,
                 #start.year = all.start.year,
                 jags.params = jags.params,
                 jags.model = jags.model,
                 MCMC.params = MCMC.params,
                 Run_Time = Run_Time,
                 Run_Date = Start_Time,
                 Sys.env = Sys.getenv())
  
  saveRDS(jm.out,
          file = out.file.name)
  
} else {
  
  jm.out <- readRDS(out.file.name)
}

LOOIC.n <- compute.LOOIC(loglik.array = jm.out$jm$sims.list$log.lkhd,
                         MCMC.params = MCMC.params)

# Compute new Rhat (rank-normalized Rhat as per Vehtari et al. 2021)
params <- "^VS\\.Fixed|^BF\\.Fixed|^Max|^S1|^S2|^P1|^P2|^K|^OBS\\.RF\\["
new.Rhat <- rank.normalized.R.hat(jm.out$jm$samples, params)
#max.Rhat <- lapply(.out[[k]]$jm$Rhat, FUN = max, na.rm = T) %>%
#  unlist()
max.new.Rhat.big <- new.Rhat[which(new.Rhat > 1.01)]

if (grepl("a", ver)){
  # mcmc_dens(jm.out$jm$samples, c("S1.alpha", "S1.beta",
  #                                "S2.alpha", "S2.beta"))
  #"P.alpha", "P.beta",
  #"K.alpha", "K.beta"))
  # P.alpha and P.beta seem to be not behaving well - the right tails are not 
  # captured. 
  p.trace.hyper.params <- mcmc_trace(jm.out$jm$samples, c("S1.alpha", "S1.beta",
                                                          "S2.alpha", "S2.beta",
                                                          #"P.alpha", "P.beta",
                                                          "K.alpha", "K.beta"))
  
}

# The following functions are in Granite_Canyon_Counts_fcns.R
K.trace <- K.trace.plots(ver = ver,
                         jm = jm.out$jm,
                         jags.data = jm.out$jags.input$jags.data,
                         new.Rhat = new.Rhat,
                         start.year = jm.out$jags.input$start.years)

P.trace <- P.trace.plots(ver = ver,
                         jm = jm.out$jm,
                         jags.data = jm.out$jags.input$jags.data,
                         new.Rhat = new.Rhat,
                         start.year = jm.out$jags.input$start.years)

S.trace <- S1.S2.trace.plots(ver = ver,
                         jm = jm.out$jm,
                         jags.data = jm.out$jags.input$jags.data,
                         new.Rhat = new.Rhat,
                         start.year = jm.out$jags.input$start.years)

p.trace.Fixed.params <- mcmc_trace(jm.out$jm$samples, c("BF.Fixed", "VS.Fixed"))

par.idx <- c(1:jm.out$jags.input$jags.data$n.year)

# Create a dataframe with all years, including unsampled years.
all.years <- data.frame(start.year = seq(min(jm.out$jags.input$start.years), 
                                         max(jm.out$jags.input$start.years))) %>%
  mutate(Season = paste0(start.year, "/", start.year + 1))

# Look at the annual abundance estimates:
Nhat. <- data.frame(Season = paste0(jm.out$jags.input$start.years, "/", jm.out$jags.input$start.years+1),
                    Nhat = jm.out$jm$q50$Corrected.Est,
                    LCL = jm.out$jm$q2.5$Corrected.Est,
                    UCL = jm.out$jm$q97.5$Corrected.Est) %>%
  right_join(all.years, by = "Season") %>%
  arrange(start.year) %>%
  mutate(Method = paste0("Eguchi ", ver))

# This is for daily estimates
N.hats.day <- data.frame(Season = rep(paste0(jm.out$jags.input$start.years, "/", 
                                             jm.out$jags.input$start.years+1), 
                                      each = nrow(jm.out$jm$mean$N)), 
                         Day = rep(1:nrow(jm.out$jm$mean$N), 
                                   times = length(jm.out$jags.input$start.years)),
                         Mean = as.vector(jm.out$jm$mean$N),
                         LCL = as.vector(jm.out$jm$q2.5$N),
                         UCL = as.vector(jm.out$jm$q97.5$N),
                         Method = "Eguchi") 

# Daily estimates plots
p.daily.Richards <- ggplot(N.hats.day %>% group_by(Season)) + 
  geom_ribbon(aes(x = Day, ymin = LCL, ymax = UCL),
              fill = "blue", alpha = 0.5) +
  geom_path(aes(x = Day, y = Mean)) + 
  #geom_point(aes(x = Day, y = Mean)) +
  facet_wrap(~ Season)

# Bring in WinBUGS results
WinBUGS.Run.Date <- "2025-04-11"
WinBugs.out <- readRDS(file = paste0("RData/WinBUGS_2007to2025_v2_min", min.dur,
                                     "_100000_",
                                     WinBUGS.Run.Date, ".rds"))

# Compute rank-normalized Rhat for WinBUGS output
BUGS.params <- "^lambda|^beta|^OBS.RF"
BUGS.samples <- WinBugs.out$BUGS.out$sims.array
BUGS.col.names <- grep(BUGS.params, dimnames(BUGS.samples)[[3]], 
                       value = T, perl = T)
subset.BUGS.samples <- BUGS.samples[,,BUGS.col.names]
subset.BUGS.mcmc.array <- as_draws_array(subset.BUGS.samples, .nchains = 5)
BUGS.Rhat <- apply(subset.BUGS.mcmc.array, MARGIN = 3, FUN = posterior::rhat)

BUGS.Rhat.big <- BUGS.Rhat[which(BUGS.Rhat > 1.01)]

Corrected.Est <- WinBugs.out$BUGS.out$sims.list$Corrected.Est

# We don't have raw data for 2006/2007 and 2007/2008 seasons
seasons <- c("2006/2007", "2007/2008", jm.out$jags.input$jags.input.new$seasons)

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
                                UCL = as.vector(UCL.mat),
                                Method = "Durban")

# Daily estimates plots
p.daily.Durban <- ggplot(N.hats.day.Durban %>% group_by(Season)) + 
  geom_ribbon(aes(x = Day, ymin = LCL, ymax = UCL),
              fill = "blue", alpha = 0.5) +
  geom_path(aes(x = Day, y = Mean)) + 
  facet_wrap(~ Season)

p.daily <- ggplot(rbind(N.hats.day.Durban, N.hats.day) %>% group_by(Season)) +
  geom_ribbon(aes(x = Day, ymin = LCL, ymax = UCL, fill = Method),
              alpha = 0.5) +
  geom_path(aes(x = Day, y = Mean, color = Method)) + 
  facet_wrap(~ Season)


Durban.abundance.df %>%
  rbind(Nhat.) -> all.estimates

Reported.estimates <- read.csv(file = "Data/Nhats_2025.csv") %>%  
  transmute(Season = Season,
            Nhat = Nhat,
            LCL = LCL,
            UCL = UCL,
            Method = paste0(Method, "-Reported")) %>%
  right_join(all.years, by = "Season") %>%
  arrange(start.year) %>%
  relocate(Method, .after = start.year)

p.Nhats <- ggplot(all.estimates) +
  geom_point(aes(x = start.year, y = Nhat,
                 color = Method),
             alpha = 0.5) +
  geom_errorbar(aes(x = start.year, ymin = LCL, ymax = UCL,
                    color = Method)) +
  ylim(2000, 40000) +
  theme(legend.position = "top")


# 
# # Create a dataframe with all years, including unsampled years.
# all.years <- data.frame(year = seq(min(all.start.year), max(all.start.year))) %>%
#   mutate(Season = paste0(year, "/", year + 1))
# 
# # Look at the annual abundance estimates:
# Nhat. <- data.frame(Season = paste0(all.start.year, "/", all.start.year+1),
#                     Nhat = jm.out$jm$q50$Corrected.Est,
#                     LCL = jm.out$jm$q2.5$Corrected.Est,
#                     UCL = jm.out$jm$q97.5$Corrected.Est) %>%
#   right_join(all.years, by = "Season") %>%
#   arrange(year) %>%
#   mutate(Method = "Richards")
# 
# # This is for daily estimates
# N.hats.day <- data.frame(Season = rep(paste0(all.start.year, "/", all.start.year+1), 
#                                       each = nrow(jm.out$jm$mean$N)), #rep(Nhat.$Season, each = nrow(jm.out$jm$mean$N)),
#                          Day = rep(1:nrow(jm.out$jm$mean$N), times = length(all.start.year)),
#                          Mean = as.vector(jm.out$jm$mean$N),
#                          LCL = as.vector(jm.out$jm$q2.5$N),
#                          UCL = as.vector(jm.out$jm$q97.5$N)) 
# 
# # Daily estimates plots
# p.daily.Richards <- ggplot(N.hats.day %>% group_by(Season)) + 
#   geom_ribbon(aes(x = Day, ymin = LCL, ymax = UCL),
#               fill = "blue", alpha = 0.5) +
#   geom_path(aes(x = Day, y = Mean)) + 
#   geom_point(data = Nhats.HT.all,
#              aes(x = Day, y = Nhat),
#              alpha = 0.3) +
#   facet_wrap(~ Season)
# 
# # These are not the best estimates because they were not updated as more data
# # were collected. I should use the output from the most recent WinBUGS run for 
# # the last x years.
# Reported.estimates <- read.csv(file = "Data/all_estimates_2024.csv") %>%
#   transmute(Season = Season,
#             Nhat = Nhat,
#             LCL = LCL,
#             UCL = UCL,
#             Method = paste0(Method, "-Reported")) %>%
#   right_join(all.years, by = "Season") %>%
#   arrange(year) %>%
#   relocate(Method, .after = year)
# 
# WinBugs.out <- readRDS(file = paste0("RData/WinBUGS_10yr_v2_min",
#                                      min.dur, ".rds"))
# 
# Corrected.Est <- WinBugs.out$BUGS_out$sims.list$Corrected.Est
# seasons <- c("2006/2007", "2007/2008", "2009/2010", "2010/2011", 
#              "2014/2015", "2015/2016", "2019/2020", "2021/2022",
#              "2022/2023", "2023/2024")
# 
# Durban.abundance.df <- data.frame(Season = seasons,
#                                   Nhat = apply(Corrected.Est,
#                                                FUN = mean,
#                                                MARGIN = 2),
#                                   # CV = apply(Corrected.Est,
#                                   #            FUN = function(x) 100*sqrt(var(x))/mean(x),
#                                   #            MARGIN = 2),
#                                   # median = apply(Corrected.Est, 
#                                   #                FUN = median, 
#                                   #                MARGIN = 2),
#                                   LCL = apply(Corrected.Est, 
#                                               MARGIN = 2, 
#                                               FUN = quantile, 0.025),
#                                   UCL = apply(Corrected.Est, 
#                                               MARGIN = 2, 
#                                               FUN = quantile, 0.975)) %>%
#   right_join(all.years, by = "Season") %>%
#   arrange(year) %>%
#   mutate(Method = "Durban")
# 
# # Create a dataframe for daily estimates:
# daily.estim <- WinBugs.out$BUGS_out$sims.list$Daily.Est
# 
# # get stats:
# mean.mat <- LCL.mat <- UCL.mat <- matrix(data = NA, 
#                                          nrow = dim(daily.estim)[2], 
#                                          ncol = dim(daily.estim)[3])
# 
# for (k1 in 1:dim(daily.estim)[2]){
#   for (k2 in 1:dim(daily.estim)[3]){
#     mean.mat[k1, k2] <- mean(daily.estim[,k1,k2])
#     LCL.mat[k1, k2] <- quantile(daily.estim[,k1,k2], 0.025)
#     UCL.mat[k1, k2] <- quantile(daily.estim[,k1,k2], 0.975)
#   }
#   
# }
# 
# N.hats.day.Durban <- data.frame(Season = rep(seasons, each = dim(daily.estim)[2]),
#                                 Day = rep(1:dim(daily.estim)[2], length(seasons)),
#                                 Mean = as.vector(mean.mat),
#                                 LCL = as.vector(LCL.mat),
#                                 UCL = as.vector(UCL.mat))
# 
# # Daily estimates plots
# p.daily.Durban <- ggplot(N.hats.day.Durban %>% group_by(Season)) + 
#   geom_ribbon(aes(x = Day, ymin = LCL, ymax = UCL),
#               fill = "blue", alpha = 0.5) +
#   geom_path(aes(x = Day, y = Mean)) + 
#   geom_point(data = Nhats.HT.all[Nhats.HT.all$Season %in% seasons,],
#              aes(x = Day, y = Nhat),
#              alpha = 0.3) + 
#   facet_wrap(~ Season)
# 
# 
# # In reported estimates, there are two 2006/2007.
# Reported.estimates %>%
#   na.omit() %>%
#   select(Season) %>% 
#   unique() -> sampled.seasons 
# 
# Laake.abundance.new %>%
#   rbind(Durban.abundance.df) %>%
#   rbind(Nhat.) -> all.estimates
# 
# p.Nhats <- ggplot(all.estimates) +
#   geom_point(aes(x = year, y = Nhat,
#                  color = Method),
#              alpha = 0.5) +
#   geom_errorbar(aes(x = year, ymin = LCL, ymax = UCL,
#                     color = Method)) +
#   ylim(0, 45000)
