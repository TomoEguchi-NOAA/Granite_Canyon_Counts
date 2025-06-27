# Jags_Richards_Since2010.R
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
# P defines where the peak is relatvie to the range of d min(d) < P < max(d)
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

source("Granite_Canyon_Counts_fcns.R")
#source("AllData2Jags_input.R")
options(mc.cores = 5)

Run.date <- Sys.Date() #"2024-12-05" #

# Minimum length of observation periods in minutes
min.dur <- 60 #85 #30 #85 #

ver <- "v4"

years <- c(2010, 2011, 2015, 2016, 
           2020, 2022, 2023, 2024, 2025)

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
                 "log.lkhd")

MCMC.params <- list(n.samples = 250000,
                    n.thin = 100,
                    n.burnin = 200000,
                    n.chains = 5)

jm.out <- Jags_Richards_Since2010_fcn(min.dur = min.dur, 
                                      max.day = 90, 
                                      ver = ver, 
                                      years = years, 
                                      data.dir = "RData/V2.1_Feb2025", 
                                      jags.params = jags.params, 
                                      MCMC.params = MCMC.params)

# model.name <- paste0("Richards_Nmixture_", ver) 
# out.file.name <- paste0("RData/JAGS_", model.name,"_min", min.dur,
#                         "_Since2010_NoBUGS_",
#                         Run.date, ".rds")
# jm.out <- readRDS(out.file.name)

# need to turn zeros into NAs when there were no second station:
jags.data <- jm.out$jags.input$jags.data
data.array <- jags.data$n
data.array[,2,which(jags.data$n.station == 1)] <- NA

LOOIC.n <- compute.LOOIC(loglik.array = jm.out$jm$sims.list$log.lkhd,
                         data.array = data.array,
                         MCMC.params = jm.out$MCMC.params)

# There are some (< 2%) bad ones. I should look at which ones are not fitting well.
 

# Horwitz-Thompson estimates. This should use sighting probability from 
# estimates. 
# Nhats.HT <- jm.out$jags.data$n[,1,]/jm.out$jags.data$watch.prop[,1,]
# day.idx <- jm.out$jags.data$day[,1,]
# 
# Nhats.HT.day <- matrix(nrow = 94, ncol = dim(Nhats.HT)[2])
# for (k1 in 1:dim(Nhats.HT)[2]){
#   for (k2 in 1:94){
#     Nhats.HT.day[k2, k1] <- sum(Nhats.HT[day.idx[,k1] == k2, k1], na.rm = T)
#   }
# }
#   
# Nhats.HT.all <- data.frame(Season = rep(paste0(all.start.year, "/", all.start.year+1), 
#                                         each = nrow(Nhats.HT.day)),
#                            Day = rep(seq(1, 94), times = ncol(Nhats.HT.day)),
#                            Nhat = as.vector(Nhats.HT.day))
# 

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

# plot.trace.dens function is in Granite_Canyon_Counts_fcns.R
ps.K <- plot.trace.dens(var.name = "\\bK\\b(?!\\.\\w*)", 
                        jm = jm.out$jm)

ps.K. <- plot.trace.dens(var.name = "K.", 
                        jm = jm.out$jm)

ps.P <- plot.trace.dens(var.name = "\\bP\\b(?!\\.\\w*)", 
                        jm = jm.out$jm)

ps.P.alpha <- plot.trace.dens(var.name = "P.alpha", 
                         jm = jm.out$jm)

ps.S1 <- plot.trace.dens(var.name = "S1", 
                         jm = jm.out$jm)

ps.S2 <- plot.trace.dens(var.name = "S2", 
                         jm = jm.out$jm)

ps.Max <- plot.trace.dens(var.name = "Max", 
                          jm = jm.out$jm)


all.start.year <- jm.out$jags.input$start.years

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
            UCL = UCL) %>%
  right_join(all.years, by = "Season") %>%
  arrange(start.year) %>%
  mutate(Method = "Reported")

all.Nhats <- rbind(Nhat., Reported.estimates)

ggplot(all.Nhats) +
  geom_point(aes(x = Season, y = Nhat, color = Method)) +
  geom_errorbar(aes(x = Season, ymin = LCL, ymax = UCL, color = Method))
