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

Run.date <- Sys.Date()
#Run.date <- "2024-12-03"

# Minimum length of observation periods in minutes
# In order to run a new minimum duration, WinBUGS needs to be run first.
min.dur <- 30 #10 #85 #

ver <- "v5"

model.name <- paste0("Richards_pois_bino_", ver) 
jags.model <- paste0("models/model_", model.name, ".txt")

out.file.name <- paste0("RData/JAGS_", model.name,"_min", min.dur,
                        "_AllYears_",
                        Run.date, ".rds")

WinBUGS.years <- c(2010, 2011, 2015, 2016, 
                   2020, 2022, 2023, 2024)
all.years <- c(2007, 2008, WinBUGS.years)

jags.input <- AllData2JagsInput(min.dur = min.dur, 
                                WinBUGS.years = WinBUGS.years, 
                                WinBUGS.n.stations = c(1, 1, 2, 2, rep(1, times = 6)), 
                                WinBUGS.out.file = paste0("RData/WinBUGS_",
                                                          min(all.years), "to", 
                                                          max(all.years), "_v2_min",
                                                          min.dur, ".rds"),
                                data.dir = "RData/V2.1_Nov2024")

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

MCMC.params <- list(n.samples = 250000,
                    n.thin = 100,
                    n.burnin = 200000,
                    n.chains = 5)

if (!file.exists(out.file.name)){
  
  Start_Time<-Sys.time()
  
  jm <- jagsUI::jags(jags.input$jags.data,
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
                 jags.input = jags.input,
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
ps.K <- plot.trace.dens(param = "K", 
                        jags.out = "jm.out$jm")

ps.P <- plot.trace.dens(param = "P", 
                        jags.out = "jm.out$jm")


ps.S1 <- plot.trace.dens(param = "S1", 
                         jags.out = "jm.out$jm")

ps.S2 <- plot.trace.dens(param = "S2", 
                         jags.out = "jm.out$jm")

ps.Max <- plot.trace.dens(param = "Max", 
                          jags.out = "jm.out$jm")

