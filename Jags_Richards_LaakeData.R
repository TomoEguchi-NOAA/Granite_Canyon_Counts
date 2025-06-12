# Jags_Richards_LaakeData.R
# 
# Uses Laake data and runs
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

Run.date <- Sys.Date()

# Minimum length of observation periods in minutes
# In order to run a new minimum duration, WinBUGS needs to be run first.
min.dur <- 30 #10 #85 #

ver <- "v5"

model.name <- paste0("Richards_pois_bino_", ver) 
jags.model <- paste0("models/model_", model.name, ".txt")

out.file.name <- paste0("RData/JAGS_", model.name,"_min", min.dur,
                        "_LaakeData_",
                        Run.date, ".rds")

jags.input.Laake <- LaakeData2JagsInput(min.dur)

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
  
  jm <- jagsUI::jags(jags.input.Laake$jags.data,
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
                 jags.input = jags.input.Laake,
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

