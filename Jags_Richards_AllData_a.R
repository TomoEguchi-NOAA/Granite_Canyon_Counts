#Jags_Richards_AllData_a.R
# It runs NoBUGS_Richards_fcn for two (or more) models and saves the results.
# 

rm(list = ls())
#library(ERAnalysis)
library(tidyverse)
library(ggplot2)
library(loo)
library(bayesplot)

source("Granite_Canyon_Counts_fcns.R")
options(mc.cores = 5)

Run.date <- Sys.Date() #"2025-04-21" #"2025-04-17" #

# Minimum length of observation periods in minutes
min.dur <- 60 #10 #85 #

# These are the ending year of each season - for example, 2022 in the following vector indicates
# for the 2021/2022 season. These data were extracted using Extract_Data_All_v2.Rmd
# Data prior to the 2009/2010 season are in Laake's ERAnalayis package. 
years <- c(2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024, 2025)
data.dir <- "RData/V2.1_Feb2025"
max.day <- 100

# MCMC.params <- list(n.samples = 200000,
#                     n.thin = 100,
#                     n.burnin = 150000,
#                     n.chains = 5)

# MCMC.params <- list(n.samples = 350000,
#                     n.thin = 100,
#                     n.burnin = 300000,
#                     n.chains = 5)

MCMC.params <- list(n.samples = 550000,
                    n.thin = 100,
                    n.burnin = 500000,
                    n.chains = 5)

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
# MCMC.params <- list(n.samples = 100,
#                     n.thin = 2,
#                     n.burnin = 50,
#                     n.chains = 5)

jags.params <- c("VS.Fixed", "BF.Fixed",
                 "Max", "K", "K1", "K2", "S1", "S2", "P",
                 "P1", "P2", "P2.U",
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

for (ver in c("v1a", "v2a", "v3a", "v4a", "v5a", "v6a", "v7a", "v8a", "v1b", "v2b", "v3b", "v4b", "v5b", "v6b", "v7b", "v8b")){
#for (ver in c("v1b", "v2b", "v3b", "v4b", "v5b", "v6b", "v7b", "v8b")){
  print(paste0("Starting ", ver, " at ", Sys.time()))
  jm.out <- NoBUGS_Richards_fcn(min.dur = min.dur, 
                                ver = ver, 
                                years = years, 
                                data.dir = "RData/V2.1_Feb2025", 
                                jags.params = jags.params, 
                                MCMC.params = MCMC.params,
                                Run.date = Run.date,
                                obs.n.min = 10,
                                max.day = 100)
}

