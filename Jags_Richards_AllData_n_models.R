#Jags_Richards_AllData_n_models.R
# Runs multiple models and saves results in .rds files
# 


rm(list = ls())

library(tidyverse)
library(loo)

source("Granite_Canyon_Counts_fcns.R")
options(mc.cores = 5)

# Minimum length of observation periods in minutes
min.dur <- 60 #10 #85 #

ver <- c("v1a1", "v2a1", "v3a1", "v4a1", "v5a1", "v6a1", "v7a1", "v8a1" )
Run.date <- Sys.Date()

# These are the ending year of each season - for example, 2022 in the following vector indicates
# for the 2021/2022 season. These data were extracted using Extract_Data_All_v2.Rmd
# Data prior to the 2009/2010 season are in Laake's ERAnalayis package. 
years <- c(2008, 2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024, 2025)
data.dir <- "RData/V2.1_Feb2025"
max.day <- 100

# MCMC.params <- list(n.samples = 550000,
#                     n.thin = 100,
#                     n.burnin = 500000,
#                     n.chains = 5)

#5000 samples
MCMC.params <- list(n.samples = 250000,
                    n.thin = 200,
                    n.burnin = 50000,
                    n.chains = 5)

# 2500 samples
MCMC.params <- list(n.samples = 125000,
                    n.thin = 100,
                    n.burnin = 75000,
                    n.chains = 5)

# 225 samples
# MCMC.params <- list(n.samples = 100,
#                     n.thin = 2,
#                     n.burnin = 10,
#                     n.chains = 5)

jags.params <- c("VS.Fixed", "BF.Fixed",
                 "Max", "K", "K1", "K2", "S1", "S2", "P",
                 "P1", "P2",
                 "mean.prob", "prob", "obs.prob",
                 "mean.N", "Corrected.Est", "N", "obs.N",
                 #"OBS.RF", "sigma.Obs",
                 "Max.alpha", "Max.beta",
                 "S1.alpha", "S2.alpha",
                 "S1.beta", "S2.beta",
                 "mu.P", "rho.P", "sd.proc.P",
                 "mu.log.Max", "rho.Max", "sd.proc.Max",
                 "Raw.Est", "beta.obs",
                 "alpha",
                 #"P.alpha", "P.beta",
                 #"K.alpha", "K.beta",
                 #"beta.1",
                 #"N.alpha", "N.obs",
                 "log.lkhd")
# 
# jags.params <- c("VS.Fixed", "BF.Fixed",
#                  "Max", "K", "K1", "K2", "S1", "S2", "P",
#                  "P1", "P2",
#                  "mean.prob", "prob", "obs.prob",
#                  "mean.N", "Corrected.Est", "N", "obs.N",
#                  "OBS.RF", "sigma.Obs",
#                  "Max.alpha", "Max.beta",
#                  "S1.alpha", "S2.alpha",
#                  "S1.beta", "S2.beta",
#                  "sigma.P", "sigma.Max",
#                  #"P.alpha", "P.beta",
#                  #"K.alpha", "K.beta",
#                  #"beta.1",
#                  #"N.alpha", "N.obs",
#                  "log.lkhd")


for (k in 1:length(ver)){
  jm.out <- NoBUGS_Richards_fcn(min.dur = min.dur, 
                                ver = ver[k], 
                                years = years, 
                                data.dir = "RData/V2.1_Feb2025", 
                                jags.params = jags.params, 
                                MCMC.params = MCMC.params,
                                Run.date = Run.date,
                                obs.n.min = 10,
                                max.day = 100,
                                N.obs = 10,
                                model.name.root = "Richards_HSSM_")
  
}
