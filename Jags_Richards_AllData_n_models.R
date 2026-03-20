#Jags_Richards_AllData_n_models.R
# Runs multiple models and saves results in .rds files
# 


rm(list = ls())

library(tidyverse)
library(loo)

source("Granite_Canyon_Counts_fcns.R")
source("Richards_HSSM_model_definition.R")
options(mc.cores = parallel::detectCores())

# Minimum length of observation periods in minutes
min.dur <- 60 #10 #85 #

model.defs <- data.frame(Lkhd = c(rep("NegBin", 4), rep("Poisson", 4)),
                         P = "time",
                         S1 = rep(c("time", "S1"), 4),
                         S2 = rep(c("time", "S2", "S2", "time"), 2))

Run.date <- Sys.Date()

# These are the ending year of each season - for example, 2022 in the following vector indicates
# for the 2021/2022 season. These data were extracted using Extract_Data_All_v2.Rmd
# Data prior to the 2009/2010 season are in Laake's ERAnalayis package. 
years <- c(2008, 2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024, 2025, 2026)
data.dir <- "RData/V2.1_Feb2026"
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
# MCMC.params <- list(n.samples = 125000,
#                     n.thin = 100,
#                     n.burnin = 75000,
#                     n.chains = 5)

# 225 samples
# MCMC.params <- list(n.samples = 100,
#                     n.thin = 2,
#                     n.burnin = 10,
#                     n.chains = 5)

# MCMC.params <- list(n.samples = 500,
#                     n.thin = 2,
#                     n.burnin = 10,
#                     n.chains = 5)

jags.params <- c("VS.Fixed", "BF.Fixed",
                 "Max", "K", "K1", "K2", "S1", "S2", "P",
                 "mean.prob", "prob", "obs.prob",
                 "mean.N", "Corrected.Est", "N", "obs.N",
                 #"OBS.RF", "sigma.Obs",
                 "Max.alpha", "Max.beta",
                 "S1.alpha", "S2.alpha",
                 "S1.beta", "S2.beta",
                 "beta0.P", "beta1.P", "sd.proc.P",
                 "beta.p",  "sd.obs",
                 "beta0.Max", "beta1.Max", "sd.proc.Max",
                 "Raw.Est", "beta.obs",
                 "alpha", "r",
                 #"P.alpha", "P.beta",
                 #"K.alpha", "K.beta",
                 #"beta.1",
                 #"N.alpha", "N.obs",
                 "log.lkhd")
model.names <- list()
for (k in 1:nrow(model.defs)){
  model.names[[k]] <- Richards_HSSM_model_definition(K = 1,
                                                     S1 = model.defs[k, "S1"],
                                                     S2 = model.defs[k, "S2"],
                                                     P = "time",
                                                     Max = "time",
                                                     lkhd = model.defs[k, "Lkhd"])
  
  model.name.no.dir <- strsplit(model.names[[k]], split = "models/")[[1]][2]
  
  jm.out <- NoBUGS_Richards_fcn(min.dur = min.dur, 
                                years = years, 
                                data.dir = data.dir, 
                                jags.params = jags.params, 
                                MCMC.params = MCMC.params,
                                Run.date = Run.date,
                                obs.n.min = 10,
                                max.day = 100,
                                N.obs = 10,
                                model.name = model.name.no.dir,
                                ext = ".jags")
  
}
