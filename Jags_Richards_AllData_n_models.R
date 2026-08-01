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
data.dir <- "RData/V2.1_May2026"
max.day <- 100

#5000 samples
MCMC.params <- list(n.samples = 250000,
                    n.thin = 200,
                    n.burnin = 50000,
                    n.chains = 5)

# 225 samples
# MCMC.params <- list(n.samples = 100,
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
                 "alpha", "r", "kappa",
                 #"P.alpha", "P.beta",
                 #"K.alpha", "K.beta",
                 #"beta.1",
                 #"N.alpha", "N.obs",
                 "log.lkhd")

jags.input <- NoBUGS_Jags_input(min.dur, years, data.dir, max.day)
jags.data <- jags.input$jags.data

model.names <- list()
for (k in 1:nrow(model.defs)){
  model.names[[k]] <- Richards_HSSM_model_definition(K = 1,
                                                     S1 = model.defs[k, "S1"],
                                                     S2 = model.defs[k, "S2"],
                                                     P = "time",
                                                     Max = "time",
                                                     lkhd = model.defs[k, "Lkhd"])
  
  model.name.no.dir <- strsplit(model.names[[k]], split = "models/")[[1]][2]
  
  if (length(grep("M1a", model.names[[k]])) > 0){
    S1.length <- jags.data$n.year
    S2.length <- jags.data$n.year
  }
  
  if (length(grep("M2a", model.names[[k]])) > 0){
    S1.length <- 1
    S2.length <- 1
  }

  if (length(grep("M3a", model.names[[k]])) > 0){
    S1.length <- jags.data$n.year
    S2.length <- 1
  }
  
  if (length(grep(c("M4a"), model.names[[k]])) > 0){
    S1.length <- 1
    S2.length <- jags.data$n.year
  } 
  
  make_inits <- function(seed) list(
    beta0.P = runif(1, 35, 55),  
    beta1.P = rnorm(1, 0, 1),  
    sd.proc.P = runif(1, 0.5, 3),
    P = runif(jags.data$n.year, 35, 55),
    beta0.Max = rnorm(1, 7.6, 0.7), 
    beta1.Max = rnorm(1, 0, 0.5), 
    sd.proc.Max = runif(1, 0.2, 1),
    log.Max = rnorm(jags.data$n.year, 7.6, 0.7),
    S1 = runif(S1.length, 1, 12), 
    S2 = runif(S2.length, 1, 12),   # wide -> spreads chains across the sharp/gradual basins
    S1.alpha = runif(1, 5, 15), 
    S1.beta = runif(1, 0.5, 2),
    S2.alpha = runif(1, 5, 15), 
    S2.beta = runif(1, 0.5, 2),
    r = runif(1, 1, 20),
    BF.Fixed = rnorm(1, 0, 1), 
    VS.Fixed = rnorm(1, 0, 1),
    sd.obs = runif(1, 0.3, 1.2), 
    alpha = rnorm(jags.data$n.obs.fixed, 1.39, 0.5),
    .RNG.name = "base::Mersenne-Twister", .RNG.seed = seed
  )
  
  inits <- lapply(1:MCMC.params$n.chains, function(i) make_inits(1000 + i))
  
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
                                ext = ".jags",
                                inits = inits)
  
}




