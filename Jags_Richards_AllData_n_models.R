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

ver <- c("v5a", "v2a", "v15a", "v16a", "v17a", "v18a", "v19a", "v20a" )
Run.date <- Sys.Date()

# These are the ending year of each season - for example, 2022 in the following vector indicates
# for the 2021/2022 season. These data were extracted using Extract_Data_All_v2.Rmd
# Data prior to the 2009/2010 season are in Laake's ERAnalayis package. 
years <- c(2008, 2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024, 2025)
data.dir <- "RData/V2.1_Feb2025"
max.day <- 100

MCMC.params <- list(n.samples = 550000,
                    n.thin = 100,
                    n.burnin = 500000,
                    n.chains = 5)
# 
# MCMC.params <- list(n.samples = 1000,
#                     n.thin = 10,
#                     n.burnin = 500,
#                     n.chains = 5)

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
                 #"K.alpha", "K.beta",
                 #"beta.1",
                 #"N.alpha", "N.obs",
                 "log.lkhd")


# The following function uses "new" data since 2010 as well as those from Laake's 
# analysis to compute abundance since the 1967/1968 season. There were two seasons
# where the survey continued beyond the 90th day. So, max.day needs to be increased
# from 90. I used 100. 
for (k in 1:length(ver)){
  jm.out <- NoBUGS_Richards_fcn(min.dur = min.dur, 
                                ver = ver[k], 
                                years = years, 
                                data.dir = "RData/V2.1_Feb2025", 
                                jags.params = jags.params, 
                                MCMC.params = MCMC.params,
                                Run.date = Run.date,
                                obs.n.min = 10,
                                max.day = 100)
  
}
