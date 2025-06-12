# Jags_Spline_Nmixture
# Uses spline fit to the gray whale counts to estimate abundance using
# N-mixture model. 
# 



rm(list = ls())

#library(ERAnalysis)
library(splines)
library(tidyverse)
library(ggplot2)
library(loo)
library(bayesplot)

source("Granite_Canyon_Counts_fcns.R")

WinBUGS.Run.Date <- "2025-04-11"

Run.date <- Sys.Date() #"2025-04-21" #"2025-04-17" #

# Minimum length of observation periods in minutes
min.dur <- 60 #10 #85 #

ver <- "v1" #  "v3" #"v5" # "v4" # 

# These are the ending year of each season - for example, 2022 in the following vector indicates
# for the 2021/2022 season. These data were extracted using Extract_Data_All_v2.Rmd
# Data prior to the 2009/2010 season are in Laake's ERAnalayis package. 
years <- c(2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024, 2025)
data.dir = "RData/V2.1_Feb2025"
max.day = 100

# MCMC.params <- list(n.samples = 250000,
#                     n.thin = 100,
#                     n.burnin = 200000,
#                     n.chains = 5)
# 
# # v3 does not converge well with the above MCMC setting so increasing samples
# MCMC.params <- list(n.samples = 750000,
#                     n.thin = 500,
#                     n.burnin = 500000,
#                     n.chains = 5)

# MCMC.params <- list(n.samples = 25000,
#                     n.thin = 10,
#                     n.burnin = 5000,
#                     n.chains = 5)

MCMC.params <- list(n.samples = 100000,
                    n.thin = 50,
                    n.burnin = 50000,
                    n.chains = 5)

jags.params <- c("VS.Fixed", "BF.Fixed", 
                 "p", "mean.prob",
                 "lambda", "Corrected.Est", "N",
                 "eps", 
                 "log.lkhd")

### Testing model using just new data  ###
# model.name <- paste0("Richards_pois_bino_", ver) 
# jags.model <- paste0("models/model_", model.name, ".txt")
# 
jags.data.list <- data2Jags_input_NoBUGS(min.dur = 60,
                                         years = years,
                                         data.dir = "RData/V2.1_Feb2025",
                                         max.day = 100)

jags.data <- jags.data.list$jags.data
jags.data$B <- as.matrix(bs(c(1:100), df = 14, intercept = TRUE) ) %>%
  labelled::remove_attributes(c("dimnames", "degree", "knots", "Boundary.knots", "intercept"))
jags.data$nbasis <- ncol(jags.data$B)

jags.data.1 <- list(nyears = jags.data$n.year,
                    ntime = jags.data$n.days,
                    B = jags.data$B,
                    nrep = jags.data$periods[,1],
                    ycount = jags.data$n[,1,],
                    vs = jags.data$vs[,1,],
                    bf = jags.data$bf[,1,],
                    obs = jags.data$obs[,1,],
                    K = ncol(jags.data$B))

jags.params.1 <- c("beta", "alpha0", "alpha1", "N")

# n <- jags.data$n
# 
# for (k in 1:dim(n)[3]){
#   for (k2 in 1:jags.data$n.station[k]){
#     n[(jags.data$periods[k, k2]-1):jags.data$periods[k, k2], k2, k] <- c(NA, NA)
#     
#   }
# }
# jags.data$n <- NULL
# jags.data$n <- n
jags.data$N <- NULL

jm <- jagsUI::jags(jags.data.1,
                   inits = NULL,
                   parameters.to.save= jags.params,
                   model.file = "models/model_Jags_Spline_v1_1.txt",
                   n.chains = MCMC.params$n.chains,
                   n.burnin = MCMC.params$n.burnin,
                   n.thin = MCMC.params$n.thin,
                   n.iter = MCMC.params$n.samples,
                   DIC = T,
                   parallel=T)
# 
# jm.out <- list(jm = jm,
#                jags.input = jags.data.list,
#                #start.year = all.start.year,
#                jags.params = jags.params,
#                jags.model = jags.model,
#                MCMC.params = MCMC.params,
#                Sys.env = Sys.getenv())


# jags.data$n[1:10, 1, 1]
# [1] 1 7 0 2 1 3 3 2 2 0
# > jags.data$watch.prop[1:6, 1, 1]
# [1] 0.1675911 0.1680544 0.1651222 0.1657411 0.1666656 0.1665111
# > sum(jags.data$watch.prop[1:6, 1, 1])
# [1] 0.9996856
# > sum(jm$mean$N[1:6, 1, 1])
# [1] 260.0753
# > sum(jags.data$watch.prop[1:6, 1, 1]) * 271
# [1] 270.9148
# > jags.data$watch.prop[1:6, 1, 1] * sum(jm$mean$N[1:6, 1, 1])
# jm$mean$prob[1:6, 1, 1]
###  End of testing model using just new data  ###