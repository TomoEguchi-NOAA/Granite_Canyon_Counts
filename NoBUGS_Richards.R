# NoBUGS_Richards.R
# 
# Runs Richards function analysis without using WinBUGS input.

rm(list=ls())
library(tidyverse)
library(abind)
library(ggplot2)
library(loo)
library(bayesplot)

source("Granite_Canyon_Counts_fcns.R")

min.dur <- 10 #85 in minutes. Need to run Extract_Data_All_v2.Rmd with the min.dur before using this. 
ver <- "v5"
Run.date <- Sys.Date()
data.dir <- "RData/V2.1_Nov2024"

model.name <- paste0("Richards_pois_bino_", ver) 
jags.model <- paste0("models/model_", model.name, ".txt")

out.file.name <- paste0("RData/JAGS_", model.name,"_min", min.dur,
                        "_NoBUGS_",
                        Run.date, ".rds")

new.years <- c(2010, 2011, 2015, 2016,
               2020, 2022, 2023, 2024)

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

jags.input.list <- AllData2JagsInput_NoBUGS(min.dur, years = new.years, data.dir)   
jags.data <- jags.input.list$jags.data

jags.input <- list(jags.data = jags.input.list$jags.data,
                   min.dur = min.dur, 
                   jags.input.Laake = jags.input.list$jags.input.Laake,
                   jags.input.new = jags.input.list$jags.input.new,
                   data.dir = data.dir)

if (!file.exists(out.file.name)){
  
  Start_Time<-Sys.time()
  
  jm <- jagsUI::jags(jags.data,
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


