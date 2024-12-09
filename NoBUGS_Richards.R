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

min.dur <- 85
ver <- "v5"
Run.date <- Sys.Date()
data.dir <- "RData/V2.1_Nov2024"

model.name <- paste0("Richards_pois_bino_", ver) 
jags.model <- paste0("models/model_", model.name, ".txt")

out.file.name <- paste0("RData/JAGS_", model.name,"_min", min.dur,
                        "_NoBUGS_",
                        Run.date, ".rds")
                        
jags.input.Laake <- LaakeData2JagsInput(min.dur)

years <- c(2010, 2011, 2015, 2016, 
           2020, 2022, 2023, 2024)

jags.input.new <- data2Jags_input_NoBUGS(min.dur = min.dur, 
                                         years = years,
                                         data.dir = data.dir)

# Combine all observers from both datasets and give them new ID numbers
obs.Laake <- jags.input.Laake$obs %>%
  mutate(data = "Laake") %>%
  rename(ID = obs.ID)

# In Laake's dataset, there is no ID for "No observers" and NAs are used
# to fill in when there were no observers. 
obs.Laake <- rbind(obs.Laake, 
                   data.frame(obs = "No obs", 
                              ID = nrow(obs.Laake)+1, 
                              data = "Laake"))

# Note the no observer ID has been changed from 36 to whatever the number
# that was assigned
obs.new <- jags.input.new$obs %>%
  mutate(data = "new")

obs.Laake %>% 
  full_join(obs.new, by = "obs") %>%
  rename(ID.Laake = ID.x) -> obs.all

obs.all %>%
  mutate(ID.all = seq(1:nrow(obs.all))) %>%
  select(obs, ID.Laake, ID.new, ID.all) -> obs.all.new

# Replace observer IDs in each dataset with new IDs (obs.all.new$ID.all)
jags.obs.Laake <- jags.input.Laake$jags.data$obs

# Replace NAs with no observer ID
jags.obs.Laake[is.na(jags.obs.Laake)] <- obs.Laake[nrow(obs.Laake), "ID"]

# use "replace" to swap the original observer IDs with new IDs
# using the look up table (obs.all.new)
# c(new, x)[match(x, c(old, x))], where x is data
# https://stackoverflow.com/questions/16228160/multiple-replacement-in-r

No.obs.ID <- obs.all.new %>% 
  filter(obs == "No obs") %>% 
  select(ID.all) %>% pull()

jags.obs.Laake.new <- array(data = No.obs.ID,
                            dim = dim(jags.obs.Laake))

jags.obs.Laake.new[,1,] <- apply(jags.obs.Laake[,1,], 
                                 FUN = function(x) c(as.vector(obs.all.new$ID.all), 
                                                     x)[match(x, as.vector(obs.all.new$ID.Laake), x)], 
                                 MARGIN = 2)

jags.obs.Laake.new[,2,] <- apply(jags.obs.Laake[,2,], 
                                 FUN = function(x) c(as.vector(obs.all.new$ID.all), 
                                                     x)[match(x, as.vector(obs.all.new$ID.Laake), x)], 
                                 MARGIN = 2)

jags.obs.new <- jags.input.new$jags.data$obs

jags.obs.new.new <- array(data = No.obs.ID,
                          dim = dim(jags.obs.new))

jags.obs.new.new[,1,] <- apply(jags.obs.new[,1,], 
                                 FUN = function(x) c(as.vector(obs.all.new$ID.all), 
                                                     x)[match(x, as.vector(obs.all.new$ID.new), x)], 
                                 MARGIN = 2)

jags.obs.new.new[,2,] <- apply(jags.obs.new[,2,], 
                                 FUN = function(x) c(as.vector(obs.all.new$ID.all), 
                                                     x)[match(x, as.vector(obs.all.new$ID.new), x)], 
                                 MARGIN = 2)

# all data arrays need to be combined between Laake and new datasets
# First make the first dimension the same size between the two arrays:
n.row.dif <- dim(jags.obs.Laake.new)[1] - dim(jags.obs.new.new)[1]

# observers:
jags.obs.new.new.1 <- abind(jags.obs.new.new,
                            array(data = No.obs.ID, 
                                  dim = c(n.row.dif, 2, 
                                          dim(jags.obs.new.new)[3])),
                            along = 1)

# Combine the arrays.
jags.obs.all <- abind(jags.obs.Laake.new, 
                      jags.obs.new.new.1,
                      along = 3)

# whale counts:
jags.n.new <- abind(jags.input.new$jags.data$n,
                    array(data = NA, 
                          dim = c(n.row.dif, 2, 
                                  dim(jags.obs.new.new)[3])),
                    along = 1)

jags.n.all <- abind(jags.input.Laake$jags.data$n,
                    jags.n.new,
                    along = 3)

# Watch length
jags.watch.prop.new <- abind(jags.input.new$jags.data$watch.prop,
                             array(data = NA, 
                                   dim = c(n.row.dif, 2, 
                                           dim(jags.obs.new.new)[3])),
                             along = 1)

jags.watch.prop.all <- abind(jags.input.Laake$jags.data$watch.prop,
                             jags.watch.prop.new,
                             along = 3)

# day
jags.day.new <- abind(jags.input.new$jags.data$day,
                      array(data = NA, 
                            dim = c(n.row.dif, 2, 
                                    dim(jags.obs.new.new)[3])),
                      along = 1)

jags.day.all <- abind(jags.input.Laake$jags.data$day,
                             jags.day.new,
                             along = 3)

# visibility
# Visibility and Beaufort do not have extra two rows. So need to recalculate
# the difference in the numbers of rows
n.row.dif <- dim(jags.input.Laake$jags.data$vs)[1] - dim(jags.input.new$jags.data$vs)[1]
jags.vs.new <- abind(jags.input.new$jags.data$vs,
                    array(data = NA, 
                          dim = c(n.row.dif, 2, 
                                  dim(jags.obs.new.new)[3])),
                    along = 1)

jags.vs.all <- abind(jags.input.Laake$jags.data$vs,
                     jags.vs.new,
                     along = 3)

# Beaufort
jags.bf.new <- abind(jags.input.new$jags.data$bf,
                     array(data = NA, 
                           dim = c(n.row.dif, 2, 
                                   dim(jags.obs.new.new)[3])),
                     along = 1)

jags.bf.all <- abind(jags.input.Laake$jags.data$bf,
                     jags.bf.new,
                     along = 3)

jags.data <- list(  n = jags.n.all,
                    n.station = c(jags.input.Laake$jags.data$n.station,
                                  jags.input.new$jags.data$n.station),
                    n.year = dim(jags.bf.all)[3],
                    n.obs = nrow(obs.all.new),
                    periods = rbind(as.matrix(jags.input.Laake$jags.data$periods), 
                                    jags.input.new$jags.data$periods),
                    n.days = max(jags.input.Laake$jags.data$n.days,
                                 jags.input.new$jags.data$n.days),
                    obs = jags.obs.all,
                    vs = jags.vs.all,
                    bf = jags.bf.all,
                    watch.prop = jags.watch.prop.all,
                    day = jags.day.all)


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

jags.input <- list(jags.data = jags.data,
                   min.dur = min.dur, 
                   jags.input.Laake = jags.input.Laake,
                   jags.input.new = jags.input.new,
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


