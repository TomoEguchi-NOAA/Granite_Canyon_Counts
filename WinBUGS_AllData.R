# WinBUGS_AllData.R
# 
# Combines Laake data and more recent data and runs
# GW_Nmix_Orig.bugs


rm(list = ls())

library(R2WinBUGS)
library(tidyverse)
library(bayesplot)

WinBUGS.dir <- paste0(Sys.getenv("HOME"), "/WinBUGS14")

source("Granite_Canyon_Counts_fcns.R")

Run.date <- Sys.Date() #"2025-04-17" #

# Minimum length of observation periods in minutes
# In order to run a new minimum duration
min.dur <- 60 #10 #85 #

# These are the ending year of each season - different from the standard, which
# uses the beginning year. So, for example, 2022 in the following vector indicates
# for the 2021/2022 season. These are the "new" data.
# min.dur is one of 10, 30, 85 or any other that was used in Extract_Data_All_v2.Rmd
years <- c(2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024, 2025)
data.dir = "RData/V2.1_Feb2025"
input.list <- AllData2JagsInput_NoBUGS(min.dur, years = years, data.dir)  

# Convert jags input to WinBUGS input. Necessary variables are:
# n = n[1:max(periods[1:x]),,1:x],
# n.com = n[1:max(periods[1:x]),,1:x],
# n.sp = n[1:max(periods[1:x]),,1:x],
#### The above three are the same.
# 
# n.station = dim(n[1:max(periods[1:x]),,1:x])[2],
# n.year = dim(n[1:max(periods[1:x]),,1:x])[3],
# n.obs = max(obs[1:max(periods[1:x]),,1:x], na.rm = T),
# periods = periods[1:x],
# obs = obs[1:max(periods[1:x]),,1:x],
# u = u[1:max(periods[1:x]),,1:x],
#### u is an index array for whether or not observers were present for each period 
# 
# vs = vs[1:max(periods[1:x]),1:x],
# bf = bf[1:max(periods[1:x]),1:x],

# day = t[1:(max(periods[1:x])+2),1:x],
# N = N[,1:x],
# N.com = N[,1:x],
# N.sp = N[,1:x],
# ### The above three lines are the same
# 
# knot = c(-1.46,-1.26,-1.02,-0.78,
#          -0.58,-0.34,-0.10,0.10,
#          0.34,0.57,0.78,1.02,1.26,1.46),
# n.knots=14,
# Watch.Length=Watch.Length[1:(max(periods[1:x])+2), 1:x]

# For this analysis, I have to change "days" array because 90 was not the maximum
# in the old (Laake's) datasets. I truncate data to maximum of day 90. 
day.mat <- input.list$jags.data$day[,1,]
day.mat[day.mat > 89] <- NA   # Remove all 90 and above
day.mat[day.mat < 2] <- NA    # Remove all 1s

periods <- colSums(!is.na(day.mat))

# Create the u array by replacing the non-observer ID (n.obs + 1) with 0s
# and observers with 1s
u <- input.list$jags.data$obs
u[u > max(input.list$jags.data$n.obs)] <- 0
u[u != 0] <- 1

#The 'data' has to be the inverse of the inits, 
# with NAs for all of the estimated Ns, and 0s for the days 1 and 90
N <- matrix(NA, 
            nrow = max(periods) + 2, 
            ncol = length(periods)) 

for(i in 1:length(periods)){
  N[(periods[i]+1):(periods[i]+2),i] <- 0 #True number of whales passing fixed at 0 for day 1 and 90
  day.mat[(periods[i]+1):(periods[i]+2),i] <- c(1, 90)
}

Watch.Length <- input.list$jags.data$watch.prop[,1,]
Watch.Length[Watch.Length < 0.0001] <- 1  # day 1 and 90 are considered full coverage

n <- input.list$jags.data$n
n[is.na(n)] <- 0

BUGS.data <- list(n = n,
                  n.com = n,
                  n.sp = n,
                  n.station = max(input.list$jags.data$n.station),
                  n.year = input.list$jags.data$n.year,
                  n.obs = input.list$jags.data$n.obs,
                  periods = periods,
                  obs = input.list$jags.data$obs,
                  #Watch.Length = 0.0625,
                  u = u,
                  vs = input.list$jags.data$vs[,1,],
                  bf = input.list$jags.data$bf[,1,],
                  #day=day,
                  day = day.mat,
                  N = N,
                  N.com = N,
                  N.sp = N,
                  knot = c(-1.46,-1.26,-1.02,-0.78,
                           -0.58,-0.34,-0.10,0.10,
                           0.34,0.57,0.78,1.02,1.26,1.46),
                  n.knots = 14,
                  #begin=begin,
                  #end=end,
                  Watch.Length = input.list$jags.data$watch.prop) 

START HERE 2025-06-02
# Check data array dimensions
t.max <- BUGS.data$n.year


# Create an initial values function:
N_inits1 <- BUGS.data$n[, 1,] * 2 + 2
N_inits2 <- BUGS.data$n[, 2,] * 2 + 2 

N_inits <- N_inits1
N_inits[N_inits1 < N_inits2] <- N_inits2[N_inits1 < N_inits2]

N_inits <- rbind(N_inits,
                 matrix(data = NA, nrow = 2, ncol = dim(BUGS.data$n)[3]))

for (k in 1:dim(BUGS.data$n)[3]){
  N_inits[(BUGS.data$periods[k]+1):nrow(N_inits), k] <- NA  
}

BUGS.inits <- function() list(mean.prob = 0.5,
                              BF.Fixed = 0,
                              VS.Fixed = 0,
                              mean.prob.sp = 0.5,
                              BF.Fixed.sp = 0,
                              VS.Fixed.sp = 0,
                              mean.prob.com = 0.5,
                              BF.Fixed.com = 0,
                              VS.Fixed.com = 0,
                              mean.beta = c(0,0,0), #mean.beta = c(5,0.14,-3.5),
                              beta.sigma = c(1,1,1),#beta.sigma = c(7,7,7),
                              BF.Switch = 1,
                              VS.Switch = 1,
                              OBS.Switch = 1,
                              sigma.Obs = 1,
                              BF.Switch.sp = 1,
                              VS.Switch.sp = 1,
                              OBS.Switch.sp = 1,
                              sigma.Obs.sp = 1,
                              BF.Switch.com = 1,
                              VS.Switch.com = 1,
                              OBS.Switch.com = 1,
                              sigma.Obs.com = 1,
                              N = N_inits,
                              N.com = N_inits,
                              N.sp = N_inits,
                              #z = matrix(1,nrow=90,ncol=6),
                              beta.sp = array(data=0, dim=c(2,BUGS.data$n.year)),
                              sd.b.sp = rep(1, times = BUGS.data$n.year), #c(1,1,1,1,1,1),
                              z = matrix(1, nrow=90, ncol= BUGS.data$n.year))


MCMC.params <- list(n.iter = 85000, #200, #
                   n.thin = 50, #2, #80
                   n.burnin = 50000, #50, #
                   n.chains = 5)
# 
# This set runs quickly. Good to use this set to see if data and model are correct.
MCMC.params <- list(n.iter = 200,
                   n.thin = 2,
                   n.burnin = 50,
                   n.chains = 5)

n.samples <- ((MCMC.params$n.iter - MCMC.params$n.burnin)/MCMC.params$n.thin)*MCMC.params$n.chains

options(scipen = 999)  # to avoid having 1e+05 for 100000
out.file.name <- paste0("RData/WinBUGS_", min(input.list$jags.input.Laake$all.start.year)+1, "to",
                        max(input.list$jags.input.new$years), "_v2_min", input.list$min.dur,
                        "_", MCMC.params$n.iter, "_", Run.date, ".rds")
options(scipen = 0)

parameters <- c("lambda","OBS.RF","OBS.Switch",
                "BF.Switch","BF.Fixed","VS.Switch",
                "VS.Fixed","mean.prob", 
                "Corrected.Est","Raw.Est","z",
                "com","sp","Daily.Est","mean.beta",
                "beta.sigma","beta","beta.sp","b.sp","sd.b.sp")

if (!file.exists(out.file.name)){
  
  #Run time: 
  Start_Time<-Sys.time()
  
  BUGS_out <- bugs(data = BUGS.data,
                   inits = BUGS.inits,
                   parameters.to.save = parameters,
                   model.file="GW_Nmix_Orig.bugs",
                   n.chains = MCMC.params$n.chains,
                   n.iter = MCMC.params$n.iter, 
                   n.burnin = MCMC.params$n.burnin, 
                   n.thin = MCMC.params$n.thin,
                   debug = T,
                   bugs.directory = WinBUGS.dir,
                   DIC = FALSE)
  
  # node dimension does not match 2025-06-02
  # 
  Run_Time <- Sys.time() - Start_Time
  Ver2.results <-list(BUGS.input = WinBUGS.input,
                      BUGS.out = BUGS_out,
                      MCMCparams = MCMCparams,
                      Run_Time = Run_Time,
                      Run_Date = Run_Date,
                      Sys.info = Sys.info())
  
  saveRDS(Ver2.results, file = out.file.name)
  
}

# WinBUGS.input <- data2WinBUGS_input(data.dir = "RData/V2.1_Feb2025",
#                                     years = years,
#                                     min.dur = 30) 

# Need to convert Laake's data into WinBUGS input.


# MCMC.params <- list(n.samples = 250000,
#                     n.thin = 100,
#                     n.burnin = 200000,
#                     n.chains = 5)

# v3 does not converge well with the above MCMC setting so increasing samples
MCMC.params <- list(n.samples = 750000,
                    n.thin = 500,
                    n.burnin = 500000,
                    n.chains = 5)

# MCMC.params <- list(n.samples = 2500,
#                     n.thin = 10,
#                     n.burnin = 200,
#                     n.chains = 5)

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

###############################
# The following function uses "new" data since 2010 as well as those from Laake's 
# analysis to compute abundance since the 1967/1968 season. 
jm.out <- NoBUGS_Richards_fcn(min.dur = min.dur, 
                              ver = ver, 
                              years = years, 
                              data.dir = "RData/V2.1_Feb2025", 
                              jags.params = jags.params, 
                              MCMC.params = MCMC.params,
                              Run.date = Run.date)

# need to turn zeros into NAs when there were no second station:
data.array <- jm.out$jags.input$jags.data$n
data.array[,2,which(jm.out$jags.input$jags.data$n.station == 1)] <- NA

LOOIC.n <- compute.LOOIC(loglik.array = jm.out$jm$sims.list$log.lkhd,
                         data.array = data.array,
                         MCMC.params = jm.out$MCMC.params)

# There are some (< 0.5%) bad ones. I should look at which ones are not fitting well.

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

# v4 has one P and one K.
par.idx <- c(1:nrow(jm.out$jm$mean$S1))

if (ver == "v4"){
  mcmc_trace(jm.out$jm$samples, c("P", "K"))
  mcmc_dens(jm.out$jm$samples, c("P", "K"))
} else if (ver == "v3"){
  mcmc_trace(jm.out$jm$samples, paste0("P[", par.idx, "]"))
  mcmc_trace(jm.out$jm$samples, paste0("K[", par.idx, "]"))
  
} else {
  mcmc_trace(jm.out$jm$samples, paste0("P[", par.idx, "]"))
  mcmc_trace(jm.out$jm$samples, "K")
}

mcmc_trace(jm.out$jm$samples, paste0("S1[", par.idx, "]"))
mcmc_trace(jm.out$jm$samples, paste0("S2[", par.idx, "]"))
mcmc_trace(jm.out$jm$samples, paste0("Max[", par.idx, "]"))


# plot.trace.dens function is in Granite_Canyon_Counts_fcns.R
# ps.K <- plot.trace.dens(param = "K", 
#                         jags.out = "jm.out$jm")
# 
# ps.P <- plot.trace.dens(param = "P", 
#                         jags.out = "jm.out$jm")
# 
# 
# ps.S1 <- plot.trace.dens(param = "S1", 
#                          jags.out = "jm.out$jm")
# 
# ps.S2 <- plot.trace.dens(param = "S2", 
#                          jags.out = "jm.out$jm")
# 
# ps.Max <- plot.trace.dens(param = "Max", 
#                           jags.out = "jm.out$jm")

# These three year-specific parameters seemed to behave fine.

all.start.year <- c(jm.out$jags.input$jags.input.Laake$all.start.year,
                    jm.out$jags.input$jags.input.new$start.years)

# Create a dataframe with all years, including unsampled years.
all.years <- data.frame(start.year = seq(min(all.start.year), max(all.start.year))) %>%
  mutate(Season = paste0(start.year, "/", start.year + 1))

# Look at the annual abundance estimates:
Nhat. <- data.frame(Season = paste0(all.start.year, "/", all.start.year+1),
                    Nhat = jm.out$jm$q50$Corrected.Est,
                    LCL = jm.out$jm$q2.5$Corrected.Est,
                    UCL = jm.out$jm$q97.5$Corrected.Est) %>%
  right_join(all.years, by = "Season") %>%
  arrange(start.year) %>%
  mutate(Method = "Eguchi")

# This is for daily estimates
N.hats.day <- data.frame(Season = rep(paste0(all.start.year, "/", all.start.year+1), 
                                      each = nrow(jm.out$jm$mean$N)), #rep(Nhat.$Season, each = nrow(jm.out$jm$mean$N)),
                         Day = rep(1:nrow(jm.out$jm$mean$N), times = length(all.start.year)),
                         Mean = as.vector(jm.out$jm$mean$N),
                         LCL = as.vector(jm.out$jm$q2.5$N),
                         UCL = as.vector(jm.out$jm$q97.5$N)) 

# Daily estimates plots
p.daily.Richards <- ggplot(N.hats.day %>% group_by(Season)) + 
  geom_ribbon(aes(x = Day, ymin = LCL, ymax = UCL),
              fill = "blue", alpha = 0.5) +
  geom_path(aes(x = Day, y = Mean)) + 
  geom_point(aes(x = Day, y))
  facet_wrap(~ Season)

# These are not the best estimates because they were not updated as more data
# were collected. I should use the output from the most recent WinBUGS run for 
# the last x years.
Reported.estimates <- read.csv(file = "Data/all_estimates_2024.csv") %>%
  transmute(Season = Season,
            Nhat = Nhat,
            LCL = LCL,
            UCL = UCL,
            Method = paste0(Method, "-Reported")) %>%
  right_join(all.years, by = "Season") %>%
  arrange(start.year) %>%
  relocate(Method, .after = start.year)

WinBugs.run.date <- "2025-04-11"
WinBugs.out <- readRDS(file = paste0("RData/WinBUGS_2007to2025_v2_min", min.dur, 
                                     "_100000_",
                                     WinBugs.run.date, ".rds"))

Corrected.Est <- WinBugs.out$BUGS.out$sims.list$Corrected.Est

# We don't have raw data for 2006/2007 and 2007/2008 seasons
seasons <- c("2006/2007", "2007/2008", jm.out$jags.input$jags.input.new$seasons)

Durban.abundance.df <- data.frame(Season = seasons,
                                  Nhat = apply(Corrected.Est,
                                               FUN = mean,
                                               MARGIN = 2),
                                  # CV = apply(Corrected.Est,
                                  #            FUN = function(x) 100*sqrt(var(x))/mean(x),
                                  #            MARGIN = 2),
                                  # median = apply(Corrected.Est, 
                                  #                FUN = median, 
                                  #                MARGIN = 2),
                                  LCL = apply(Corrected.Est, 
                                              MARGIN = 2, 
                                              FUN = quantile, 0.025),
                                  UCL = apply(Corrected.Est, 
                                              MARGIN = 2, 
                                              FUN = quantile, 0.975)) %>%
  right_join(all.years, by = "Season") %>%
  arrange(start.year) %>%
  mutate(Method = "Durban")

# Create a dataframe for daily estimates:
daily.estim <- WinBugs.out$BUGS.out$sims.list$Daily.Est

# get stats:
mean.mat <- LCL.mat <- UCL.mat <- matrix(data = NA, 
                                         nrow = dim(daily.estim)[2], 
                                         ncol = dim(daily.estim)[3])

for (k1 in 1:dim(daily.estim)[2]){
  for (k2 in 1:dim(daily.estim)[3]){
    mean.mat[k1, k2] <- mean(daily.estim[,k1,k2])
    LCL.mat[k1, k2] <- quantile(daily.estim[,k1,k2], 0.025)
    UCL.mat[k1, k2] <- quantile(daily.estim[,k1,k2], 0.975)
  }
  
}

N.hats.day.Durban <- data.frame(Season = rep(seasons, each = dim(daily.estim)[2]),
                                Day = rep(1:dim(daily.estim)[2], length(seasons)),
                                Mean = as.vector(mean.mat),
                                LCL = as.vector(LCL.mat),
                                UCL = as.vector(UCL.mat))

# Daily estimates plots
p.daily.Durban <- ggplot(N.hats.day.Durban %>% group_by(Season)) + 
  geom_ribbon(aes(x = Day, ymin = LCL, ymax = UCL),
              fill = "blue", alpha = 0.5) +
  geom_path(aes(x = Day, y = Mean)) + 
  facet_wrap(~ Season)

# Include non-survey years - no estimates for 2007/2008 because I don't have
# raw data for that year. Only the WinBUGS inputs. 
Laake.abundance.new <- read.csv(file = "Data/all_estimates_Laake_2025.csv") %>%
  mutate(LCL = CL.low,
         UCL = CL.high) %>%
  select(c(Season, Nhat, LCL, UCL)) %>%
  right_join(all.years, by = "Season") %>%
  arrange(start.year) %>%
  mutate(Method = "Laake")

#Laake.output <- read_rds(file = "RData/Laake_abundance_estimates_2024.rds")

# In reported estimates, there are two 2006/2007.
Reported.estimates %>%
  na.omit() %>%
  select(Season) %>% 
  unique() -> sampled.seasons 

Laake.abundance.new %>%
  rbind(Durban.abundance.df) %>%
  rbind(Nhat.) -> all.estimates
#  rbind(spline.Nhat) 
#  rbind(Reported.estimates %>% na.omit()) -> all.estimates

p.Nhats <- ggplot(all.estimates) +
  geom_point(aes(x = start.year, y = Nhat,
                 color = Method),
             alpha = 0.5) +
  geom_errorbar(aes(x = start.year, ymin = LCL, ymax = UCL,
                    color = Method)) +
  ylim(5000, 35000)

# ggsave(plot = p.Nhats,
#        filename = paste0("figures/Nhats_", ver, "_", min.dur, "min.png"),
#        device = "png",
#        dpi = 600)

Nhat. %>% 
  select(Season, start.year, Nhat, LCL, UCL) %>%
  rename(Nhat.Eguchi = Nhat,
         LCL.Eguchi = LCL,
         UCL.Eguchi = UCL) %>%
  cbind(Laake.abundance.new %>%
          select(Nhat, LCL, UCL) %>%
          rename(Nhat.Laake = Nhat,
                 LCL.Laake = LCL,
                 UCL.Laake = UCL)) %>%
  cbind(Durban.abundance.df %>%
          select(Nhat, LCL, UCL) %>%
          rename(Nhat.Durban = Nhat,
                 LCL.Durban = LCL,
                 UCL.Durban = UCL)) %>%
  mutate(d.Laake.Eguchi = Nhat.Laake - Nhat.Eguchi,
         d.Durban.Eguchi = Nhat.Durban - Nhat.Eguchi) -> Nhat.all.wide

# Check convergence
high.Rhat <- function(x){
  return(data.frame(idx = which(x > 1.1),
                    start.year = all.start.year[which(x > 1.1)],
                    Rhat = x[which(x > 1.1)]))
}

high.Rhat.Max <- high.Rhat(jm.out$jm$Rhat$Max)

high.Rhat.K <- high.Rhat(jm.out$jm$Rhat$K)
high.Rhat.S1 <- high.Rhat(jm.out$jm$Rhat$S1)
high.Rhat.S2 <- high.Rhat(jm.out$jm$Rhat$S2)
high.Rhat.P <- high.Rhat(jm.out$jm$Rhat$P)


