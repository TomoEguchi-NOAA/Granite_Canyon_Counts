# Jags_Richards_AllData.R
# 
# Combines Laake data and more recent data and runs
# model_Richards_pois_bino_vX, where X is version number. See below.  
# Some diagnostics are conducted.

# M1 <- (1 + (2 * exp(K) - 1) * exp((1/S1) * (P - d))) ^ (-1/exp(K))
# M2 <- (1 + (2 * exp(K) - 1) * exp((1/S2) * (P - d))) ^ (-1/exp(K))
# N <- min + (max - min) * (M1 * M2)
#
# d is the number of days from the beginning of nesting season
# S1 < 0 and S2 > 0 define the "fatness" of the function
# K > 0 defines the "flatness" at the peak of the function
# P defines where the peak is relative to the range of d; min(d) < P < max(d)
# min is "the basal level of the number of whales outside the migration season"
# max > min. min == 0.

# Model versions - depending on which parameters are constant/year-specific
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
options(mc.cores = 5)
save.plot <- F

WinBUGS.Run.Date <- "2025-04-11"
#WinBUGS.Run.Date <- "2025-06-06"

Run.date <- "2025-06-24" #Sys.Date() #"2025-04-21" #"2025-04-17" #

# Minimum length of observation periods in minutes
min.dur <- 60 #10 #85 #

ver <- "v5a" #  "v3" #"v5" # "v4" # 

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

# MCMC.params <- list(n.samples = 200000,
#                     n.thin = 100,
#                     n.burnin = 150000,
#                     n.chains = 5)

# MCMC.params <- list(n.samples = 250000,
#                     n.thin = 100,
#                     n.burnin = 200000,
#                     n.chains = 5)
# 
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
# MCMC.params <- list(n.samples = 10000,
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
                 "K.alpha", "K.beta",
                 #"beta.1",
                 #"N.alpha", "N.obs",
                 "log.lkhd")


# The following function uses "new" data since 2010 as well as those from Laake's 
# analysis to compute abundance since the 1967/1968 season. There were two seasons
# where the survey continued beyond the 90th day. So, max.day needs to be increased
# from 90. I used 100. 
jm.out <- NoBUGS_Richards_fcn(min.dur = min.dur, 
                              ver = ver, 
                              years = years, 
                              data.dir = "RData/V2.1_Feb2025", 
                              jags.params = jags.params, 
                              MCMC.params = MCMC.params,
                              Run.date = Run.date,
                              obs.n.min = 10,
                              max.day = 100)

if (!jm.out$new.run){
  jags.params <- jm.out$jags.params
  MCMC.params <- jm.out$MCMC.params
}

# need to turn zeros into NAs when there were no second station:
data.array <- jm.out$jags.input$jags.data$n
data.array[,2,which(jm.out$jags.input$jags.data$n.station == 1)] <- NA
data.array[,2,which(jm.out$jags.input$jags.data$n.station == 1)] <- NA

LOOIC.n <- compute.LOOIC(loglik.array = jm.out$jm$sims.list$log.lkhd,
                         data.array = data.array,
                         MCMC.params = MCMC.params)

# There are some (< 0.5%) bad ones. I should look at which ones are not fitting well.

# Compute new Rhat (rank-normalized Rhat as per Vehtari et al. 2021)
params <- "^VS\\.Fixed|^BF\\.Fixed|^Max|^S1|^S2|^P1|^P2|^K|^OBS\\.RF\\["
new.Rhat <- rank.normalized.R.hat(jm.out$jm$samples, params)
#max.Rhat <- lapply(.out[[k]]$jm$Rhat, FUN = max, na.rm = T) %>%
#  unlist()
max.new.Rhat.big <- new.Rhat[which(new.Rhat > 1.01)]

# Look at Rhat statistics
# max.Rhat <- lapply(jm.out$jm$Rhat, FUN = max, na.rm = T) %>%
#   unlist()
# max.Rhat.big <- max.Rhat[which(max.Rhat > 1.1)]

if (grepl("a", ver)){
  # mcmc_dens(jm.out$jm$samples, c("S1.alpha", "S1.beta",
  #                                "S2.alpha", "S2.beta"))
                                 #"P.alpha", "P.beta",
                                 #"K.alpha", "K.beta"))
  # P.alpha and P.beta seem to be not behaving well - the right tails are not 
  # captured. 
  p.trace.hyper.params <- mcmc_trace(jm.out$jm$samples, c("S1.alpha", "S1.beta",
                                                          "S2.alpha", "S2.beta",
                                                          #"P.alpha", "P.beta",
                                                          "K.alpha", "K.beta"))
  
}

p.trace.Fixed.params <- mcmc_trace(jm.out$jm$samples, c("BF.Fixed", "VS.Fixed"))

all.start.year <- c(jm.out$jags.input$jags.input.Laake$all.start.year,
                    jm.out$jags.input$jags.input.new$start.years)

# v4 has one P and one K. S1 is always year specific
par.idx <- c(1:jm.out$jags.input$jags.data$n.year)

# K trace plots
if (grepl("v2(?!\\d[a-zA-Z])", ver, perl = TRUE) | 
    grepl("v5", ver) | 
    grepl("v15", ver) | 
    grepl("v16", ver) |
    grepl("v17", ver) |
    grepl("v18", ver) |
    grepl("v19", ver) |
    grepl("v20", ver)){
  p.trace.K <- mcmc_trace(jm.out$jm$samples, c("K"))
} else if (grepl("v1(?!\\d[a-zA-Z])", ver, perl = TRUE) | 
           grepl("v3", ver) | 
           grepl("v12", ver) | 
           grepl("v13", ver) | 
           grepl("v21", ver) | 
           grepl("v22", ver) | 
           grepl("v23", ver) | 
           grepl("v24", ver)){
  p.trace.K <- mcmc_trace(jm.out$jm$samples, paste0("K[", par.idx, "]"))
  high.Rhat.K <- high.Rhat(new.Rhat[grep("K", names(new.Rhat))])  
} else if (grepl("v9", ver) |
           grepl("v10", ver) |
           grepl("v11", ver) |
           grepl("v14", ver) |
           grepl("v25", ver) |
           grepl("v26", ver) |
           grepl("v27", ver) |
           grepl("v28", ver)){
  p.trace.K1 <- mcmc_trace(jm.out$jm$samples, paste0("K1[", par.idx, "]"))
  p.trace.K2 <- mcmc_trace(jm.out$jm$samples, paste0("K2[", par.idx, "]"))
  high.Rhat.K1 <- high.Rhat(new.Rhat[grep("K1", names(new.Rhat))])  
  high.Rhat.K2 <- high.Rhat(new.Rhat[grep("K2", names(new.Rhat))])    
}

# P1 and P2 trace plots
if (grepl("v1(?!\\d[a-zA-Z])", ver, perl = TRUE) | 
    grepl("v2(?!\\d[a-zA-Z])", ver, perl = TRUE) | 
    grepl("v10", ver) | 
    grepl("v11", ver) | 
    grepl("v13", ver) | 
    grepl("v15", ver) | 
    grepl("v17", ver) | 
    grepl("v18", ver) | 
    grepl("v21", ver)| 
    grepl("v22", ver)| 
    grepl("v25", ver)| 
    grepl("v26", ver)){
  p.trace.P1 <- mcmc_trace(jm.out$jm$samples, paste0("P1[", par.idx, "]"))
  p.trace.P2 <- mcmc_trace(jm.out$jm$samples, paste0("P2[", par.idx, "]"))
  high.Rhat.P1 <- high.Rhat(new.Rhat[grep("P1", names(new.Rhat))])
  high.Rhat.P2 <- high.Rhat(new.Rhat[grep("P2", names(new.Rhat))])
}

# P trace plots
if (grepl("v3", ver) | 
    grepl("v5", ver) | 
    grepl("v9", ver) | 
    grepl("v12", ver)| 
    grepl("v14", ver)| 
    grepl("v16", ver)| 
    grepl("v19", ver)| 
    grepl("v20", ver)| 
    grepl("v16", ver)| 
    grepl("v23", ver)| 
    grepl("v24", ver)| 
    grepl("v27", ver)| 
    grepl("v28", ver)){
  p.trace.P <- mcmc_trace(jm.out$jm$samples, paste0("P[", par.idx, "]"))
  high.Rhat.P <- high.Rhat(new.Rhat[grep("P", names(new.Rhat))])
  
}

# Max trace plots
p.trace.Max <- mcmc_trace(jm.out$jm$samples, paste0("Max[", par.idx, "]"))
high.Rhat.Max <- high.Rhat(new.Rhat[grep("Max", names(new.Rhat))])

# S1 and S2 trace plots
if (grepl("v1(?!\\d[a-zA-Z])", ver, perl = TRUE) | 
    grepl("v2(?!\\d[a-zA-Z])", ver, perl = TRUE) | 
    grepl("v3", ver) | 
    grepl("v5", ver) | 
    grepl("v9", ver) | 
    grepl("v10", ver)){
  p.trace.S1 <- mcmc_trace(jm.out$jm$samples, paste0("S1[", par.idx, "]"))
  p.trace.S2 <- mcmc_trace(jm.out$jm$samples, paste0("S2[", par.idx, "]"))
  high.Rhat.S1 <- high.Rhat(new.Rhat[grep("S1\\[", names(new.Rhat))])
  high.Rhat.S2 <- high.Rhat(new.Rhat[grep("S2\\[", names(new.Rhat))])
  
} else if (grepl("v11", ver) | 
           grepl("v12", ver) | 
           grepl("v13", ver) | 
           grepl("v14", ver) | 
           grepl("v15", ver) | 
           grepl("v16", ver)){
  p.trace.S1.S2 <- mcmc_trace(jm.out$jm$samples, c("S1", "S2"))
  high.Rhat.S1 <- high.Rhat(new.Rhat[grep("S1", names(new.Rhat))])
  high.Rhat.S2 <- high.Rhat(new.Rhat[grep("S2", names(new.Rhat))])
  
} else if (grepl("v17", ver) |
           grepl("v19", ver) |
           grepl("v21", ver) |
           grepl("v23", ver) |
           grepl("v25", ver) |
           grepl("v27", ver)){
  p.trace.S1 <- mcmc_trace(jm.out$jm$samples, paste0("S1[", par.idx, "]"))
  p.trace.S2 <- mcmc_trace(jm.out$jm$samples, c("S2"))
} else if (grepl("v18", ver) |
           grepl("v20", ver) |
           grepl("v22", ver) |
           grepl("v24", ver) |
           grepl("v26", ver) |
           grepl("v28", ver)){
  p.trace.S2 <- mcmc_trace(jm.out$jm$samples, paste0("S2[", par.idx, "]"))
  p.trace.S1 <- mcmc_trace(jm.out$jm$samples, c("S1"))
}

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
  mutate(Method = paste0("Eguchi ", ver))

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
  #geom_point(aes(x = Day, y = Mean)) +
  facet_wrap(~ Season)

# These are not the best estimates because they were not updated as more data
# were collected. I should use the output from the most recent WinBUGS run for 
# the last x years.
#Reported.estimates <- read.csv(file = "Data/all_estimates_2024.csv") %>%
Reported.estimates <- read.csv(file = "Data/Nhats_2025.csv") %>%  
  transmute(Season = Season,
            Nhat = Nhat,
            LCL = LCL,
            UCL = UCL,
            Method = paste0(Method, "-Reported")) %>%
  right_join(all.years, by = "Season") %>%
  arrange(start.year) %>%
  relocate(Method, .after = start.year)

#WinBugs.run.date <- "2025-04-11"
WinBugs.out <- readRDS(file = paste0("RData/WinBUGS_2007to2025_v2_min", min.dur,
                                     "_100000_",
                                     WinBUGS.Run.Date, ".rds"))

# Compute rank-normalized Rhat for WinBUGS output
BUGS.params <- "^lambda|^beta|^OBS.RF"
BUGS.samples <- WinBugs.out$BUGS.out$sims.array
BUGS.col.names <- grep(BUGS.params, dimnames(BUGS.samples)[[3]], 
                  value = T, perl = T)
subset.BUGS.samples <- BUGS.samples[,,BUGS.col.names]
subset.BUGS.mcmc.array <- as_draws_array(subset.BUGS.samples, .nchains = 5)
BUGS.Rhat <- apply(subset.BUGS.mcmc.array, MARGIN = 3, FUN = rhat)

#max.Rhat <- lapply(.out[[k]]$jm$Rhat, FUN = max, na.rm = T) %>%
#  unlist()
BUGS.Rhat.big <- BUGS.Rhat[which(BUGS.Rhat > 1.01)]

# WinBugs.out <- readRDS(file = paste0("RData/WinBUGS_1968to2025_v2_min", min.dur, 
#                                      "_85000_",
#                                      WinBUGS.Run.Date, ".rds"))

Corrected.Est <- WinBugs.out$BUGS.out$sims.list$Corrected.Est

# We don't have raw data for 2006/2007 and 2007/2008 seasons
seasons <- c("2006/2007", "2007/2008", jm.out$jags.input$jags.input.new$seasons)

all.season <- paste0(all.start.year, "/", all.start.year+1)
Durban.abundance.df <- data.frame(Season = WinBugs.out$BUGS.input$seasons,
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

N.hats.day.Durban <- data.frame(Season = rep(WinBugs.out$BUGS.input$seasons, 
                                             each = dim(daily.estim)[2]),
                                Day = rep(1:dim(daily.estim)[2], 
                                          length(WinBugs.out$BUGS.input$seasons)),
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

# Reported estimates are identical to the reanalysis so remove. 
Laake.abundance.new %>%
  rbind(Durban.abundance.df) %>%
  rbind(Nhat.) -> all.estimates
#  rbind(spline.Nhat) 
  #rbind(Reported.estimates %>% na.omit()) -> all.estimates

p.Nhats <- ggplot(all.estimates) +
  geom_point(aes(x = start.year, y = Nhat,
                 color = Method),
             alpha = 0.5) +
  geom_errorbar(aes(x = start.year, ymin = LCL, ymax = UCL,
                    color = Method)) +
  ylim(2000, 40000) +
  theme(legend.position = "top")

if (save.plot)
  ggsave(plot = p.Nhats,
         filename = paste0("figures/Nhats_", ver, "_", min.dur, "min.png"),
         device = "png",
         dpi = 600)

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


# Compare how daily sums among years
obsd.periods.primary <- jm.out$jags.input$jags.data$periods[,1]
watch.prop.primary <- jm.out$jags.input$jags.data$watch.prop[,1,]
obsd.effort.primary <- rbind(rep(0, times = dim(watch.prop.primary)[2]), 
                             watch.prop.primary, 
                             rep(0, times = dim(watch.prop.primary)[2]))

obsd.n.primary <- jm.out$jags.input$jags.data$n[,1,]
obsd.day.primary <- jm.out$jags.input$jags.data$day[,1,]
obsd.n.prop <- obsd.n.primary[,] * obsd.effort.primary

obsd.n.df <- data.frame(Season = rep(all.season, each = dim(obsd.n.prop)[1]),
                        obsd.n = as.vector(obsd.n.prop),
                        day = as.vector(obsd.day.primary),
                        effort = as.vector(obsd.effort.primary)) %>%
  na.omit()

# ggplot(obsd.n.df) +
#   geom_point(aes(x = day, y = obsd.n)) +
#   facet_wrap(~ Season)
# 
# ggplot(obsd.n.df) +
#   geom_point(aes(x = day, y = effort)) +
#   facet_wrap(~ Season)

p.RF.Obs <- plot.trace.dens(jm.out$jm, "OBS.RF")

# Simple comparison between observed counts per hour vs. estimated abundance
# obsd.n.prop.sum <- data.frame(Season = all.season,
#                               n.sum = colSums(obsd.n.prop, na.rm = T))
# 
# Nhat.all.wide %>% 
#   left_join(obsd.n.prop.sum, by = "Season") -> Nhat.all.wide 
# 
# ggplot(Nhat.all.wide) +
#   geom_point(aes(y = Nhat.Laake, x = n.sum), color = "blue") +
#   geom_point(aes(y = Nhat.Eguchi, x = n.sum), color = "red")
# 
# Nhat.all.wide %>%
#   filter(n.sum < 1000) -> Nhat.all.wide.1000
# 
# ggplot(Nhat.all.wide.1000) +
#   geom_point(aes(y = Nhat.Laake, x = n.sum), color = "blue") +
#   geom_point(aes(y = Nhat.Eguchi, x = n.sum), color = "red")
