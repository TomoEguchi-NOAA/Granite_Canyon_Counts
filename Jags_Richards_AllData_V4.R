# Jags_Richards_AllData_V4
# 
# Combines Laake data and more recent data and runs
# model_Richards_pois_bino_v4.txt. 

# The v4 model has constant P and K. It seems that these two parameters are
# a bit of problem. In v5, P is year-specific and K is constant over years.
# 
 

rm(list = ls())

#library(ERAnalysis)
library(tidyverse)
library(ggplot2)
library(loo)
library(bayesplot)

source("Granite_Canyon_Counts_fcns.R")
source("LaakeData2Jags.R")

Run.date <- Sys.Date()
#Run.date <- "2024-11-01"
# I use the "Laake" model because the number of periods per year
# can be different. This should be the same as non-Laake version
# as of 2024-07-09

model.name <- "Richards_pois_bino_v4"
jags.model <- paste0("models/model_", model.name, ".txt")

out.file.name <- paste0("RData/JAGS_", model.name, "_AllYears_",
                        Run.date, ".rds")

# This comes from LaakeData2Jags.R above
Laake.data <- jags.data

# More recent data from the output of JAGS Richards Ver1.Rmd:
run.date <- "2024-11-01"
file.name <- paste0("RData/JAGS_Richards_pois_bino_v4_10yr_v2_", 
                           run.date, ".rds")

jm.out <- readRDS(file.name)

.data <- jm.out$jags.data
.start.year <- c(2007, 2009, 2010, 2014, 2015, 2019, 2021, 2022, 2023)

all.start.year <- c(Laake.start.year, .start.year)

# Both datasets contain 2006. Take it out from the recent one.
bf <- vs <- all.n <- all.obs <- array(dim = c(max(dim(Laake.data$n)[1], dim(.data$n)[1]), 
                                              2, (dim(Laake.data$n)[3] + dim(.data$n)[3]-1)))

day <- watch.prop <- array(dim = c(max(dim(Laake.data$n)[1], dim(.data$n)[1]), 
                                   2, (dim(Laake.data$n)[3] + dim(.data$n)[3]-1)))

c2 <- 2
c1 <- c <- 1
for (k in 1:length(all.start.year)){
  if (all.start.year[k] < 2007){
    all.n[1:dim(Laake.data$n)[1], 1, c] <- Laake.data$n[, 1, c1]
    all.n[1:dim(Laake.data$n)[1], 2, c] <- Laake.data$n[, 2, c1]
    
    all.obs[1:dim(Laake.data$obs)[1], 1, c] <- Laake.data$obs[, 1, c1]
    all.obs[1:dim(Laake.data$obs)[1], 2, c] <- Laake.data$obs[, 2, c1]
    
    watch.prop[1:dim(Laake.data$watch.prop)[1], 1, c] <- Laake.data$watch.prop[, 1, c1]
    watch.prop[1:dim(Laake.data$watch.prop)[1], 2, c] <- Laake.data$watch.prop[, 2, c1]
    
    day[1:dim(Laake.data$day)[1], 1, c] <- Laake.data$day[, 1, c1]
    day[1:dim(Laake.data$day)[1], 2, c] <- Laake.data$day[, 2, c1]
    
    bf[1:dim(Laake.data$bf)[1], 1, c] <- Laake.data$bf[, 1, c1]
    bf[1:dim(Laake.data$bf)[1], 2, c] <- Laake.data$bf[, 2, c1]
    
    vs[1:dim(Laake.data$vs)[1], 1, c] <- Laake.data$vs[, 1, c1]
    vs[1:dim(Laake.data$vs)[1], 2, c] <- Laake.data$vs[, 2, c1]
    
    c <- c + 1
    c1 <- c1 + 1
  } else {
    all.n[1:1:dim(.data$n)[1], 1, c] <- .data$n[, 1, c2]
    all.n[1:1:dim(.data$n)[1], 2, c] <- .data$n[, 2, c2]
    all.obs[1:1:dim(.data$obs)[1], 1, c] <- .data$obs[, 1, c2]
    all.obs[1:1:dim(.data$obs)[1], 2, c] <- .data$obs[, 2, c2]
    watch.prop[1:1:dim(.data$watch.prop)[1], 1, c] <- .data$watch.prop[, 1, c2]
    watch.prop[1:1:dim(.data$watch.prop)[1], 2, c] <- .data$watch.prop[, 2, c2]
    day[1:1:dim(.data$day)[1], 1, c] <- .data$day[, 1, c2]
    day[1:1:dim(.data$day)[1], 2, c] <- .data$day[, 2, c2]
    bf[1:1:dim(.data$bf)[1], 1, c] <- .data$bf[, 1, c2]
    bf[1:1:dim(.data$bf)[1], 2, c] <- .data$bf[, 2, c2]
    vs[1:1:dim(.data$vs)[1], 1, c] <- .data$vs[, 1, c2]
    vs[1:1:dim(.data$vs)[1], 2, c] <- .data$vs[, 2, c2]
    c <- c + 1
    c2 <- c2 + 1
  }
  
}

n.station <- c(Laake.data$n.station, .data$n.station[2:length(.data$n.station)])
n.year <- length(all.start.year)
n.obs <- length(unique(as.vector(all.obs))) 

# tmp <- matrix(nrow = length(.data$periods)-1, ncol = 2)
# tmp[,1] <- .data$periods[2:length(.data$periods)]
# tmp[.data$n.station[2:9] == 2, 2] <- tmp[.data$n.station[2:9] == 2, 1]
periods <- rbind(Laake.data$periods, .data$periods[2:length(.data$n.station),])

# vs <- cbind(Laake.data$vs[,1,], 
#             rbind(.data$vs[, 2:ncol(.data$vs)], 
#                   matrix(nrow = dim(Laake.data$vs)[1] - nrow(.data$vs), 
#                          ncol = ncol(.data$vs)-1)))
# 
# bf <- cbind(Laake.data$bf[,1,], 
#              rbind(.data$bf[, 2:ncol(.data$bf)], 
#                    matrix(nrow = dim(Laake.data$bf)[1] - nrow(.data$bf), 
#                           ncol = ncol(.data$bf)-1)))

jags.data <- list(n = all.n,
                  n.station = n.station,
                  n.year = n.year,
                  n.obs = n.obs,
                  periods = periods,
                  obs = all.obs,
                  vs = vs,
                  bf = bf,
                  watch.prop = watch.prop,
                  day = day,
                  n.days = 94)   # 94 comes from Laake's data. Their maximum observation period was 94 days

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
                 jags.data = jags.data,
                 start.year = all.start.year,
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
data.array <- jm.out$jags.data$n
data.array[,2,which(jm.out$jags.data$n.station == 1)] <- NA

LOOIC.n <- compute.LOOIC(loglik.array = jm.out$jm$sims.list$log.lkhd,
                         data.array = data.array,
                         MCMC.params = jm.out$MCMC.params)

# There are some (< 2%) bad ones. I should look at which ones are not fitting well.
 

# Horwitz-Thompson estimates:
Nhats.HT <- jm.out$jags.data$n[,1,]/jm.out$jags.data$watch.prop[,1,]
day.idx <- jm.out$jags.data$day[,1,]

Nhats.HT.day <- matrix(nrow = 94, ncol = dim(Nhats.HT)[2])
for (k1 in 1:dim(Nhats.HT)[2]){
  for (k2 in 1:94){
    Nhats.HT.day[k2, k1] <- sum(Nhats.HT[day.idx[,k1] == k2, k1], na.rm = T)
  }
}
  
Nhats.HT.all <- data.frame(Season = rep(paste0(all.start.year, "/", all.start.year+1), 
                                        each = nrow(Nhats.HT.day)),
                           Day = rep(seq(1, 94), times = ncol(Nhats.HT.day)),
                           Nhat = as.vector(Nhats.HT.day))


# Look at Rhat statistics
max.Rhat <- lapply(jm.out$jm$Rhat, FUN = max, na.rm = T) %>%
  unlist()
max.Rhat.big <- max.Rhat[which(max.Rhat > 1.1)]

bayesplot::mcmc_dens(jm.out$jm$samples, c("S1.alpha", "S1.beta",
                                          "S2.alpha", "S2.beta",
                                          "P.alpha", "P.beta",
                                          "K.alpha", "K.beta"))
# P.alpha and P.beta seem to be not behaving well - the right tails are not 
# captured. 
bayesplot::mcmc_trace(jm.out$jm$samples, c("S1.alpha", "S1.beta",
                                          "S2.alpha", "S2.beta",
                                          "P.alpha", "P.beta",
                                          "K.alpha", "K.beta"))

bayesplot::mcmc_dens(jm.out$jm$samples, c("BF.Fixed", "VS.Fixed"))

# v4 has one P and one K.
par.idx <- c(1:nrow(jm.out$jm$mean$S1))

mcmc_trace(jm.out$jm$samples, c("P", "K"))
mcmc_dens(jm.out$jm$samples, c("P", "K"))

# mcmc_trace(jm.out$jm$samples, paste0("P[", par.idx, "]"))
# mcmc_trace(jm.out$jm$samples, paste0("K[", par.idx, "]"))
bayesplot::mcmc_trace(jm.out$jm$samples, paste0("S1[", par.idx, "]"))
bayesplot::mcmc_trace(jm.out$jm$samples, paste0("S2[", par.idx, "]"))
bayesplot::mcmc_trace(jm.out$jm$samples, paste0("Max[", par.idx, "]"))
# These three year-specific parameters seemed to behave fine.

# Create a dataframe with all years, including unsampled years.
all.years <- data.frame(year = seq(min(all.start.year), max(all.start.year))) %>%
  mutate(Season = paste0(year, "/", year + 1))

# Look at the annual abundance estimates:
Nhat. <- data.frame(Season = paste0(all.start.year, "/", all.start.year+1),
                    Nhat = jm.out$jm$q50$Corrected.Est,
                    LCL = jm.out$jm$q2.5$Corrected.Est,
                    UCL = jm.out$jm$q97.5$Corrected.Est) %>%
  right_join(all.years, by = "Season") %>%
  arrange(year) %>%
  mutate(Method = "Richards")

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
  geom_point(data = Nhats.HT.all,
             aes(x = Day, y = Nhat),
             alpha = 0.3) +
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
  arrange(year) %>%
  relocate(Method, .after = year)
  

WinBugs.out <- readRDS(file = "RData/WinBUGS_10yr_v2_min85.rds")
Corrected.Est <- WinBugs.out$BUGS_out$sims.list$Corrected.Est
seasons <- c("2006/2007", "2007/2008", "2009/2010", "2010/2011", 
             "2014/2015", "2015/2016", "2019/2020", "2021/2022",
             "2022/2023", "2023/2024")

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
  arrange(year) %>%
  mutate(Method = "Durban")

# Create a dataframe for daily estimates:
daily.estim <- WinBugs.out$BUGS_out$sims.list$Daily.Est

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
  geom_point(data = Nhats.HT.all[Nhats.HT.all$Season %in% seasons,],
             aes(x = Day, y = Nhat),
             alpha = 0.3) + 
  facet_wrap(~ Season)

# Include non-survey years - no estimates for 2007/2008 because I don't have
# raw data for that year. Only the WinBUGS inputs. 
Laake.abundance.new <- read.csv(file = "Data/all_estimates_Laake_2024.csv") %>%
  mutate(LCL = CL.low,
         UCL = CL.high) %>%
  select(c(Season, Nhat, LCL, UCL)) %>%
  right_join(all.years, by = "Season") %>%
  arrange(year) %>%
  mutate(Method = "Laake")

#Laake.output <- read_rds(file = "RData/Laake_abundance_estimates_2024.rds")

# In reported estimates, there are two 2006/2007.
Reported.estimates %>%
  na.omit() %>%
  select(Season) %>% 
  unique() -> sampled.seasons 

# read in spline results - Not very good so remove
# JAGS spline Ver1.Rmd
# spline.out <- read_rds("RData/JAGS_Spline_v2_results_All_Data_2024-07-09.rds")
# 
# spline.Nhat <- data.frame(Season = sampled.seasons$Season,
#                           Nhat = spline.out$jm$q50$Corrected.Est,
#                           LCL = spline.out$jm$q2.5$Corrected.Est,
#                           UCL = spline.out$jm$q97.5$Corrected.Est) %>%
#   right_join(all.years, by = "Season") %>%
#   arrange(year) %>%
#   mutate(Method = "Bayesian Spline")
# 
# daily.estim.spline <- exp(spline.out$jm$sims.list$sp)

# get stats:
# mean.mat.spline <- LCL.mat.spline <- UCL.mat.spline <- matrix(data = NA, 
#                                          nrow = dim(daily.estim.spline)[2], 
#                                          ncol = dim(daily.estim.spline)[3])
# 
# for (k1 in 1:dim(daily.estim.spline)[2]){
#   for (k2 in 1:dim(daily.estim.spline)[3]){
#     mean.mat.spline[k1, k2] <- mean(daily.estim.spline[,k1,k2])
#     LCL.mat.spline[k1, k2] <- quantile(daily.estim.spline[,k1,k2], 0.025)
#     UCL.mat.spline[k1, k2] <- quantile(daily.estim.spline[,k1,k2], 0.975)
#   }
#   
# }
# 
# N.hats.day.spline <- data.frame(Season = rep(paste0(all.start.year, "/", all.start.year+1), 
#                                                        each = nrow(spline.out$jm$mean$sp)),
#                                 Day = rep(1:dim(daily.estim.spline)[2], dim(daily.estim.spline)[3]),
#                                 Mean = as.vector(mean.mat.spline),
#                                 LCL = as.vector(LCL.mat.spline),
#                                 UCL = as.vector(UCL.mat.spline))
# 
# # Daily estimates plots
# p.daily.spline <- ggplot(N.hats.day.spline %>% group_by(Season)) + 
#   geom_ribbon(aes(x = Day, ymin = LCL, ymax = UCL),
#               fill = "blue", alpha = 0.5) +
#   geom_path(aes(x = Day, y = Mean)) + 
#   geom_point(data = Nhats.HT.all,
#              aes(x = Day, y = Nhat),
#              alpha = 0.3) + 
#   facet_wrap(~ Season)

# Laake.abundance.new %>%
#   rbind(Durban.abundance.df) %>%
#   rbind(Nhat.) %>%
#   rbind(spline.Nhat) %>%
#   rbind(Reported.estimates) -> all.estimates
# Reported.estimates %>%
#   filter(Method == "Laake") %>%
#   rbind(Durban.abundance.df) -> previous.estimates

# ggplot(all.estimates) +
#   geom_point(aes(x = year, y = Nhat,
#                  color = Method)) +
#   geom_errorbar(aes(x = year, ymin = LCL, ymax = UCL,
#                     color = Method)) +
#   ylim(0, 50000)

# Bayesian spline is the worst, so remove and re-plot
Laake.abundance.new %>%
  rbind(Durban.abundance.df) %>%
  rbind(Nhat.) -> all.estimates
#  rbind(spline.Nhat) 
#  rbind(Reported.estimates %>% na.omit()) -> all.estimates

p.Nhats <- ggplot(all.estimates) +
  geom_point(aes(x = year, y = Nhat,
                 color = Method),
             alpha = 0.5) +
  geom_errorbar(aes(x = year, ymin = LCL, ymax = UCL,
                    color = Method)) +
  ylim(0, 45000)

# Precision of the estimates using Richards function is a lot better than other 
# methods. They are a bit lower (negatively biased) than other two. I think the
# problem is detection probabilities. mean.prob gets very small (mean = 0.04) 
# even though the prior is beta(0.95, 0.05). 

# Changed the prior to unif(0.9, 1.0) - 2024-10-21

# This problem has been fixed as of November 1, 2024. The problem was how I set
# up the proportion of watch period. It was set at hours observed over maximum
# possible. But, this should have been over 24 hrs (or one day), which was 
# already computed when data were extracted. By dividing that (times 60 min times
# 24 hrs) by 540 minutes made those numbers bigger by 2.667 times (24*60/540), 
# resulting in a smaller estimated abundance. 

# Although the estimates are a lot better now they are a bit off, especially in 
# some years. These need to be looked at carefully. 2024-11-13



