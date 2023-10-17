# Jags_Richards_AllData_V3
# 
# Combines Laake data and more recent data and runs
# model_Richards_pois_bino_v3.txt. 
# 
 

rm(list = ls())

#library(ERAnalysis)
library(tidyverse)
library(ggplot2)
library(loo)
library(bayesplot)

source("Granite_Canyon_Counts_fcns.R")

Run.date <- Sys.Date()
#Run.date <- "2023-10-13"
out.file.name <- paste0("RData/JAGS_Richards_pois_bino_v3_AllYears_",
                        Run.date, ".rds")

# I use the "Laake" model because the number of periods per year
# can be different.
jags.model <- "models/model_Richards_pois_bino_v3_Laake.txt"

# Bring in the output from the most recent Jags run:
run.date.Laake <- "2023-10-06" #Sys.Date() # # "2023-08-11"
Laake.file.name <- paste0("RData/JAGS_Richards_v3_Laake_", 
                        run.date.Laake, ".rds")

Laake.jm.out <- readRDS(Laake.file.name)

# effort
Laake_PrimaryEffort <- read.csv(file = "Data/Laake_PrimaryEffort.csv") %>%
  filter(Use)
Laake.start.year <- unique(Laake_PrimaryEffort$Start.year)
Laake.data <- Laake.jm.out$jags.data

# More recent data:
run.date <- "2023-10-13"
file.name <- paste0("RData/JAGS_Richards_pois_bino_v3_9yr_v2_", 
                           run.date, ".rds")

jm.out <- readRDS(file.name)

.data <- jm.out$jags.data
.start.year <- c(2007, 2009, 2010, 2014, 2015, 2019, 2021, 2022)

# seasons <- c("2006/2007", "2007/2008", "2009/2010", "2010/2011",
#              "2014/2015", "2015/2016", "2019/2020", "2021/2022",
#              "2022/2023")

all.start.year <- c(Laake.start.year, .start.year)

# Both datasets contain 2006. Take it out from the recent one.
bf <- vs <- day <- watch.prop <- all.n <- all.obs <- array(dim = c(max(dim(Laake.data$n)[1], dim(.data$n)[1]), 
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
    watch.prop[1:1:dim(.data$watch.prop)[1], 1, c] <- .data$watch.prop[, c2]
    watch.prop[1:1:dim(.data$watch.prop)[1], 2, c] <- .data$watch.prop[, c2]
    day[1:1:dim(.data$day)[1], 1, c] <- .data$day[, c2]
    day[1:1:dim(.data$day)[1], 2, c] <- .data$day[, c2]
    bf[1:1:dim(.data$bf)[1], 1, c] <- .data$bf[, c2]
    bf[1:1:dim(.data$bf)[1], 2, c] <- .data$bf[, c2]
    vs[1:1:dim(.data$vs)[1], 1, c] <- .data$vs[, c2]
    vs[1:1:dim(.data$vs)[1], 2, c] <- .data$vs[, c2]
    c <- c + 1
    c2 <- c2 + 1
  }
  
}

n.station <- c(Laake.data$n.station, .data$n.station[2:9])
n.year <- length(all.start.year)
n.obs <- length(unique(as.vector(all.obs))) - 1  # minus NA

tmp <- matrix(nrow = length(.data$periods)-1, ncol = 2)
tmp[,1] <- .data$periods[2:length(.data$periods)]
tmp[.data$n.station[2:9] == 2, 2] <- tmp[.data$n.station[2:9] == 2, 1]
periods <- rbind(Laake.data$periods, tmp)

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
                  day = day)

jags.params <- c("OBS.RF", "OBS.Switch",
                 "BF.Switch", "BF.Fixed",
                 "VS.Switch", "VS.Fixed",
                 "mean.prob", "mean.N", "Max",
                 "Corrected.Est", "Raw.Est", "N",
                 "K", "S1", "S2", "P",
                 "Max.alpha", "Max.beta",
                 "S1.alpha", "S2.alpha",
                 "S1.beta", "S2.beta",
                 "P.alpha", "P.beta",
                 "K.alpha", "K.beta",
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
 


# Look at Rhat statistics
max.Rhat <- lapply(jm.out$jm$Rhat, FUN = max, na.rm = T) %>%
  unlist()
max.Rhat.big <- max.Rhat[which(max.Rhat > 1.1)]


mcmc_dens(jm.out$jm$samples, c("Max.alpha", "Max.beta",
                               "S1.alpha", "S1.beta",
                               "S2.alpha", "S2.beta",
                               "P.alpha", "P.beta",
                               "K.alpha", "K.beta"))

mcmc_dens(jm.out$jm$samples, c("BF.Fixed", "VS.Fixed"))

par.idx <- c(1:nrow(jm.out$jm$mean$P))

mcmc_trace(jm.out$jm$samples, paste0("P[", par.idx, "]"))
mcmc_trace(jm.out$jm$samples, paste0("K[", par.idx, "]"))
mcmc_trace(jm.out$jm$samples, paste0("S1[", par.idx, "]"))
mcmc_trace(jm.out$jm$samples, paste0("S2[", par.idx, "]"))
mcmc_trace(jm.out$jm$samples, paste0("Max[", par.idx, "]"))


# Look at the abudance estimates:
Nhat. <- data.frame(median = jm.out$jm$q50$Corrected.Est,
                    LCL = jm.out$jm$q2.5$Corrected.Est,
                    UCL = jm.out$jm$q97.5$Corrected.Est,
                    mean = jm.out$jm$mean$Corrected.Est,
                    Year =  all.start.year + 1,
                    Season = paste0(all.start.year, "/", all.start.year+1))

N.hats <- data.frame(Season = rep(Nhat.$Season, each = nrow(jm.out$jm$mean$N)),
                     Day = rep(1:nrow(jm.out$jm$mean$N), times = length(Nhat.$Season)),
                     Mean = as.vector(jm.out$jm$mean$N),
                     LCL = as.vector(jm.out$jm$q2.5$N),
                     UCL = as.vector(jm.out$jm$q97.5$N))

ggplot(N.hats %>% group_by(Season)) + 
  geom_ribbon(aes(x = Day, ymin = LCL, ymax = UCL),
              fill = "blue", alpha = 0.5) +
  geom_path(aes(x = Day, y = Mean)) + 
  facet_wrap(~ Season)

# These are not the best estimates because they were not updated as more data
# were collected. I should use the output from the most recent WinBUGS run for 
# the last 9 years.
Reported.estimates <- read.csv(file = "Data/all_estimates_2023.csv") %>%
  select(-c(X, Year))

WinBugs.9yr.out <- readRDS(file = "RData/WinBUGS_9yr_v2_min85.rds")
Corrected.Est <- WinBugs.9yr.out$BUGS_out$sims.list$Corrected.Est
seasons <- c("2006/2007", "2007/2008", "2009/2010", "2010/2011", 
             "2014/2015", "2015/2016", "2019/2020", "2021/2022",
             "2022/2023")

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
                                              FUN = quantile, 0.975),
                                  Method = "Durban")

Reported.estimates %>%
  filter(Method == "Laake") %>%
  rbind(Durban.abundance.df) -> previous.estimates


Nhat. %>%
  left_join(previous.estimates, by = "Season") %>%
  mutate(delta_Nhat = mean - Nhat) -> all.estimates
# 
ggplot(all.estimates) +
  geom_point(aes(x = Season, y = mean ),
             color = "blue") +
  geom_errorbar(aes(x = Season, ymin = LCL.x, ymax = UCL.x),
                color = "blue") +
  geom_point(aes(x = Season, y = Nhat ),
             color = "green") +
  geom_errorbar(aes(x = Season, ymin = LCL.y, ymax = UCL.y),
                color = "green")
