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

out.file.name <- paste0("RData/JAGS_Richards_pois_bino_v3_AllYears_",
                        Sys.Date(), ".rds")

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

# Look at Rhat statistics
max.Rhat <- lapply(jm.out$jm$Rhat, FUN = max, na.rm = T) %>%
  unlist()
max.Rhat.big <- max.Rhat[which(max.Rhat > 1.1)]



Nhat. <- data.frame(median = jm.out$jm$q50$Corrected.Est,
                      LCL = jm.out$jm$q2.5$Corrected.Est,
                      UCL = jm.out$jm$q97.5$Corrected.Est,
                      mean = jm.out$jm$mean$Corrected.Est,
                      Season =  all.start.year + 1)

# Nhat.df %>% 
#   left_join(Laake.estimates, by = "Season") %>%
#   mutate(delta_Nhat = median - Nhat) -> all.estimates
# 
# ggplot(all.estimates) +
#   geom_point(aes(x = Season, y = median ),
#              color = "blue") +
#   geom_errorbar(aes(x = Season, ymin = LCL.x, ymax = UCL.x),
#                 color = "blue") +
#   geom_point(aes(x = Season, y = Nhat ),
#              color = "green") +
#   geom_errorbar(aes(x = Season, ymin = LCL.y, ymax = UCL.y),
#                 color = "green")

# mcmc_dens(jm.out.Laake$jm$samples, c("Max.alpha", "Max.beta", 
#                                      "S1.alpha", "S1.beta", 
#                                      "S2.alpha", "S2.beta",
#                                      "P.alpha", "P.beta"))

P1 <- c("P[1]", "P[2]", "P[3]", "P[4]", "P[5]", "P[6]", "P[7]", "P[8]", "P[9]")
P2 <- c("P[10]", "P[11]", "P[12]", "P[13]", "P[14]", "P[15]", "P[16]", "P[17]", "P[18]")
P3 <- c("P[19]", "P[20]", "P[21]", "P[22]", "P[23]")

K3 <- c("K[19]", "K[20]", "K[21]", "K[22]", "K[23]")
K2 <- c("K[10]", "K[11]", "K[12]", "K[13]", "K[14]", "K[15]", "K[16]", "K[17]", "K[18]")
K1 <- c("K[1]", "K[2]", "K[3]", "K[4]", "K[5]", "K[6]", "K[7]", "K[8]", "K[9]")

S1.1 <- c("S1[1]", "S1[2]", "S1[3]", "S1[4]", "S1[5]", "S1[6]", "S1[7]", "S1[8]", "S1[9]")
S2.1 <- c("S2[1]", "S2[2]", "S2[3]", "S2[4]", "S2[5]", "S2[6]", "S2[7]", "S2[8]", "S2[9]")

S1.2 <- c("S1[10]", "S1[11]", "S1[12]", "S1[13]", "S1[14]", "S1[15]", "S1[16]", "S1[17]", "S1[18]")
S2.2 <- c("S2[10]", "S2[11]", "S2[12]", "S2[13]", "S2[14]", "S2[15]", "S2[16]", "S2[17]", "S2[18]")

S1.3 <- c("S1[19]", "S1[20]", "S1[21]", "S1[22]", "S1[23]")
S2.3 <- c("S2[19]", "S2[20]", "S2[21]", "S2[22]", "S2[23]")

