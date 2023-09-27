#SplineJags_2016to2023
# Analyzed recent data (from 2016) with spline only model

rm(list = ls())
library(tidyverse)
library(jagsUI)
library(ggplot2)

seasons <- c("2006/2007", "2007/2008", "2009/2010", "2010/2011", 
             "2014/2015", "2015/2016", "2019/2020", "2021/2022",
             "2022/2023")

# Bring in the output from a successful run of WinBUGS (a rare event...)
# BUGS input and output
x <- 9
min.duration <- 85
out.file.name <- paste0("RData/WinBUGS_", x, "yr_v2_min", min.duration, ".rds")
BUGS.9yr.results <- readRDS(out.file.name)
BUGS.data <- BUGS.9yr.results$BUGS.data

n.sp <- BUGS.data$n.sp
for (k in 1:dim(BUGS.data$n.sp)[3]){
  if (BUGS.data$periods[k] < nrow(BUGS.data$n.sp))
    n.sp[(BUGS.data$periods[k]+1):nrow(BUGS.data$n.sp), , k] <- NA
}

# Create jags input
jags.data <- list(n.sp = abind::abind(n.sp, 
                                      array(NA, dim = c(2, 2, dim(n.sp)[3])), 
                                      along = 1), 
                  N.sp = BUGS.data$N.sp,
                  n.station = c(1,1,2,2,1,1,1,1,1),
                  n.year = BUGS.data$n.year,
                  n.obs = BUGS.data$n.obs,
                  periods = BUGS.data$periods,
                  obs = BUGS.data$obs,
                  vs = BUGS.data$vs,
                  bf = BUGS.data$bf,
                  Watch.Length = BUGS.data$Watch.Length,
                  day = BUGS.data$day,
                  #u = u,
                  n.days = rep(90, times = dim(n.sp)[3]),
                  knot = c(-1.46,-1.26,-1.02,-0.78,
                           -0.58,-0.34,-0.10,0.10,
                           0.34,0.57,0.78,1.02,1.26,1.46),
                  n.knots = 14)

u <- array(data = 0, dim = dim(n.sp))

for (y in 1:BUGS.data$n.year){
  for (s in 1:jags.data$n.station[y]){
    for (p in 1:BUGS.data$periods[y]){
      
      u[p,s,y] <- 1
    }
  }
}

jags.data$u <- u

jags.params <- c("lambda.sp",
                 "OBS.RF.sp",
                 "OBS.Switch.sp",
                 "BF.Switch.sp",
                 "BF.Fixed.sp",
                 "VS.Switch.sp",
                 "VS.Fixed.sp",
                 "mean.prob.sp",
                 "BF.Fixed.sp",
                 "VS.Fixed.sp",
                 "Corrected.Est",
                 "Raw.Est",
                 "z",
                 "sp",
                 "Daily.Est","mean.beta",
                 "beta.sigma","beta","beta.sp","b.sp","sd.b.sp")

MCMC.params <- list(n.samples = 250000,
                    n.thin = 100,
                    n.burnin = 200000,
                    n.chains = 5)

# Needs initial values for N to avoid errors. N has to be significantly larger
# than observed n.
jags.inits <- function() list(N.sp = (jags.data$n.sp[,1,] * 2) + 10)

out.file.name <- "RData/JAGS_Spline_results_2006_2023.rds"
jags.model <- paste0("models/Model_Nmix_Spline.txt")


if (!file.exists(out.file.name)){
  Start_Time<-Sys.time()
  
  jm <- jagsUI::jags(jags.data,
                     inits = jags.inits,
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
                 jags.params = jags.params,
                 jags.model = jags.model,
                 MCMC.params = MCMC.params,
                 Run_Time = Run_Time,
                 System = Sys.getenv())
  
  saveRDS(jm.out,
          file = out.file.name)
  
} else {
  jm.out <- readRDS(out.file.name)
}

Nhat.df <- data.frame(median = jm.out$jm$q50$Corrected.Est,
                      LCL = jm.out$jm$q2.5$Corrected.Est,
                      UCL = jm.out$jm$q97.5$Corrected.Est,
                      mean = jm.out$jm$mean$Corrected.Est,
                      saesons = seasons)

ggplot(Nhat.df) +
  geom_point(aes(x = seasons, y = median )) +
  geom_errorbar(aes(x = seasons, ymin = LCL, ymax = UCL))
