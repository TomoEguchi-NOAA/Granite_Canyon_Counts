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
max.day <- 90

# These are the ending year of each season - different from the standard, which
# uses the beginning year. So, for example, 2022 in the following vector indicates
# for the 2021/2022 season. These are the "new" data.
# min.dur is one of 10, 30, 85 or any other that was used in Extract_Data_All_v2.Rmd
years <- c(2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024, 2025)
data.dir = "RData/V2.1_Feb2025"
input.list <- AllData2JagsInput_NoBUGS(min.dur, years = years, data.dir, max.day = max.day)  

jags.data <- input.list$jags.data

obs.vec <- as.vector(jags.data$obs) 
data.frame(obs = obs.vec) %>%
  mutate(obs.f = as.factor(obs)) %>%
  group_by(obs.f) %>%
  summarize(n = n(),
            obs = first(obs)) -> obs.summary

obs.n.min <- 10
obs.too.few <- obs.summary %>% filter(n < obs.n.min)  

obs.to.keep <- obs.summary %>% filter(n >= obs.n.min) 
obs.to.keep$new.ID <- seq(1, dim(obs.to.keep)[1])
obs.others <- max(obs.to.keep$new.ID)

obs <- jags.data$obs
new.no.obs <- obs.others + 1
old.no.obs <- max(obs.to.keep$obs)
for (k in 1:nrow(obs.too.few)){
  obs[obs == obs.too.few$obs[k]] <- NA
}

for (k in 1:(nrow(obs.to.keep)-1)){
  obs[obs == obs.to.keep$obs[k]]  <- obs.to.keep$new.ID[k] 
}

obs[is.na(obs)] <- obs.others
obs[obs == old.no.obs] <- new.no.obs



n <- input.list$jags.data$n[2:(max(input.list$jags.data$periods[,1]) - 2), ,]
day <- input.list$jags.data$day[2:(max(input.list$jags.data$periods[,1]) - 2), ,]
day[day == max.day] <- NA
day.new <- array(dim = c(max(input.list$jags.data$periods[,1]), 2, dim(input.list$jags.data$n)[3]))
k1 <- k2 <- 1
for (k1 in 1:2){
  for (k2 in 1:dim(input.list$jags.data$n)[3]){
    day.1 <- na.omit(day[, k1, k2])
    if (length(day.1) > 1)
      day.new[1:(length(day.1)+2),k1,k2] <- c(day.1, 1, 90)    
  }
}

# jags data includes day 1 and max.day
periods <- input.list$jags.data$periods[,1] - 2



WinBUGS.input <- data2WinBUGS_input(data.dir = "RData/V2.1_Feb2025",
                                    years = c(2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024, 2025),
                                    min.dur = 60) 

# In June 2025, I changed the jags input to move day 1 and max.day to the beginning
# and the end of day and n arrays. So, that needs to be changed for WinBUGS input. 



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


periods <- colSums(!is.na(day.mat))

# Create the u array by replacing the non-observer ID (n.obs + 1) with 0s
# and observers with 1s
u <- input.list$jags.data$obs
u[u > max(input.list$jags.data$n.obs)] <- 0
u[u != 0] <- 1
u <- u[1:max(periods),,]

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

# need to add 1s where days are 1 and 90:
Watch.Length[day.mat == 1] <- 1
Watch.Length[day.mat == 90] <- 1


obs <- input.list$jags.data$obs
obs <- obs[1:max(periods),,]

BUGS.data <- list(n = n,
                  n.com = n,
                  n.sp = n,
                  n.station = max(input.list$jags.data$n.station),
                  n.year = input.list$jags.data$n.year,
                  n.obs = input.list$jags.data$n.obs + 1, # +1 for no-observer ID
                  periods = periods,
                  obs = obs,
                  u = u,
                  vs = input.list$jags.data$vs[,1,],
                  bf = input.list$jags.data$bf[,1,],
                  day = day.mat,
                  N = N,
                  N.com = N,
                  N.sp = N,
                  knot = c(-1.46,-1.26,-1.02,-0.78,
                           -0.58,-0.34,-0.10,0.10,
                           0.34,0.57,0.78,1.02,1.26,1.46),
                  n.knots = 14,
                  Watch.Length = Watch.Length) 

# Check data array dimensions
t.max <- BUGS.data$n.year
# periods.max <- BUGS.data$periods + 2
# dim.day <- dim(BUGS.data$day)
# dim.Watch.Length <- dim(BUGS.data$Watch.Length)
# dim.n <- dim(BUGS.data$n)
# dim.N <- dim(BUGS.data$N)
# for (k in 1:t.max) print(BUGS.data$Watch.Length[periods.max[k]-1, k])

n.obs <- BUGS.data$n.obs
dim.obs <- dim(BUGS.data$obs)
periods <- BUGS.data$periods
# for (k in 1:t.max) print(BUGS.data$obs[1:periods[k], 1, k])
# for (k in 1:t.max) print(BUGS.data$obs[1:periods[k], 2, k])

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

x <- BUGS.data$n.year
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
                              beta.sp = array(data=0, dim=c(2,x)),
                              sd.b.sp = rep(1, times = x), #c(1,1,1,1,1,1),
                              z = matrix(1, nrow=90, ncol= x))


MCMC.params <- list(n.iter = 200, #85000, #
                   n.thin = 2, #50, #80
                   n.burnin = 50, #50000, #
                   n.chains = 5)
# 
# This set runs quickly. Good to use this set to see if data and model are correct.
# MCMC.params <- list(n.iter = 200,
#                    n.thin = 2,
#                    n.burnin = 50,
#                    n.chains = 5)

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
  
  # node dimension does not match 2025-06-02 - solved Watch.Length had an extra dimension.
  # in the WinBUGS approach, primary and secondary stations are assumed to have the 
  # same watch length, even though they are not the same sometimes... 
  # 
  # array index is greater than array upper bound for OBS.RF.sp 2025-06-03 - solved
  # in the JAGS approach, I don't use "no-observer" ID. So, n.obs in the WinBUGS
  # approach needs n.obs to be 1 more than n.obs from the JAGS input
  # 
  # made use of undefined node Watch.Length 2025-06-03 - solved. Watch.Length
  # for Laake's data did not have 1 in day 1 and day 90. Once they were entered,
  # it runs fine. 
  Run_Time <- Sys.time() - Start_Time
  Ver2.results <-list(BUGS.input = input.list,
                      BUGS.data = BUGS.data,
                      BUGS.out = BUGS_out,
                      MCMCparams = MCMC.params,
                      Run_Time = Run_Time,
                      Run_Date = Run.date,
                      Sys.info = Sys.info())
  
  saveRDS(Ver2.results, file = out.file.name)
  
}


# Extract estimated counts
Daily.Est <- Ver2.results$BUGS.out$sims.list$Daily.Est
sp <- Ver2.results$BUGS.out$sims.list$sp
com <- Ver2.results$BUGS.out$sims.list$com
Corrected.Est <- Ver2.results$BUGS.out$sims.list$Corrected.Est

# Each one of them is (# samples) x (90 days) x (# years)
# To plot them using ggplot's facet, I need to convert
# these into 2D dataframes of statistics (upper and lower 
# CIs, median, etc.)
# Daily.Est.list <- sp.list <- com.list <- vector(mode = "list", 
#                                                 length = dim(Daily.Est)[3])
# 
# Daily.Est.UCIs <- Daily.Est.LCIs <- vector(mode = "list",
#                                            length = dim(Daily.Est)[3])

all.seasons <- c(input.list$jags.input.Laake$seasons,
                 input.list$jags.input.new$seasons)

stats.list <- vector(mode = "list",
                     length = dim(Daily.Est)[3])

for (k in 1:dim(Daily.Est)[3]){
  # Daily.Est.list[[k]] <- Daily.Est[,,k]
  # Daily.Est.UCIs[[k]] <- apply(Daily.Est[,,k],2,quantile,0.975)
  # Daily.Est.LCIs[[k]] <- apply(Daily.Est[,,k],2,quantile,0.275)
  # 
  # sp.list[[k]] <- sp[,,k]
  # com.list[[k]] <- com[,,k]
  
  stats.list[[k]] <- data.frame(Daily.Est.median = apply(Daily.Est[,,k], 2,
                                                         median),
                                Daily.Est.LCL = apply(Daily.Est[,,k], 2,
                                                      quantile,0.275),
                                Daily.Est.UCL = apply(Daily.Est[,,k], 2,
                                                      quantile,0.975),
                                sp.median = apply(exp(sp[,,k]), 2,
                                                  median),
                                sp.LCL = apply(exp(sp[,,k]), 2,
                                               quantile,0.025),
                                sp.UCL = apply(exp(sp[,,k]), 2,
                                               quantile,0.975),
                                com.median = apply(exp(com[,,k]), 2,
                                                   median),
                                com.LCL = apply(exp(com[,,k]), 2,
                                                quantile,0.025),
                                com.UCL = apply(exp(com[,,k]), 2,
                                                quantile,0.975),
                                #total.median = apply(exp(sp[,,k]), 1, sum),
                                days = 1:dim(Daily.Est)[2],
                                year = all.seasons[k])
}

all.stats <- do.call("rbind", stats.list) %>% group_by(year)

obs.day <- BUGS.data$day
n.obs.primary <- BUGS.data$n.com[,1,]
Nhats <- list()
k <- 1
for (k in 1:ncol(obs.day)){
  tmp.day <-  BUGS.data$day[,k] %>% 
    na.omit() 
  tmp.n <- n.obs.primary[1:(length(tmp.day)-2),k]
  tmp.watch <- BUGS.data$Watch.Length[1:(length(tmp.day)-2),k] %>% na.omit()
  Nhats[[k]] <- data.frame(day = tmp.day[tmp.day > 1 & tmp.day < 90],
                           n = tmp.n,
                           watch = tmp.watch) %>%
    group_by(day) %>%
    summarize(day = first(day),
              n = sum(n),
              watch = sum(watch),
              Nhat = n/watch,
              year = all.seasons[k])
  
}

Nhats.df <- do.call("rbind", Nhats) %>% group_by(year)

ggplot(data = all.stats) + 
  geom_line(aes(x = days, y = sp.median),
            color = "darkorange") + 
  geom_ribbon(aes(x = days, 
                  ymin = sp.LCL, 
                  ymax = sp.UCL),
              fill = "orange", 
              alpha = 0.5) +
  geom_line(aes(x = days, y = com.median),
            color = "darkred") +
  geom_ribbon(aes(x = days, 
                  ymin = com.LCL, 
                  ymax = com.UCL),
              fill = "red", 
              alpha = 0.5) +
  geom_point(aes(x = days, y = Daily.Est.median),
             shape = 16, fill = "blue", alpha = 0.3) +
  geom_point(data = Nhats.df,
             aes(x = day, y = Nhat),
             shape = 4, color = "darkgreen", alpha = 0.5) +
  facet_wrap(vars(year)) +
  xlab("Days since December 1") + 
  ylab("Whales per day")


#Com vs Sp
ggplot(data = all.stats) +
  geom_line(aes(x = days, y = sp.median)) +
  geom_ribbon(aes(x = days, 
                  ymin = sp.LCL, 
                  ymax = sp.UCL),
              fill = "orange", alpha = 0.5) + 
  facet_wrap(vars(year))+
  xlab("Days since December 1") + 
  ylab("Whales per day (spline)")


abundance.df <- data.frame(Season = all.seasons,
                           total.mean = apply(Corrected.Est,
                                              FUN = mean,
                                              MARGIN = 2),
                           total.CV = apply(Corrected.Est,
                                            FUN = function(x) 100*sqrt(var(x))/mean(x),
                                            MARGIN = 2),
                           total.median = apply(Corrected.Est, 
                                                FUN = median, 
                                                MARGIN = 2),
                           total.LCL = apply(Corrected.Est, 
                                             MARGIN = 2, 
                                             FUN = quantile, 0.025),
                           total.UCL = apply(Corrected.Est, 
                                             MARGIN = 2, 
                                             FUN = quantile, 0.975))


# write.csv(abundance.df,
#           file = paste0("Data/WinBUGS_abundance_", min(WinBUGS.input$all.years), "to",
#                         max(WinBUGS.input$all.years), "_v2_min", WinBUGS.input$min.dur,
#                         "_", MCMCparams$n.iter, "_", Run_Date, ".csv"))

ggplot(data = abundance.df) + 
  geom_point(aes(x = Season, y = total.median)) + 
  geom_errorbar(aes(x = Season, ymin = total.LCL, ymax = total.UCL))

