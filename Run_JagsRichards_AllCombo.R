# Run_JagsRichards_AllCombo.R
# 
# This runs all combinations of minimum watch duration and datasets and saves results.
# 

rm(list=ls())

source("Granite_Canyon_Counts_fcns.R")

data.dir <- "RData/V2.1_Nov2024"

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

MCMC.params <- list(n.samples = 250,
                    n.thin = 2,
                    n.burnin = 200,
                    n.chains = 5)

new.years <- c(2010, 2011, 2015, 2016, 
               2020, 2022, 2023, 2024)

WinBUGS.n.stations <- c(1, 1, 2, 2, rep(1, times = 6))



min.durs <- c(10, 30, 85)
vers <- c("v5", "v4", "v3", "v1")

k1 <- k2 <- 1
for (k1 in 3:length(min.durs)){
  WinBUGS.outfile <- paste0("RData/", list.files(path = "RData", 
                                pattern = paste0("WinBUGS_2007to2024_v2_min",
                                                 min.durs[k1], "_")))
  # if (min.durs[k1] < 85){
  #   WinBUGS.outfile <- "RData/WinBUGS_2007to2024_v2_min30_2024-11-23.rds"
  # } else {
  #   WinBUGS.outfile <- "RData/WinBUGS_2007to2024_v2_min85_2024-11-23.rds"
  # }
  
  for (k2 in 1:length(vers)){
    Jags_Richards_LaakeData_fcn(min.durs[k1], vers[k2], jags.params, MCMC.params)   
    
    Jags_Richards_AllData_fcn(min.durs[k1], vers[k2], 
                              new.years, WinBUGS.out.file = WinBUGS.outfile,
                              WinBUGS.n.stations, 
                              data.dir, jags.params, MCMC.params)
    
    Jags_Richards_NoLaakeData_fcn(min.durs[k1], vers[k2], 
                                  WinBUGS.outfile, years = new.years, 
                                  n.stations = WinBUGS.n.stations, 
                                  data.dir, jags.params, MCMC.params)
    
    NoBUGS_Richards_fcn(min.durs[k1], vers[k2], years = new.years, 
                        data.dir, jags.params, MCMC.params)
    
    Jags_Richards_Since2010_fcn(min.durs[k1], vers[k2], years = new.years, 
                                data.dir, jags.params, MCMC.params)
  }

}
