# Run_JagsRichards_AllCombo.R
# 
# This runs all combinations of minimum watch duration and datasets and saves results.
# 
# To see analysis results, use Results_JagsRichards_AllCombos.R

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

# MCMC.params <- list(n.samples = 250,
#                     n.thin = 2,
#                     n.burnin = 200,
#                     n.chains = 5)

new.years <- c(2010, 2011, 2015, 2016, 
               2020, 2022, 2023, 2024)

WinBUGS.n.stations <- c(1, 1, 2, 2, rep(1, times = 6))

min.durs <- c(10, 30, 85)
vers <- c("v5", "v4", "v3", "v1")

WinBUGS.outfile <- out.file.name <- list()
c1 <- c2 <- k1 <- k2 <- 1
for (k1 in 1:length(min.durs)){
  WinBUGS.outfile[[c1]] <- paste0("RData/", list.files(path = "RData", 
                                                       pattern = paste0("WinBUGS_2007to2024_v2_min",
                                                                        min.durs[k1], "_85000_")))
  c1 <- c1 + 1
 
  min.dur <- min.durs[k1]
  
  for (k2 in 1:length(vers)){
    model.name <- paste0("Richards_pois_bino_", vers[k2])

    # Runs just Laake data with Richards function
    out.file.name[[c2]] <- list.files(path = "RData/", 
                                      pattern = paste0("JAGS_", model.name,"_min", min.dur,
                                                       "_LaakeData_"))
    
    if (length(out.file.name[[c2]]) == 0){
      Jags_Richards_LaakeData_fcn(min.durs[k1], vers[k2], jags.params, MCMC.params)   
    }
    
    c2 <- c2 + 1
    # Runs all data including Laake's and new data
    out.file.name[[c2]] <- list.files(path = "RData/",
                                      pattern = paste0("JAGS_", model.name,"_min", min.dur,
                                                       "_AllYears_"))
    if(length(out.file.name[[c2]]) == 0){
      Jags_Richards_AllData_fcn(min.durs[k1], vers[k2], 
                                new.years, 
                                WinBUGS.out.file = WinBUGS.outfile,
                                WinBUGS.n.stations, 
                                data.dir, jags.params, MCMC.params)
      
    }
    
    c2 <- c2 + 1
    # Runs just WinBUGS and new data without Laake's data
    out.file.name[[c2]] <- list.files(path = "RData/", 
                                      pattern = paste0("JAGS_", model.name,
                                                       "_min", min.dur,
                                                       "_Since2006_"))
    if (length(out.file.name[[c2]]) == 0){
      Jags_Richards_NoLaakeData_fcn(min.durs[k1], vers[k2], years = new.years,
                                    data.dir, jags.params, MCMC.params)
      
    }
    
    c2 <- c2 + 1
    # Runs over all data but not using WinBUGS input - all newly extracted
    out.file.name[[c2]] <- list.files(path = "RData/", 
                                      pattern = paste0("JAGS_", model.name,"_min", min.dur,
                                                       "_NoBUGS_"))
    if (length(out.file.name[[c2]]) == 0){
      NoBUGS_Richards_fcn(min.durs[k1], vers[k2], years = new.years, 
                          data.dir, jags.params, MCMC.params)}
    
    c2 <- c2 + 1
    # Runs data since 2010 without using WinBUGS input
    out.file.name[[c2]] <- list.files(path = "RData/",
                                paste0("JAGS_", model.name,"_min", min.dur,
                                       "_Since2010_NoBUGS_"))
    if(length(out.file.name[[c2]]) == 0){
      Jags_Richards_Since2010_fcn(min.durs[k1], vers[k2], years = new.years, 
                                  data.dir, jags.params, MCMC.params)}
  }
  c2 <- c2 + 1
}
