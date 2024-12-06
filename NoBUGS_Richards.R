# NoBUGS_Richards.R
# 
# Runs Richards function analysis without using WinBUGS input.

library(tidyverse)
library(abind)
source("Granite_Canyon_Counts_fcns.R")

min.dur <- 30
jags.input.Laake <- LaakeData2JagsInput(min.dur)

years <- c(2010, 2011, 2015, 2016, 
           2020, 2022, 2023, 2024)

jags.input.new <- data2Jags_input_NoBUGS(min.dur = min.dur, 
                                         years = years,
                                         n.stations = c(2, 2, rep(1, times = 6)), 
                                         data.dir = "RData/V2.1_Nov2024",
                                         run.date = Sys.Date())
Need to figure out how to merge observers between the two datasets. 2024-12-05
