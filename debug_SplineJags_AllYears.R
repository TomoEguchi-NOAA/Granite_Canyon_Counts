# debug_SplineJags_AllYears.R
# This script is used to debug problems with SplineJags_AllYears.R, where the 
# estimated abundance seems to be a lot less than Laake's estimates when there
# were two observation stations.
# 

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(readr)
library(jagsUI)
library(lubridate)
library(loo)
library(bayesplot)
library(ERAnalysis)

# Bring in the output file from SplineJags_AllYears.R
out.file.name <- "RData/JAGS_Spline_results_All_Data.rds"

jm.out <- readRDS(out.file.name)

# Find sampling days for each year
find.sampled <- function(x, varname){
  # Renaming the variable name with the content of a variable:
  # https://forum.posit.co/t/have-rename-recognize-a-variable-value-as-a-column-name/169548
  
  x %>%
    as.data.frame() %>%
    pivot_longer(cols = everything(), names_to = "name") %>% 
    mutate(variable = value) %>%
    rename("{varname}" := variable) -> tmp.1 
  
  tmp.1$Season <-  lapply(strsplit(tmp.1$name, "V"), 
                          FUN = function(x) x[2]) %>% 
    unlist() %>% 
    as.numeric()
  
  tmp.1 %>% 
    select(-c(name, value)) %>%
    arrange(Season) -> sampled.day 

  return(sampled.day)  
}

# Pull out data and change the format to "long"
sampled.day.1 <- find.sampled(jm.out$jags.data$day.1, varname = "day")
sampled.day.2 <- find.sampled(jm.out$jags.data$day.2, varname = "day")

sampled.n.1 <- find.sampled(jm.out$jags.data$n.1, varname = "n")  
sampled.n.2 <- find.sampled(jm.out$jags.data$n.2, varname = "n")  

# the number of whales on day 1 and day 90 are NAs rather than zeros. 
# Could this be a problem? Turns out this was not the problem 2024-04-22

sampled.bf.1 <- find.sampled(jm.out$jags.data$bf.1, varname = "bf")  
sampled.bf.2 <- find.sampled(jm.out$jags.data$bf.2, varname = "bf")  

sampled.vs.1 <- find.sampled(jm.out$jags.data$vs.1, varname = "vs")  
sampled.vs.2 <- find.sampled(jm.out$jags.data$vs.2, varname = "vs")  

sampled.obs.1 <- find.sampled(jm.out$jags.data$obs.1, varname = "obs")  
sampled.obs.2 <- find.sampled(jm.out$jags.data$obs.2, varname = "obs")  

sampled.watch.1 <- find.sampled(jm.out$jags.data$Watch.Length.1, varname = "watch")  
sampled.watch.2 <- find.sampled(jm.out$jags.data$Watch.Length.2, varname = "watch")  

# Problem seems to be having wrong 1s in the u.2 matrix
sampled.u.1 <- find.sampled(jm.out$jags.data$u.1, varname = "u")  
sampled.u.2 <- find.sampled(jm.out$jags.data$u.2, varname = "u")  
