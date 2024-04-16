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


out.file.name <- "RData/JAGS_Spline_results_All_Data.rds"

jm.out <- readRDS(out.file.name)

jm.out$jags.data$day.1 %>% 
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = "name") %>% 
  mutaet(day = value) -> tmp.1 

tmp.1$Season <-  lapply(strsplit(tmp.1$name, "V"), 
                        FUN = function(x) x[2]) %>% unlist()

tmp.1 %>% select(Season, day) %>%
  arrange(Season, day) -> sampled.day.1 
