#CI_decline
# Computes 95% CI for the decline. Requested for Scordino's paper
# 2023-09-22

rm(list=ls())
library(tidyverse)

# Bring in the most recent analysis output from WinBUGS
# 8yr_v1: the median value for the 2021/2022 (16,650) was used in the report
#BUGS.out <- readRDS("RData/WinBUGS_8yr_v2.rds")  # v2 extraction using mine
# v1 extraction using Stewart's
BUGS.out <- readRDS("RData/WinBUGS_8yr_v1.rds")   

Corrected.Est.samples <- BUGS.out$BUGS_out$sims.list$Corrected.Est
median.Nhat <- apply(Corrected.Est.samples, MARGIN = 2, FUN = median)
CI.Nhat <- apply(Corrected.Est.samples, MARGIN = 2, FUN = quantile, c(0.025, 0.5, 0.975))

decline.rate <- (Corrected.Est.samples[,8] - Corrected.Est.samples[,6])/Corrected.Est.samples[,6]
quantile(decline.rate, c(0.025, 0.5, 0.975))

mean(decline.rate)
     