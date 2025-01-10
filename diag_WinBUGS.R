# diag_WinBUGS.R
# Runs diagnostics of WinBUGS runs using the coda package


rm(list=ls())
# Run the Run_ script to get all input information
source("Run_JagsRichards_AllCombo.R")

library(loo)
library(tidyverse)
library(R2WinBUGS)
library(coda)

# Bring in WinBUGS results:
k1 <- 1
WinBUGS.out.file <- paste0("RData/", list.files(path = "RData", 
                                                pattern = paste0("WinBUGS_2007to2024_v2_min",
                                                                 min.durs[k1], "_85000_")))

WinBUGS.out <- readRDS(WinBUGS.out.file)
max.Rhat <- max(WinBUGS.out$BUGS.out$summary[,"Rhat"])

as.data.frame(WinBUGS.out$BUGS.out$summary) %>% 
  rownames_to_column(var = "parameter") %>%
  filter(Rhat > 1.2) -> large.Rhat
