#Data2LaakeAnalysis
# Creates input data objects for analayis in ERAnalysis pacakge. 
# This script is for data from the 2006/2007 season. Data extraction is done 
# using Extract_Data_All_v2.Rmd. This script uses the output of the markdown,
# which is an rds file for each year. 


rm(list = ls())
library(ERAnalysis)
library(tidyverse)

# I don't have raw data for 2007/2008, 2009/2010, 2010/2011. 
begin.years <- c(2014, 2015, 2019, 2021, 2022)

in.dir <- "RData/V2.1_Sep2023"

files <- list.files(in.dir, pattern = "_min85_Tomo_v2.rds")

data.list <- list()
k <- 1
for (k in 1:length(files)){
  file.parts <- strsplit(files[k], split = "_") %>% unlist()
  year <- as.numeric(file.parts[2])
  tmp <- readRDS(paste0(in.dir, "/", files[k])) 
  tmp$Final_Dat %>%
    transmute(Start.year = year - 1,
              begin = begin,
              end = end,
              nwhales = n,
              effort = dur,
              vis = vs,
              beaufort = bf,
              Observer = obs,
              watch = i,
              Date = BeginDay + as.Date(paste0(Start.year, "-11-30"))) -> data.list[[k]] 
}

all.data <- do.call(rbind, data.list)

