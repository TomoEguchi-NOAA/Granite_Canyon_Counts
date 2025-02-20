#create_ObserverIDs.R
#
# Creates observer IDs from Laake data to most recent year.
# 

source("Granite_Canyon_Counts_fcns.R")
library(ggplot2)
library(tidyverse)

# Obserer list from Laake's data
data(Observer)

# Bring in WinBUGS data from Durban
# this file contains all necessary inputs for 2006 - 2019 from Josh Stewart
# They coorespond to 2007, 2008, 2010, 2011, 2015, 2016, and 2020
data.0 <- readRDS("RData/2006-2019_GC_Formatted_Data.RDS")
obs.v0 <- data.0$obs[,1,]

# There are 35 unique observers for 2007 and 2008 data. 36 was used for no observer
# place holder. 
obs.v0.ID_2007to2008 <- obs.v0[,1:2]

obs.v0_2010to2020 <- cbind(obs.v0[,3:7], seq(1:nrow(obs.v0)))

# Bring in the data from edited files
years <- c(2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024, 2025)
all.years <- c(2007, 2008, years)

data.dir <- "RData/V2.1_Feb2025"
min.dur <- 85
out.v2 <- lapply(years, 
                 FUN = function(x) readRDS(paste0(data.dir, "/out_", x,
                                                  "_min", min.dur, 
                                                  "_Tomo_v2.rds")))

obs.v2.list <- lapply(out.v2, 
                 FUN = function(x) x$Final_Data$obs)

obs.v2_2010to2025 <- array(dim = c(lapply(obs.v2.list, 
                                          FUN = "length") %>% 
                                     unlist() %>% max(), 
                                   length(obs.v2.list)))

for (k in 1:length(obs.v2.list)){
  obs.v2_2010to2025[1:length(obs.v2.list[[k]]), k] <- obs.v2.list[[k]]
}

# Create a look up table for the V2, then compare their IDs to the v0 IDs 
obs.v2.ID.table <- data.frame(initial = str_sort(unique(as.vector(obs.v2_2010to2025)))) %>%
  rowid_to_column(var = "ID.v2")

# make "No obs" to be ID 36.
obs.v2.ID.table[obs.v2.ID.table$ID.v2 == 36, "ID.v2"] <- nrow(obs.v2.ID.table)

obs.v2.ID.table[is.na(obs.v2.ID.table$initial),"ID.v2"] <- 36
obs.v2.ID.table[is.na(obs.v2.ID.table$initial),"initial"] <- "No obs"

# Use this look up table to create the observer tables with corresponding IDs
obs.v2.ID <- apply(obs.v2_2010to2025, MARGIN = 2,
                   FUN = function(x, table){
                     y <- vector(mode = "numeric", length = length(x))
                     for (k in 1:length(table$initial)){
                       y[x == table$initial[k]] <- table$ID.v2[k]
                     }
                     y[y == 0] <- NA
                     return(y)
                   }, 
                   obs.v2.ID.table) 

# Add shift ID
obs.v2.ID_2010to2020 <- cbind(obs.v2.ID[,1:5], seq(1, nrow(obs.v2.ID)))
obs.v2.ID_2010to2020.long <- data.frame(obs.ID = as.vector(obs.v2.ID_2010to2020[,1:5]),
                                        year = rep(c(2010, 2011, 2015, 2016, 2020), 
                                                   each = nrow(obs.v2.ID_2010to2020)),
                                        shift = rep(obs.v2.ID_2010to2020[,6], 5),
                                        v = "v2") %>% na.omit() 

initial.vec <- vector(mode = "character", length = nrow(obs.v2.ID_2010to2020.long))
for (k in 1:nrow(obs.v2.ID.table)){
  initial.vec[obs.v2.ID_2010to2020.long$obs.ID == obs.v2.ID.table$ID.v2[k]] <- obs.v2.ID.table$initial[k]
}

obs.v2.ID_2010to2020.long$initial <- initial.vec

# Compare these new IDs to IDs from v.0
obs.v0.ID_2010to2020.long <- data.frame(obs.ID = as.vector(obs.v0_2010to2020[,1:5]),
                                        year = rep(c(2010, 2011, 2015, 2016, 2020), 
                                                   each = nrow(obs.v0_2010to2020)),
                                        shift = rep(obs.v0_2010to2020[,6], 5),
                                        v = "v0") %>% na.omit()
obs.v0.ID_2010to2020.long$initial <- NA

obs.ID_2010to2020 <- rbind(obs.v0.ID_2010to2020.long, 
                           obs.v2.ID_2010to2020.long)

ggplot(obs.ID_2010to2020) +
  geom_point(aes(x = shift, y = obs.ID, color = v)) +
  facet_wrap(~year)
