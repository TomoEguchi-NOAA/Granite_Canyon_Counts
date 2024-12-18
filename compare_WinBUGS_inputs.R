# compare WinBUGS inputs
# 
# Abundance estimates from WinBUGS in Feb 2024 are quite different from those in 
# December 2024. The differences must come from the input data because I updated
# functions to extract raw data and create WinBUGS inputs. In this script, I try
# to compare the input data for the two periods and see what changed.

rm(list = ls())
library(tidyverse)
library(ggplot2)

source("Granite_Canyon_Counts_fcns.R")

# Results from Feb 2024 data with minimum of 85 minutes observation periods
BUGS.out.min85.Feb2024 <- readRDS("RData/WinBUGS_10yr_v2_min85.rds")
BUGS.data.Feb2024 <- BUGS.out.min85.Feb2024$BUGS.data 
n.min85.Feb2024 <- BUGS.data.Feb2024$n
periods.Feb2024 <- BUGS.data.Feb2024$periods


# Results from Dec 2024 data with minimum of 85 minutes observation periods
# No data were included for the secondary counts
BUGS.out.min85.Dec2024 <- readRDS("RData/WinBUGS_2007to2024_v2_min85_100000_2024-11-23.rds")
BUGS.data.Dec2024 <- BUGS.out.min85.Dec2024$BUGS.input$data
n.Dec2024 <- BUGS.data.Dec2024$n
periods.Dec2024 <- BUGS.data.Dec2024$periods
years.Dec2024 <- BUGS.out.min85.Dec2024$BUGS.input$all.years

# the number of periods for 2010 and 2011 are a lot larger for the Dec 2024 dataset
# than for the Feb 2024 dataset - that may be the problem.

annual.n.min85.Feb2024 <- colSums(n.min85.Feb2024, na.rm = T)
annual.n.min85.Dec2024 <- colSums(n.Dec2024, na.rm = T)

# I see a problem. For Dec 2024 dataset, no observations were recorded for the
# secondary stations for 2010 and 2011. Also, there were more observations recorded
# for those two years

# fixed the function. Compare again:
input.data.min85 <- data2WinBUGS_input(data.dir = "RData/V2.1_Nov2024",
                                       years = c(2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024),
                                       min.dur = 85)

# This looks better. Total counts for the secondary stations were less than those
# in the Feb dataset because I removed some secondary observations that did not 
# have matching primary counts. This was because the WinBUGS code does not allow
# counts only from the secondary station. The watch length is the same between
# primary and secondary for the j-th period in t-th year.

extract.n <- function(BUGS.data, all.years){
  n.P <- BUGS.data$n[,1,] %>% as.data.frame()
  colnames(n.P) <- all.years
  
  n.P %>% 
    pivot_longer(everything(), 
                 names_to = "year", values_to = "n") %>%
    arrange(year) %>%
    mutate(day = as.vector(BUGS.data$day[1:nrow(n.P),]),
           station = "P") %>%
    na.omit() -> n.P.long
  
  n.S <- BUGS.data$n[,2,]  %>% as.data.frame()
  colnames(n.S) <- all.years
  
  n.S %>% 
    pivot_longer(cols = c("2010", "2011"), 
                 names_to = "year", values_to = "n") %>%
    select(year, n) %>%
    arrange(year) %>%
    mutate(day = as.vector(BUGS.data$day[1:nrow(n.S),3:4]),
           station = "S") %>%
    na.omit() -> n.S.long
  
  n.long <- rbind(n.P.long, n.S.long)
  
  p.daily <- ggplot(n.long) +
    geom_point(aes(x = day, y = n, color = station)) +
    facet_wrap(~year)  
  
  return(list(data = n.long,
              plot = p.daily))
}

Feb.2024 <- extract.n(BUGS.data.Feb2024, all.years = years.Dec2024)
Dec.2024 <- extract.n(BUGS.data.Dec2024, all.years = years.Dec2024)
new.data <- extract.n(input.data.min85$data, all.years = years.Dec2024)
