# Data2WinBUGS_input.R
# Converts Granite Canyon count data to WinBUGS inputs. All raw data files 
# should be treated by Extract_Data_All_v2.Rmd. All output files should be in
# one directory, e.g., V2.1_Nov2024. 

rm(list=ls())
library(R2jags)
library(abind)
library(R2WinBUGS)
library(tidyverse)

# this file contains all necessary inputs for 2006 - 2019:
data.0 <- readRDS("RData/2006-2019_GC_Formatted_Data.RDS")

# data directory name
data.dir <- "RData/V2.1_Nov2024"

all.years <- c(2007, 2008, 2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024)

# output from Ver2.0 extraction - end year
years <- c("2010", "2011", "2015", "2016", "2020", "2022", "2023", "2024")

# can be 85, 30, 10, or anything that was assigned at the time of running
# Extract_Data_All_v2.Rmd
min.duration <- 10# 85 #30

# Change file names accordingly:
Nhats.filename <- paste0("Data/abundance_", years[length(years)], "_", 
                         min.duration, "min.csv")

out.v2 <- lapply(years, 
                 FUN = function(x) readRDS(paste0(data.dir, "/out_", x,
                                                  "_min", min.duration, "_Tomo_v2.rds")))

begin. <- lapply(out.v2, FUN = function(x) x$Final_Data$begin)
end. <- lapply(out.v2, FUN = function(x) x$Final_Data$end)

# Number of watch periods in each year's survey - before the 2019/2020 season
# plus the new ones
# This has been changed for 2024. I now have edited data for 2010 - 2024 seasons.
# So, I can just use the first two (2006/2007 and 2007/2008). The two numbers for 
# 2009/2010 and 2010/2011 don't match. In Durban's analysis, they were 164 and 178,
# respectively. For the new extraction, they are 180 and 145. 
#periods <-c(136, 135, 164, 178,
periods <-c(136, 135,
            lapply(begin., FUN = function(x) length(x)) %>% unlist)

# Shorten data to first x years only to replicate 
# the analysis in Durban et al 2016. 
# Or use it for other purposes.
x <- length(periods)

out.file.name <- paste0("RData/WinBUGS_", x, "yr_v2_min", min.duration, ".rds")
#out.file.name <- paste0("RData/WinBUGS_", x, "yr_v2_min", min.duration, "_8weeks_2023.rds")

Watch.Length. <- list()

for (k in 1:length(begin.)){
  Watch.Length.[[k]] <- end.[[k]] - begin.[[k]]
}

# I don't have edited data for 2006/2007 and 2007/2008. So, they need to be
# used as they were given in the WinBUGS input file from Durban. That's why
# I take the first two columns. It used to be four columns but I received edited
# data for 2010/2011 and 2011/2012. 
# 
# #whale counts
n <- labelled::remove_attributes(data.0$n, "dimnames")
n <- n[,,1:2]  # take the first two years (2007/2008 and 2008/2009)
n <- abind(n, array(0, replace(dim(n), 1, max(periods) - nrow(n))), along = 1)

# the u data is whether there were observers on watch. 
# 0 counts are often associated with years/shifts with 
# no second observer. So if u=0, it will fix observation probability at 0
# the second column for each year is for the second station - not the second
# observer.
u <- data.0$u[,,1:2]
u <- abind(u, 
           array(0, replace(dim(u), 1, max(periods) - nrow(u))), along = 1)

# #visibility
vs. <- lapply(out.v2, FUN = function(x) x$Final_Data$vs)

vs <- rbind(data.0$vs[,1:2],
                  array(NA, dim = c(max(periods) - nrow(data.0$vs), 2))) %>%
    labelled::remove_attributes("dimnames")

# Beaufort
bf. <- lapply(out.v2, FUN = function(x) x$Final_Data$bf)
bf <- rbind(data.0$bf[,1:2],
                  array(NA, dim = c(max(periods) - nrow(data.0$bf), 2)))%>%
    labelled::remove_attributes("dimnames")

# observers
# need to convert observer initials into numbers for all years
# This needs to be redone... 
obs.list <- read.csv(file = "Data/ObserverList2023.csv")

obs.2024 <- unique(out.v2[[length(years)]]$Complete_Data$obs)
new.obs <- obs.2024[!c(obs.2024 %in% obs.list$obs)]
obs.list <- rbind(obs.list,
                  data.frame(obs = new.obs,
                             ID = seq(max(obs.list$ID)+1, 
                                      max(obs.list$ID) + length(new.obs))))

obs <- data.0$obs[,,1:2]
obs <- abind(obs, 
           array(36, replace(dim(obs), 1, max(periods) - nrow(obs))), along = 1)

Watch.Length <- rbind(data.0$Watch.Length[,1:2],
                      array(NA, dim = c(max(periods) -
                                          nrow(data.0$Watch.Length), 2)))%>%
    labelled::remove_attributes("dimnames")

day <- rbind(data.0$day[,1:2],
                   array(NA, dim = c(max(periods) - nrow(data.0$day), 2)))%>%
    labelled::remove_attributes("dimnames")

k <- 1
for (k in 1:length(begin.)){
  # Final_Data is correct length and correct BF and VS
  n.stations <- length(unique(out.v2[[k]]$Final_Data$station))
  
  if (n.stations == 1){
    n.k <- cbind(c(out.v2[[k]]$Final_Data$n, 
                   rep(0, times = (dim(n)[1] - length(out.v2[[k]]$Final_Data$n)))), 
                 rep(0, times = dim(n)[1]))
  } else {
    n.1 <- out.v2[[k]]$Final_Data %>% filter(station == "P")# %>% select(n)
    n.2 <- out.v2[[k]]$Final_Data %>% filter(station == "S")# %>% select(n)
    n.k <- cbind(c())  
  }
  # n <- abind(n, 
  #            cbind(c(out.v2[[k]]$Final_Data$n, 
  #                    rep(0, times = (dim(n)[1] - length(out.v2[[k]]$Final_Data$n)))), 
  #                  rep(0, times = dim(n)[1]))) %>%
  #   labelled::remove_attributes("dimnames")
  
  n <- abind(n, n.k) %>%
    labelled::remove_attributes("dimnames")

  u <- abind(u, 
             cbind(c(rep(1, times = length(out.v2[[k]]$Final_Data$n)), 
                     rep(0, times = (dim(u)[1] - length(out.v2[[k]]$Final_Data$n)))), 
                   rep(0, times = dim(u)[1])))
  
  vs <- cbind(vs, c(out.v2[[k]]$Final_Data$vs, 
                    rep(NA, times = max(periods - length(out.v2[[k]]$Final_Data$vs)))))  %>%
    labelled::remove_attributes("dimnames")
  
  
  bf <- cbind(bf, c(out.v2[[k]]$Final_Data$bf,
                    rep(NA, times = max(periods - length(out.v2[[k]]$Final_Data$bf)))))  %>%
    labelled::remove_attributes("dimnames")
  
  obs.year <- data.frame(obs = out.v2[[k]]$Final_Data$obs) %>% 
    left_join(obs.list, by = "obs")
  
  obs <- abind(obs, 
               cbind(c(obs.year$ID, 
                       rep(36, times = (max(periods) - length(obs.year$ID)))), 
                     rep(36, times = max(periods))))

  Watch.Length <- cbind(Watch.Length,
                        c(Watch.Length.[[k]], 
                          rep(NA, 
                              times = max(periods) - length(Watch.Length.[[k]])))) %>%
    labelled::remove_attributes("dimnames")
  
  day <- cbind(day, 
               c(floor(begin.[[k]]), 
               rep(NA, times = max(periods) - length(begin.[[k]])))) %>%
    labelled::remove_attributes("dimnames")
  
}


#we're going to make N a partially observed data object with anchor points at day 1 and 90
# TE: I don't know how these numbers were created... they are generally 2x n (not all)
# N_inits <- as.matrix(read.table("Data/Initial Values/N_inits.txt",
#                                 header=T))

N_inits1 <- n[, 1,] * 2 + 2
N_inits2 <- n[, 2,] * 2 + 2 
            
N_inits <- N_inits1
N_inits[N_inits1 < N_inits2] <- N_inits2[N_inits1 < N_inits2]

N_inits <- rbind(N_inits,
                 matrix(data = NA, nrow = 2, ncol = length(periods)))

for (k in 1:length(periods)){
  N_inits[(periods[k]+1):nrow(N_inits), k] <- NA  
}

#The 'data' has to be the inverse of the inits, 
# with NAs for all of the estimated Ns, and 0s for the days 1 and 90
N <- matrix(NA, nrow=max(periods)+2, ncol=length(periods)) 

for(i in 1:length(periods)){
  N[(periods[i]+1):(periods[i]+2),i] <- 0 #True number of whales passing fixed at 0 for day 1 and 90
}

end <- Watch.Length + day

#t <- round((begin+end)/2)
t <- round((day+end)/2)

# #Add a couple of extra rows of NAs to the end of the day index reference to match up with the fixed 0s in N (above), assigning them to days 1 and 90
day <- rbind(as.matrix(day),
             matrix(NA, nrow=2, ncol=length(periods)))
#
#
for(i in 1:length(periods)){ #Set the anchor points: days 1 and 90
  day[(periods[i]+1):(periods[i]+2),i] <- c(1,90)
}

t <- rbind(as.matrix(t),
           matrix(NA,nrow=2,ncol=length(periods)))
for(i in 1:length(periods)){ #Set the anchor points: days 1 and 90
  t[(periods[i]+1):(periods[i]+2),i] <- c(1,90)
}

#Place 36s for 'no observer' for the two periods following the end of true watches (this is for the day 1 and day 90 zero-whale anchor points)
#this will force it to the mean observation probability with no observer effect

Watch.Length <- rbind(as.matrix(Watch.Length),
                      matrix(NA, nrow=2, ncol=length(periods)))

for(i in 1:length(periods)){
  Watch.Length[(periods[i]+1):(periods[i]+2),i] <- 1
}

BUGS.data <- list(n = n[1:max(periods[1:x]),,1:x],
                  n.com = n[1:max(periods[1:x]),,1:x],
                  n.sp = n[1:max(periods[1:x]),,1:x],
                  n.station = dim(n[1:max(periods[1:x]),,1:x])[2],
                  n.year = dim(n[1:max(periods[1:x]),,1:x])[3],
                  n.obs = max(obs[1:max(periods[1:x]),,1:x]),
                  periods = periods[1:x],
                  obs = obs[1:max(periods[1:x]),,1:x],
                  #Watch.Length = 0.0625,
                  u = u[1:max(periods[1:x]),,1:x],
                  vs = vs[1:max(periods[1:x]),1:x],
                  bf = bf[1:max(periods[1:x]),1:x],
                  #day=day,
                  day = t[1:(max(periods[1:x])+2),1:x],
                  N = N[,1:x],
                  N.com = N[,1:x],
                  N.sp = N[,1:x],
                  knot = c(-1.46,-1.26,-1.02,-0.78,
                         -0.58,-0.34,-0.10,0.10,
                         0.34,0.57,0.78,1.02,1.26,1.46),
                  n.knots=14,
                  #begin=begin,
                  #end=end,
                  Watch.Length=Watch.Length[1:(max(periods[1:x])+2), 1:x])

BUGS.inits <- function() list(mean.prob = 0.5,
                              BF.Fixed = 0,
                              VS.Fixed = 0,
                              mean.prob.sp = 0.5,
                              BF.Fixed.sp = 0,
                              VS.Fixed.sp = 0,
                              mean.prob.com = 0.5,
                              BF.Fixed.com = 0,
                              VS.Fixed.com = 0,
                              mean.beta = c(0,0,0), #mean.beta = c(5,0.14,-3.5),
                              beta.sigma = c(1,1,1),#beta.sigma = c(7,7,7),
                              BF.Switch = 1,
                              VS.Switch = 1,
                              OBS.Switch = 1,
                              sigma.Obs = 1,
                              BF.Switch.sp = 1,
                              VS.Switch.sp = 1,
                              OBS.Switch.sp = 1,
                              sigma.Obs.sp = 1,
                              BF.Switch.com = 1,
                              VS.Switch.com = 1,
                              OBS.Switch.com = 1,
                              sigma.Obs.com = 1,
                              N = N_inits,
                              N.com = N_inits,
                              N.sp = N_inits,
                              #z = matrix(1,nrow=90,ncol=6),
                              beta.sp = array(data=0, dim=c(2,x)),
                              sd.b.sp = rep(1, times = x), #c(1,1,1,1,1,1),
                              z = matrix(1, nrow=90, ncol= x))

