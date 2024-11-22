# plot_Jags_results.R
# Creates plots from results of Jags analysis
# It needs to run Jags_Richards_AllData.R or Jags_Richards_NoLaakeData.R.
# 
# Make sure to edit these scripts to match what is needed
# 

rm(list=ls())

.data <- "all" # or "no Laake"

if (.data == "all"){
  source("Jags_Richards_AllData.R")
} else {
  source("Jags_Richards_NoLaakeData.R")
}

# make all start years in numeric
start.years <- as.numeric(jags.input$start.years)

# Create a dataframe with all years, including unsampled years.
all.years <- data.frame(start.year = seq(min(start.years), 
                                   max(start.years))) %>%
  mutate(Season = paste0(start.year, "/", start.year + 1))

# Look at the annual abundance estimates:
Nhat. <- data.frame(Season = paste0(start.years, "/", 
                                    start.years+1),
                    Nhat = jm.out$jm$q50$Corrected.Est,
                    LCL = jm.out$jm$q2.5$Corrected.Est,
                    UCL = jm.out$jm$q97.5$Corrected.Est) %>%
  right_join(all.years, by = "Season") %>%
  arrange(start.year) %>%
  mutate(Method = "Richards")

# This is for daily estimates
N.hats.day <- data.frame(Season = rep(paste0(start.years, "/", start.years+1), 
                                      each = nrow(jm.out$jm$mean$N)), #rep(Nhat.$Season, each = nrow(jm.out$jm$mean$N)),
                         Day = rep(1:nrow(jm.out$jm$mean$N), times = length(start.years)),
                         Mean = as.vector(jm.out$jm$mean$N),
                         LCL = as.vector(jm.out$jm$q2.5$N),
                         UCL = as.vector(jm.out$jm$q97.5$N)) 

# Daily estimates plots
p.daily.Richards <- ggplot(N.hats.day %>% group_by(Season)) + 
  geom_ribbon(aes(x = Day, ymin = LCL, ymax = UCL),
              fill = "blue", alpha = 0.5) +
  geom_path(aes(x = Day, y = Mean)) + 
  facet_wrap(~ Season)


# These are not the best estimates because they were not updated as more data
# were collected. I should use the output from the most recent WinBUGS run for 
# the last x years.
Reported.estimates <- read.csv(file = "Data/all_estimates_2024.csv") %>%
  transmute(Season = Season,
            Nhat = Nhat,
            LCL = LCL,
            UCL = UCL,
            Method = paste0(Method, "-Reported")) %>%
  right_join(all.years, by = "Season") %>%
  arrange(start.year) %>%
  relocate(Method, .after = start.year)

WinBugs.out <- readRDS(file = "RData/WinBUGS_10yr_v2_min85.rds")

Corrected.Est <- WinBugs.out$BUGS_out$sims.list$Corrected.Est
WinBUGS.seasons <- c("2006/2007", "2007/2008", "2009/2010", "2010/2011", 
             "2014/2015", "2015/2016", "2019/2020", "2021/2022",
             "2022/2023", "2023/2024")

WinBUGS.abundance.df <- data.frame(Season = WinBUGS.seasons,
                                  Nhat = apply(Corrected.Est,
                                               FUN = mean,
                                               MARGIN = 2),
                                  # CV = apply(Corrected.Est,
                                  #            FUN = function(x) 100*sqrt(var(x))/mean(x),
                                  #            MARGIN = 2),
                                  # median = apply(Corrected.Est, 
                                  #                FUN = median, 
                                  #                MARGIN = 2),
                                  LCL = apply(Corrected.Est, 
                                              MARGIN = 2, 
                                              FUN = quantile, 0.025),
                                  UCL = apply(Corrected.Est, 
                                              MARGIN = 2, 
                                              FUN = quantile, 0.975)) %>%
  right_join(all.years, by = "Season") %>%
  arrange(year) %>%
  mutate(Method = "Durban")

# Create a dataframe for daily estimates:
daily.estim <- WinBugs.out$BUGS_out$sims.list$Daily.Est

# get stats:
mean.mat <- LCL.mat <- UCL.mat <- matrix(data = NA, 
                                         nrow = dim(daily.estim)[2], 
                                         ncol = dim(daily.estim)[3])

for (k1 in 1:dim(daily.estim)[2]){
  for (k2 in 1:dim(daily.estim)[3]){
    mean.mat[k1, k2] <- mean(daily.estim[,k1,k2])
    LCL.mat[k1, k2] <- quantile(daily.estim[,k1,k2], 0.025)
    UCL.mat[k1, k2] <- quantile(daily.estim[,k1,k2], 0.975)
  }
  
}

N.hats.day.Durban <- data.frame(Season = rep(seasons, each = dim(daily.estim)[2]),
                                Day = rep(1:dim(daily.estim)[2], length(seasons)),
                                Mean = as.vector(mean.mat),
                                LCL = as.vector(LCL.mat),
                                UCL = as.vector(UCL.mat))

# Daily estimates plots
p.daily.Durban <- ggplot(N.hats.day.Durban %>% group_by(Season)) + 
  geom_ribbon(aes(x = Day, ymin = LCL, ymax = UCL),
              fill = "blue", alpha = 0.5) +
  geom_path(aes(x = Day, y = Mean)) + 
  geom_point(data = Nhats.HT.all[Nhats.HT.all$Season %in% seasons,],
             aes(x = Day, y = Nhat),
             alpha = 0.3) + 
  facet_wrap(~ Season)

# Include non-survey years - no estimates for 2007/2008 because I don't have
# raw data for that year. Only the WinBUGS inputs. 
Laake.abundance.new <- read.csv(file = "Data/all_estimates_Laake_2024.csv") %>%
  mutate(LCL = CL.low,
         UCL = CL.high) %>%
  select(c(Season, Nhat, LCL, UCL)) %>%
  right_join(all.years, by = "Season") %>%
  arrange(year) %>%
  mutate(Method = "Laake")

# In reported estimates, there are two 2006/2007.
Reported.estimates %>%
  na.omit() %>%
  select(Season) %>% 
  unique() -> sampled.seasons 

Laake.abundance.new %>%
  rbind(Durban.abundance.df) %>%
  rbind(Nhat.) -> all.estimates

p.Nhats <- ggplot(all.estimates) +
  geom_point(aes(x = year, y = Nhat,
                 color = Method),
             alpha = 0.5) +
  geom_errorbar(aes(x = year, ymin = LCL, ymax = UCL,
                    color = Method)) +
  ylim(0, 45000)

#### Some observations in results over time. ###########################
# Precision of the estimates using Richards function is a lot better than other 
# methods. They are a bit lower (negatively biased) than other two. I think the
# problem is detection probabilities. mean.prob gets very small (mean = 0.04) 
# even though the prior is beta(0.95, 0.05). 

# Changed the prior to unif(0.9, 1.0) - 2024-10-21

# This problem has been fixed as of November 1, 2024. The problem was how I set
# up the proportion of watch period. It was set at hours observed over maximum
# possible. But, this should have been over 24 hrs (or one day), which was 
# already computed when data were extracted. By dividing that (times 60 min times
# 24 hrs) by 540 minutes made those numbers bigger by 2.667 times (24*60/540), 
# resulting in a smaller estimated abundance. 

# Although the estimates are a lot better now they are a bit off, especially in 
# some years. These need to be looked at carefully. 2024-11-13

####    ###########################

