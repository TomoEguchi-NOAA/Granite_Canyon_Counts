# Summarizes the sensitivity analysis using Stan (Stan_sensitiviy.R)



#  0. Parity to Jags. (sen0)
#  1. use_pooling_Max = 0 — removes hierarchical shrinkage on the season levels
#                           altogether: each log_Max[y] gets an independent diffuse
#                           prior. Use this to test whether shrinkage toward the 
#                           trend inflates low seasons. When 0, beta0_Max / 
#                           beta1_Max / sigma_proc_Max are not identified and 
#                           simply sample their priors (sen1)
#  2. anchor_mu = qlogis(0.70) (sen2)
#  3. anchor_mu = qlogis(0.90) — detection sensitivity (sen3)
#  4. use_trend_P = 0 - is there a trend in P? (sen4)
#  5. gamma_prior_mu = 0.0
#  6. gamma_prior_sd = 2.0
#  7. Base - all default values
#  #  
rm(list = ls())
library(tidyverse)
library(posterior)
library(cmdstanr)

source("Granite_Canyon_Counts_fcns.R")

# Sensitivity analyses are done only on one model (M1a2_1gamma)
# The default values for S1, S2, and likelihood in create.stan.data are set
# to run the M1a2 model as of 2026-07-29

sensitivity <- c("sen0", "sen1", "sen2", "sen3", "sen4", "sen5", "sen6", "sen7", "sen8")

# // ---- MODEL STRUCTURE (Table 1) ---------------------------------------
#   //   S1_by_season  S2_by_season  likelihood_NB      model
# //        1             1              0/1          M1a1 / M1a2
# //        0             0              0/1          M2a1 / M2a2
# //        1             0              0/1          M3a1 / M3a2
# //        0             1              0/1          M4a1 / M4a2

model <- "M1a2_1gamma"

S1_by_season <- S2_by_season <- likelihood_NB <- 1

# Create input data:
min.dur <- 60 #10 #85 #

# These are the ending year of each season - for example, 2022 in the following vector indicates
# for the 2021/2022 season. These data were extracted using Extract_Data_All_v2.Rmd
# Data prior to the 2009/2010 season are in Laake's ERAnalayis package. 
years <- c(2008, 2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024, 2025, 2026)
data.dir <- "RData/V2.1_May2026"
max.day <- 100

#jags.input.list <- AllData2JagsInput_NoBUGS(min.dur, years = years, data.dir, max.day)         

jags.input <- NoBUGS_Jags_input(min.dur = min.dur, 
                                years = years, 
                                data.dir = data.dir, 
                                max.day = max.day, 
                                obs.n.min = 10, N.obs = 10)

jags.data <- jags.input$jags.data


sensitivity.table <- data.frame(ID = paste0("sen", seq(0, 8)),
                                Sensitivity = c("Parity",
                                                "use_pooling_Max = 0",
                                                "anchor_mu = qlogis(0.7)",
                                                "anchor_mu = qlogis(0.9)",
                                                "use_trend_P = 0",
                                                "gamma_prior_mu = 0",
                                                "gamma_prios_sd = 2.0",
                                                "Base",
                                                "anchor_mu = qlogis(0.825)"),
                                Sens_abb = c("Parity",
                                             "PoolMax_0",
                                             "anchor_mu_7",
                                             "anchor_mu_9",
                                             "trendP_0",
                                             "gamma_mu_0",
                                             "gamma_sd_2",
                                             "Base",
                                             "anchor_mu_825"),
                                ID.2 = c("Parity", "A", "B", "C", "D", "E", "F", "G", "Base"))

peak.day <- gamma.hat <- conv.stats <- Nhats <- list()
for (k in 1:length(sensitivity)){
  out.file <- paste0("Richards_HSSM_", model, "_mod3_stan_", sensitivity[k])
  fit_stan <- readRDS(paste0("RData//", out.file, ".rds"))
  
  Nhats[[k]] <- fit_stan$summary("Corrected_Est")$mean
  if (k > 1){
    conv.stats[[k]] <- fit_stan$summary(
      variables = c("beta0_Max", "sigma_proc_Max", "Corrected_Est", "P", "S1", "S2", "gamma_free"),
      default_summary_measures(), 
      default_convergence_measures(),
      extra_quantiles = ~quantile2(., probs = c(0.025, 0.975)))
    gamma.hat[[k]] <- fit_stan$summary("gamma_free")
    
  } else {
    conv.stats[[k]] <- fit_stan$summary(
      variables = c("beta0_Max", "sigma_proc_Max", "Corrected_Est", "P", "S1", "S2"),
      default_summary_measures(), 
      default_convergence_measures(),
      extra_quantiles = ~quantile2(., probs = c(0.025, 0.975)))
    gamma.hat[[k]] <- NULL
    
  }
  peak.day[[k]] <- fit_stan$summary("peak_day_decade")
}

Nhats.df <- do.call(cbind, Nhats) %>% data.frame()
colnames(Nhats.df) <- sensitivity.table$Sens_abb

Laake.Run.Date <- "2026-02-23"
Laake.abundance.new <- read.csv(file = paste0("Data//all_estimates_Laake_2026_", Laake.Run.Date, ".csv")) %>%
  mutate(LCL = CL.low,
         UCL = CL.high) %>%
  na.omit()

Nhats.df$Laake <- Laake.abundance.new$Nhat
start.years <- Laake.abundance.new$Start.Year
Nhats.df$Start.year <- start.years

#write.csv(Nhats.df, file = "Data/Sensitivity_Nhats_M1a2_1gamma.csv")

conv.stats.df <- do.call(rbind, conv.stats) %>% as.data.frame()
conv.stats.df$sens <- rep(sensitivity.table$Sens_abb, 
                          times = c(lapply(conv.stats, FUN = function(x) nrow(x))%>%unlist()))

#write.csv(conv.stats.df, file = "Data/Sensitivity_conv_M1a2_1gamma.csv")
gamma.free <- do.call(rbind, gamma.hat) %>% 
  data.frame() %>%
  mutate(Model = sensitivity.table$Sens_abb[2:nrow(sensitivity.table)])
#write.csv(gamma.free, file = "Data/Sensitivity_gamma_M1a2_1gamma.csv")

conv.stats.df %>%
  filter(sens == "Base") -> base.conv.stats.df

Nhats.base <- base.conv.stats.df[grep("^Corrected", base.conv.stats.df$variable),]

Nhats.base %>%
  mutate(width_95 = q97.5 - q2.5,
         start.year = start.years) -> Nhats.base
#write.csv(Nhats.base, file = "Data/Nhats_M1a2_1gamma.csv")

ggplot(Nhats.base) +
  geom_point(aes(x = start.year, y = width_95))

# How many of Laake's estimates were within 95% CI of mine with meanP = 0.825?

conv.stats.df %>%
  filter(sens == "anchor_mu_825") -> tmp
Nhats.anchor_825 <- tmp[grep("^Corrected", tmp$variable),]

Nhats.anchor_825 %>%
  mutate(Laake = Laake.abundance.new$Nhat,
         Start.year = start.years) %>%
  mutate(Laake_in = ifelse((Laake > q2.5 & Laake < q97.5), 1, 0)) -> Nhats.anchor_825

#write.csv(Nhats.anchor_825, file = "Data/Sensitivity_anchor_825_M1a2_1gamma.csv")

ggplot(Nhats.anchor_825) +
  geom_point(aes(x = Start.year, y = mean)) +
  geom_errorbar(aes(x = Start.year, ymin = q2.5, ymax = q97.5)) +
  geom_point(aes(x = Start.year, y = Laake, color = as.factor(Laake_in))) +
  theme(legend.location = "top")
