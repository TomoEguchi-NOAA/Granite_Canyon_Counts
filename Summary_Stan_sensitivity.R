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
#  8. anchor_mu = qlogis(0.825)
#  9. use_shape_dev = 1 Periodic deviation from the curve is shared
#  #  
rm(list = ls())
library(tidyverse)
library(posterior)
library(cmdstanr)

source("Granite_Canyon_Counts_fcns.R")
source("ppc_richards_hssm.R")

# Sensitivity analyses are done only on one model (M1a2_1gamma)
# The default values for S1, S2, and likelihood in create.stan.data are set
# to run the M1a2 model as of 2026-07-29

sensitivity <- c("sen0", #  0. Parity to Jags. (sen0)
                 "sen1", #  1. use_pooling_Max = 0 
                 "sen2", #  2. anchor_mu = qlogis(0.70)
                 "sen3", #  3. anchor_mu = qlogis(0.90)
                 "sen4", #  4. use_trend_P = 0
                 "sen5", #  5. gamma_prior_mu = 0.0
                 "sen6", #  6. gamma_prior_sd = 2.0
                 "sen7", #  7. Base - all default values
                 "sen8", #  8. anchor_mu = qlogis(0.825)
                 "sen9") #  9. use_shape_dev = 1

# // ---- MODEL STRUCTURE (Table 1) ---------------------------------------
#   //   S1_by_season  S2_by_season  likelihood_NB      model
# //        1             1              0/1          M1a1 / M1a2
# //        0             0              0/1          M2a1 / M2a2
# //        1             0              0/1          M3a1 / M3a2
# //        0             1              0/1          M4a1 / M4a2

model <- "M1a2_1gamma"

S1_by_season <- S2_by_season <- likelihood_NB <- 1

# Create input data:
min.dur <- 60 

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
stan.data <- create.stan.data(jags.data = jags.data)

sensitivity.table <- data.frame(ID = paste0("sen", seq(0, 9)),
                                Sensitivity = c("Parity",
                                                "use_pooling_Max = 0",
                                                "anchor_mu = qlogis(0.7)",
                                                "anchor_mu = qlogis(0.9)",
                                                "use_trend_P = 0",
                                                "gamma_prior_mu = 0",
                                                "gamma_prios_sd = 2.0",
                                                "Base",
                                                "anchor_mu = qlogis(0.825)",
                                                "use_shape_dev = 1"),
                                Sens_abb = c("Parity",
                                             "PoolMax_0",
                                             "anchor_mu_7",
                                             "anchor_mu_9",
                                             "trendP_0",
                                             "gamma_mu_0",
                                             "gamma_sd_2",
                                             "Base",
                                             "anchor_mu_825",
                                             "use_shape_dev"),
                                ID.2 = c("Parity", "A", "B", "C", "D", "E", "F", "Base", "G", "H"))

peak.day <- gamma.hat <- conv.stats <- Nhats <- LOO.fit <- list()
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
  LOO.fit[[k]] <- fit_stan$loo()
}

## Look at LOOIC and auto correlation between the base and another model:
fit_base <- readRDS("RData//Richards_HSSM_M1a2_1gamma_mod3_stan_sen7.rds")
info_base <- readRDS("RData//Richards_HSSM_M1a2_1gamma_mod3_stan_sen7_info.rds")

M <- 9 # use_shape_dev = 1 to make the curve more flexible
fit_M  <- readRDS(paste0("RData//Richards_HSSM_M1a2_1gamma_mod3_stan_sen", M, ".rds"))
info_M  <- readRDS(paste0("RData//Richards_HSSM_M1a2_1gamma_mod3_stan_sen", M, "_info.rds"))

LOOIC_base <- LOO.fit[[7]]
LOOIC_M <- LOO.fit[[M]]

PPC_base <- ppc_setup(fit_base, info_base$stan.data, n_draws = 500, seed = 1)
AutoCorr_base <- ppc_autocorr(PPC_base, n_rep = 200)

PPC_M <- ppc_setup(fit_M, info_M$stan.data, n_draws = 500, seed = 1)
AutoCorr_M <- ppc_autocorr(PPC_M, n_rep = 200)

# σ_shape = 0.08 [0.05, 0.13] — small but not zero. On the log scale that's about 
# an 8% multiplicative deviation, so pentad-level departures from the Richards 
# curve of roughly ±8–17%. Real, but modest against counts that vary by orders 
# of magnitude across a season.
# 
# ΔLOOIC = 2.8, so elpd_diff ≈ 1.4 for 20 extra parameters. That's well below 
# the ~4× se_diff threshold — worth confirming with loo_compare, but at 1.4 elpd 
# the se_diff would have to be under 0.35 to matter, which won't happen. 
# No meaningful predictive gain.
# 
loo::loo_compare(LOOIC_base, LOOIC_M) -> LOO.comp
dv <- fit_base$sampler_diagnostics(format = "df")$divergent__
dr <- fit_base$draws(format = "df")
rbind(div = colMeans(dr[dv == 1, c("gamma_free[1]","sigma_proc_P","sigma_Obs")]),
      ok  = colMeans(dr[dv == 0, c("gamma_free[1]","sigma_proc_P","sigma_Obs")]))


# 
# Autocorrelation 0.198 → 0.182, an 8% reduction, still seven-plus SD outside 
# the replicate interval of [−0.072, 0.026]. The deviation absorbed a sliver 
# and left the phenomenon intact.
# 
# The residual correlation is not a curve-flexibility problem. You gave the model 
# 20 free parameters explicitly designed to depart from the Richards shape, and 
# it took almost none of them and barely moved the autocorrelation. Combined with 
# the flat loess — no mean-structure bias anywhere in the season — the conclusion 
# is that day-to-day passage rates are genuinely correlated, which no smooth 
# curve of any family would capture. Laake's year-specific splines show the same 
# scatter in his Figure 1.
# 
# Fewer effective observations than assumed means the posterior is more 
# concentrated than it should be. For 2025/2026 you report 16,032 [13,392, 19,475]; 
# with correlation accounted for that interval would be somewhat wider. The point 
# estimate wouldn't move much, since correlation affects precision rather than location.

# The crude AR(1) heuristic gives a variance inflation of (1+ρ)/(1−ρ) ≈ 1.48, 
# so an SE factor around 1.22 — roughly 20% wider. your abundance is a sum over 
# a fitted curve informed by every observation in the season, not a simple mean 
# of correlated draws, so the effective inflation is likely smaller.

find.large.pareto <- function(LOO.fit.list, stan.data){
  pareto.base <- as.data.frame(LOO.fit.list$pointwise)
  pareto.base %>%
    mutate(year.idx = stan.data$stan.data$year_idx,
           day.idx = stan.data$stan.data$day_idx,
           row.id = 1:length(stan.data$stan.data$n)) -> pareto.base
  
  p.all.pareto <- ggplot(pareto.base) +
    geom_point(aes(x = day.idx, y = influence_pareto_k)) +
    facet_wrap(~year.idx)
  
  stan.data.df <- data.frame(year.idx = stan.data$stan.data$year_idx,
                             day.idx = stan.data$stan.data$day_idx,
                             counts = stan.data$stan.data$n)
  
  p.all.counts <- ggplot(stan.data.df) +
    geom_point(aes(x = day.idx, y = counts)) +
    facet_wrap(~year.idx)

  large.pareto <- pareto.base %>% filter(influence_pareto_k > 0.7)
  large.pareto.years <- large.pareto$year.idx
  
  out.list <- list()
  for (k in 1:length(large.pareto.years)){
    # Year 24 (2007) has a large Pareto-k data point
    pareto.base %>% filter(year.idx == large.pareto.years[k]) -> pareto.base.Y
    stan.data.df %>% filter(year.idx == large.pareto.years[k]) -> counts.Y
    
    data.pareto.Y <- data.frame(day = counts.Y$day.idx,
                                counts = counts.Y$counts,
                                pareto = pareto.base.Y$influence_pareto_k)
    
    p.pareto.Y <- ggplot(data.pareto.Y) +
      geom_point(aes(x = day, y = counts, color = pareto))
    out.list[[k]] <- list(data = data.pareto.Y,
                          plot = p.pareto.Y)
  }
  
  return(list(out.list = out.list,
              large.pareto = large.pareto,
              p.all.pareto = p.all.pareto,
              p.all.counts = p.all.counts))
}

base.pareto <- find.large.pareto(LOO.fit[[7]], stan.data)
sen9.pareto <- find.large.pareto(LOO.fit[[10]], stan.data)

# Given one data point had a large pareto-k value (year = 24 (2007), day = 33), 
# remove it to see if the results
# change much. Remove the data point from the stan data
idx.data.1 <- c(1:(base.pareto$large.pareto$row.id-1), 
                (base.pareto$large.pareto$row.id+1):length(stan.data$stan.data$n))
stan.data.1 <- stan.data$stan.data
stan.data.1$N_flat <- length(idx.data.1)
stan.data.1$n <- stan.data.1$n[idx.data.1]
stan.data.1$bf <- stan.data.1$bf[idx.data.1]
stan.data.1$vs <- stan.data.1$vs[idx.data.1]
stan.data.1$observer_idx <- stan.data.1$observer_idx[idx.data.1]
stan.data.1$watch_length <- stan.data.1$watch_length[idx.data.1]
stan.data.1$day_idx <- stan.data.1$day_idx[idx.data.1]
stan.data.1$year_idx <- stan.data.1$year_idx[idx.data.1]

model <- "M1a2_1gamma"
out.file <- paste0("Richards_HSSM_", model, "_stan_No", base.pareto$large.pareto$row.id)
model.file <- file.path("models//model_Richards_HSSM_mod3.stan")

mod <- cmdstan_model(model.file, 
                     cpp_options = list(stan_threads = TRUE, 
                                        O = 3))
n_year <- stan.data$jags.data$n.year
n_observer <- stan.data$jags.data$n.obs.fixed

if (!file.exists(paste0("RData//", out.file, ".rds"))){
  fit_stan <- mod$sample(
    data            = stan.data.1,
    init            = stan_init_fn,
    chains          = 4,
    parallel_chains = 4,
    threads_per_chain = 2,
    iter_warmup     = 1500,
    iter_sampling   = 2000,
    adapt_delta     = 0.90  
  )
  
  fit_stan$save_object(file = paste0("RData//", out.file, ".rds"))

  info_stan <-  list(stan.data = stan.data$stan.data,
                     jags.data = stan.data$jags.data,
                     init_fn = stan_init_fn()) 
  
  saveRDS(info_stan,
          file = paste0("RData//", out.file, "_info.rds"))
  
} else {
  fit_stan <- readRDS(paste0("RData//", out.file, ".rds"))
}

LOO.out <- fit_stan$loo()
params.1.stan <- c("S1", "S2", "P", 
                   "sigma_proc_P", "Corrected_Est", "Max", "log_N_latent",
                   "gamma_free", "peak_day_decade", "beta1_P")

stan.global.summary <- fit_stan$summary(
  variables = params.1.stan,
  default_summary_measures(), 
  default_convergence_measures(),
  extra_quantiles = ~quantile2(., probs = c(0.025, 0.975)))

Nhats.1 <- fit_stan$summary("Corrected_Est")$mean

Nhats.df <- do.call(cbind, Nhats) %>% data.frame()
colnames(Nhats.df) <- sensitivity.table$Sens_abb

Nhats.df$No4074 <- Nhats.1
Nhats.df$year.idx <- c(1:nrow(Nhats.df))

ggplot(Nhats.df) +
  geom_point(aes(x = year.idx, y = Base), 
             color = "blue") +
  geom_point(aes(x = year.idx, y = No4074), 
             color = "gold", alpha = 0.5)

d.Y2007 <- Nhats.df[24, "Base"] - Nhats.df[24, "No4074"]

fit_stan$diagnostic_summary()

dv <- fit_stan$sampler_diagnostics(format = "df")$divergent__
dr <- fit_stan$draws(format = "df")
rbind(div = colMeans(dr[dv == 1, c("gamma_free[1]","sigma_proc_P","sigma_Obs")]),
      ok  = colMeans(dr[dv == 0, c("gamma_free[1]","sigma_proc_P","sigma_Obs")]))

### ### ### 


Laake.Run.Date <- "2026-02-23"
Laake.abundance.new <- read.csv(file = paste0("Data//all_estimates_Laake_2026_", 
                                              Laake.Run.Date, ".csv")) %>%
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
  geom_point(aes(x = start.year, y = width_95/mean))

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
