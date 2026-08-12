# Sensitivity analysis of Stan models as per Claude's suggestions

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
#  10. use_shape_dev = 1 AND n_period = 10. Periodic grouping is changed from 5 to 10 days

rm(list = ls())
library(tidyverse)
library(posterior)
library(cmdstanr)

source("Granite_Canyon_Counts_fcns.R")

# Sensitivity analyses are done only on one model (M1a2_1gamma)
# The default values for S1, S2, and likelihood in create.stan.data are set
# to run the M1a2 model as of 2026-07-29
sensitivity <- "sen9"

# // ---- MODEL STRUCTURE (Table 1) ---------------------------------------
#   //   S1_by_season  S2_by_season  likelihood_NB      model
# //        1             1              0/1          M1a1 / M1a2
# //        0             0              0/1          M2a1 / M2a2
# //        1             0              0/1          M3a1 / M3a2
# //        0             1              0/1          M4a1 / M4a2

model <- "M1a2_1gamma"

S1_by_season <- S2_by_season <- likelihood_NB <- 1

if (length(grep("a1", model)) == 0) likelihood_NB <- 0
if (length(grep("M2", model)) == 1) {
  S1_by_season <- 0
  S2_by_season <- 0
}

if (length(grep("M3", model)) == 1) {
  S2_by_season <- 0
}

if (length(grep("M4", model)) == 1) {
  S1_by_season <- 0
}

# Create input data:
min.dur <- 60 #10 #85 #
Run.date <- Sys.Date()

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

sensitivity.table <- data.frame(ID = paste0("sen", seq(0, 10)),
                                Sensitivity = c("Parity",
                                                "use_pooling_Max = 0",
                                                "anchor_mu = qlogis(0.7)",
                                                "anchor_mu = qlogis(0.9)",
                                                "use_trend_P = 0",
                                                "gamma_prior_mu = 0",
                                                "gamma_prios_sd = 2.0",
                                                "Base",
                                                "anchor_mu = qlogis(0.825)",
                                                "use_shape_dev = 1",
                                                "n_period = 10"),
                                ID.2 = c("Parity", "A", "B", "C", "D", "E", "F", "Base", "G", "H", "I"))

if (sensitivity == "sen1"){
  stan.data <- create.stan.data(jags.data = jags.data, 
                                use_pooling_Max = 0)
} else if (sensitivity == "sen2"){
  stan.data <- create.stan.data(jags.data = jags.data, 
                                anchor_mu = qlogis(0.7))
} else if (sensitivity == "sen3"){
  stan.data <- create.stan.data(jags.data = jags.data, 
                                anchor_mu = qlogis(0.9))
} else if (sensitivity == "sen4"){
  stan.data <- create.stan.data(jags.data = jags.data, 
                                use_trend_P = 0)
} else if (sensitivity == "sen0"){
  stan.data <- create.stan.data(jags.data = jags.data, 
                                anchor_sd = 0.01, estimate_gamma = 0)
} else if (sensitivity == "sen5"){
  stan.data <- create.stan.data(jags.data = jags.data, 
                                gamma_prior_mu = 0)
} else if (sensitivity == "sen6"){
  stan.data <- create.stan.data(jags.data = jags.data, 
                                gamma_prior_sd = 2.0)
} else if (sensitivity == "sen7"){
  stan.data <- create.stan.data(jags.data = jags.data)
} else if (sensitivity == "sen8"){
  stan.data <- create.stan.data(jags.data = jags.data, 
                                anchor_mu = qlogis(0.825))
} else if (sensitivity == "sen9"){
  stan.data <- create.stan.data(jags.data = jags.data, 
                                use_shape_dev = 1)
} else if (sensitivity == "sen10"){
  stan.data <- create.stan.data(jags.data = jags.data, 
                                use_shape_dev = 1,
                                n_period = 10)
}

# Create an inits function
n_year <- stan.data$jags.data$n.year
n_observer <- stan.data$jags.data$n.obs.fixed

file <- file.path("models//model_Richards_HSSM_mod4.stan")
out.file <- paste0("Richards_HSSM_", model, "_mod4_", sensitivity, "_stan")
# 
# # --- 5. Inspect Results ---
# params.1.stan <- c("S1", "beta_S1", "S2", "P", "phi",
#                    "sigma_proc_P", "Corrected_Est", "Max", "gamma_free[1]")
# 
# # Jags output:
# jm.out <- readRDS("RData/JAGS_Richards_HSSM_M1a2_1968to2026_min60_2026-07-27_NoBUGS.rds")
# 
# # From Jags using regular expression syntax:
# params.1.jags <- c("^S1\\[", "^S1.beta", "^S2\\[", "^P\\[", "^r",
#                    "^sd.proc.P", "^Corrected.Est", "^Max\\[")
# 
# jags.params.idx <- lapply(params.1.jags,
#                           FUN = function(x) grep(x, jm.out$posterior.summary$variable)) %>% unlist()
# 
# summarise_draws(jm.out$jm$samples, 
#                 default_summary_measures(),
#                 default_convergence_measures(),
#                 extra_quantiles = ~quantile2(., probs = c(0.025, 0.975))) %>%
#   as.data.frame() -> jags.convergence.summary
# 
# jags.global.summary <- jags.convergence.summary[jags.params.idx,]
# 
# jags.Nhats <- jags.global.summary[grep("^Corrected.Est", jags.global.summary$variable),] %>%
#   mutate(MCMC = "Jags",
#          year.idx = seq(1, n_year))
# jags.Nhats$start.year <- jags.data$start.years
# 
# # Use samples to compute SE of Nhats
# Nhats_jags <- jm.out$jm$sims.list$Corrected.Est      # or S2
# mean.Nhats_jags <- apply(Nhats_jags, 2, mean) 
# se.Nhats_jags <- apply(Nhats_jags, 2, sd)/sqrt(coda::effectiveSize(Nhats_jags))


#mod <- cmdstan_model(file)
if (!file.exists(paste0("RData//", out.file, ".rds"))){
  # Compile with aggressive C++ optimization flags
  mod <- cmdstan_model(file, 
                       cpp_options = list(stan_threads = TRUE, 
                                          O = 3))
  fit_stan <- mod$sample(
    data            = stan.data$stan.data,
    init            = stan_init_fn,
    chains          = 4,
    parallel_chains = 4,
    threads_per_chain = 2,
    iter_warmup     = 1500,
    iter_sampling   = 2000,
    adapt_delta     = 0.90  
  )
  
  fit_stan$save_object(file = paste0("RData//", out.file, ".rds"))
  saveRDS(list(stan.data = stan.data$stan.data,
               jags.data = stan.data$jags.data,
               init_fn = stan_init_fn()),
          file = paste0("RData//", out.file, "_info.rds"))
  
} else {
  fit_stan <- readRDS(paste0("RData//", out.file, ".rds"))
}
# 
# # --- Get Summaries for Specific Global Parameters ---
# stan_global_summary <- fit_stan$summary(
#   variables = params.1.stan,
#   default_summary_measures(), 
#   default_convergence_measures(),
#   extra_quantiles = ~quantile2(., probs = c(0.025, 0.975)))
# 
# stan.Nhats <- stan_global_summary[grep("^Corrected_Est", stan_global_summary$variable),] %>%
#   mutate(MCMC = "Stan",
#          year.idx = seq(1, n_year)) 
# stan.Nhats$start.year <- jags.data$start.years
# 
# Nhats.stan.jags <- rbind(stan.Nhats, jags.Nhats)
# 
# ggplot(Nhats.stan.jags) +
#   geom_point(aes(x = year.idx, y = mean, color = MCMC),
#              alpha = 0.6) +
#   geom_errorbar(aes(x = year.idx, ymin = q2.5, ymax = q97.5)) +
#   theme(legend.position = "top")
# 
# Nhats_stan_summary  <- fit_stan$summary("Corrected_Est")
# mean.Nhats_stan <- Nhats_stan_summary$mean
# se.Nhats_stan <- Nhats_stan_summary$sd / sqrt(Nhats_stan_summary$ess_bulk)
# 
# z <- (mean.Nhats_stan - mean.Nhats_jags) / sqrt(se.Nhats_jags^2 + se.Nhats_stan^2)     # standardised difference
# rel <- (mean.Nhats_stan - mean.Nhats_jags) / mean.Nhats_jags
# 
# round(summary(z), 2)
# sum(abs(z) > 2)                # how many exceed 2 MCSE
# 
# #systematic bias
# mean.Nhats.df <- data.frame(Jags = mean.Nhats_jags,
#                             Stan = mean.Nhats_stan,
#                             start.year = jags.data$start.years)
# 
# title.txt <- paste0(sensitivity, "(", 
#                     sensitivity.table[grep(sensitivity, sensitivity.table$ID), "ID.2"],  "): ",
#                     sensitivity.table[grep(sensitivity, sensitivity.table$ID), "Sensitivity"])
# 
# ggplot(mean.Nhats.df) +
#   geom_point(aes(x = Jags, y = Stan, color = start.year)) +
#   geom_abline(slope = 1, color = "red") +
#   labs(title = title.txt) +
#   xlab("Jags Mean Corrected Est") +
#   ylab("Stan Mean Corrected Est") 
# 
# sum(mean.Nhats_stan > mean.Nhats_jags)                                     # want ~17, not 34 or 0
# binom.test(sum(mean.Nhats_stan > mean.Nhats_jags), length(mean.Nhats_stan))$p.value    # want non-significant
# summary((mean.Nhats_stan - mean.Nhats_jags) / sqrt(se.Nhats_jags^2 + se.Nhats_stan^2))       # want |z| mostly < 2

# If there were divergence - see below:
# dv <- fit_stan$sampler_diagnostics(format = "df")$divergent__
# dr <- fit_stan$draws(format = "df")
# colMeans(dr[dv == 1, c("sigma_proc_P","sigma_proc_Max","sigma_Obs","phi")])
# colMeans(dr[dv == 0, c("sigma_proc_P","sigma_proc_Max","sigma_Obs","phi")])

# colMeans(dr[dv == 1, c("delta[1]","sigma_proc_P","sigma_proc_Max","sigma_Obs","phi")])
# colMeans(dr[dv == 0, c("delta[1]","sigma_proc_P","sigma_proc_Max","sigma_Obs","phi")])
# 
# fit_stan$summary(c("peak_day_decade", "peak_day_slope"))
# fit_stan$summary("peak_day")
# fit_stan$summary("peak_width")
# # 
# pd <- fit_stan$draws("peak_day", format = "matrix")
# P  <- fit_stan$draws("P",        format = "matrix")
# off <- pd - P                                     # mode minus midpoint, per season
# sl <- apply(off, 1, \(o) coef(lm(o ~ stan.data$stan.data$year_values))[2])
# quantile(sl * 10, c(0.025, 0.5, 0.975))           # days per decade of drift
