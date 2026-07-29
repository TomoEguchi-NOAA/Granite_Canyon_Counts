# Running all models in STAN
# 
#  Collects Stan runs and compares results. Select the best model and do 
#  posterior predictive checks. 
#  
#  Run all models in Run_StanRichards_AllModels.R. It will save results
#  in .rds files. 
# 

rm(list = ls())
library(tidyverse)
library(posterior)
library(bayesplot)

source("Granite_Canyon_Counts_fcns.R")
source("ppc_richards_hssm.R")

min.dur <- 60 #10 #85 #
Run.date <- Sys.Date()

# These are the ending year of each season - for example, 2022 in the following vector indicates
# for the 2021/2022 season. These data were extracted using Extract_Data_All_v2.Rmd
# Data prior to the 2009/2010 season are in Laake's ERAnalayis package. 
years <- c(2008, 2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024, 2025, 2026)
data.dir <- "RData/V2.1_May2026"
max.day <- 100


jags.input <- NoBUGS_Jags_input(min.dur = min.dur, 
                                years = years, 
                                data.dir = data.dir, 
                                max.day = max.day, 
                                obs.n.min = 10, N.obs = 10)

jags.data <- jags.input$jags.data

stan.data <- create.stan.data(jags.data = jags.data)

# // ---- MODEL STRUCTURE (Table 1) ---------------------------------------
#   //   S1_by_season  S2_by_season  likelihood_NB      model
# //        1             1              0/1          M1a1 / M1a2
# //        0             0              0/1          M2a1 / M2a2
# //        1             0              0/1          M3a1 / M3a2
# //        0             1              0/1          M4a1 / M4a2
models <- c("M1a1", "M2a1", "M3a1", "M4a1",
            "M1a2", "M2a2", "M3a2", "M4a2")

params.1.stan <- c("S1", "beta_S1", "S2", "P", "phi", "alpha",
                   "sigma_proc_P", "Corrected_Est", "Max", "log_N_latent")

LOO.out <- stan.global.summary <- ppc.res <- list()
k <- 1
for (k in 1:length(models)) {
  out.file <- paste0("Richards_HSSM_", models[k], "_stan")
  
  #mod <- cmdstan_model(file)
  fit_stan <- readRDS(paste0("RData/", out.file, ".rds"))
  
  LOO.out[[k]] <- fit_stan$loo()
  stan.global.summary[[k]] <- fit_stan$summary(
    variables = params.1.stan,
    default_summary_measures(), 
    default_convergence_measures(),
    extra_quantiles = ~quantile2(., probs = c(0.025, 0.975)))
  
  # --- Posterior Predictive Simulation Loop ---
  # res$autocorr is the key one. Within-season lag-1 autocorrelation of daily-averaged residuals, compared against the same statistic computed on replicate datasets. If the Richards curve were too rigid to track the migration, residuals would come in same-signed runs along the season and the observed autocorrelation would exceed the replicate distribution. If the observed value sits inside the replicate interval, you have direct evidence the curve is flexible enough — an empirical answer where you currently have an argument.
  
  # plots$resid_day is the same thing visually, and it's the figure I'd put in the paper. A flat loess through residuals against day-of-season says the curve captures the shape; systematic curvature — a dip at the peak, say — would say it doesn't.

# pvalues covers proportion of zeros (the check from your calf memo), maximum, SD, mean, and the upper tail. The flag column marks anything outside [0.025, 0.975].

# plots$coverage gives the empirical coverage of 95% predictive intervals. This is the cleanest single number for the Results — if it lands near 95%, the observation model is calibrated. It also detects Poisson vs NB automatically by looking for phi[1], so you can run it unchanged across all eight models. Worth doing: comparing coverage and the zeros p-value between M1a1 and M1a2 is a much more direct justification for the negative binomial than ΔLOOIC is.

# plots$resid_season faceted QQ plots will show whether any particular season fits badly. Worth checking against the ones with sparse effort, and especially 2025/2026.
 
   ppc.res[[k]] <- ppc_richards(fit_stan, stan.data$stan.data, n_draws = 500)
  
  # res$pvalues
  # res$autocorr
  # res$plots$resid_day
}
  
# --- Model comparison ---





#jm.out <- readRDS()
# The following plot can be useful: Set it aside for the manuscript. Predicted-versus-observed counts across all watch periods and all 34 seasons, correlation 0.88, with the scatter fanning out proportionally to the mean exactly as a negative binomial should — that is a direct, visual answer to reviewer #2. It shows the Richards curve tracking the observed counts across the full range without any spline flexibility, and the mean-variance relationship supporting the NB choice.

# It's also most of the posterior predictive check your Methods promise and never deliver. Regenerate it from posterior draws rather than plug-in medians, add the 95% predictive band, and colour by season — that plus the within-season residual autocorrelation check gives you the whole GOF section.
# 
# reconstruct kappa at the JAGS posterior medians, using the STAN data list
# Pm <- apply(jm.out$jm$sims.list$P, 2, median)
# Mm <- apply(jm.out$jm$sims.list$Max, 2, median)
# S1m <- apply(jm.out$jm$sims.list$S1, 2, median)
# S2m <- apply(jm.out$jm$sims.list$S2, 2, median)
# g <- exp(1)
# pmt <- Pm[stan_data$year_idx] - stan_data$day_idx
# C1 <- (1 + (2*g-1)*exp(-pmt/S1m[stan_data$year_idx]))^(-1/g)
# C2 <- (1 + (2*g-1)*exp( pmt/S2m[stan_data$year_idx]))^(-1/g)
# pred <- Mm[stan_data$year_idx]*C1*C2 * stan_data$watch_length * 0.8
# plot(pred, stan_data$n); abline(0,1,col="red")
# cor(pred, stan_data$n)
# params.1.jags <- c("^S1\\[", "^S1.beta", "^S2\\[", "^P\\[", "^r",
#                    "^sd.proc.P", "^Corrected.Est", "^Max\\[")
# 
# jags.params.idx <- lapply(params.1.jags,
#                      FUN = function(x) grep(x, jm.out$posterior.summary$variable)) %>% unlist()
# 
# summarise_draws(jm.out$jm$samples, 
#                 default_summary_measures(),
#                 default_convergence_measures(),
#                 extra_quantiles = ~quantile2(., probs = c(0.025, 0.975))) %>%
#   as.data.frame() -> jags.convergence.summary
# 
# jags.global.summary <- jags.convergence.summary[jags.params.idx,]
# 
# stan.S1 <- stan_global_summary[grep("^S1\\[", stan_global_summary$variable),] %>%
#   mutate(MCMC = "Stan",
#          year.idx = seq(1, n_year))
# jags.S1 <- jags.global.summary[grep("^S1\\[", jags.global.summary$variable),] %>%
#   mutate(MCMC = "Jags",
#          year.idx = seq(1, n_year))
# 
# S1.stan.jags <- rbind(stan.S1, jags.S1)
# 
# ggplot(S1.stan.jags) +
#   geom_point(aes(x = year.idx, y = mean, color = MCMC)) +
#   geom_errorbar(aes(x = year.idx, ymin = q2.5, ymax = q97.5)) +
#   theme(legend.position = "top")
# 
# s_jags <- jm.out$jm$sims.list$S1      # or S2
# m_j <- apply(s_jags, 2, mean); se_j <- apply(s_jags, 2, sd)/sqrt(coda::effectiveSize(s_jags))
# sj  <- fit_stan$summary("S1")
# m_s <- sj$mean; se_s <- sj$sd / sqrt(sj$ess_bulk)
# 
# z <- (m_s - m_j) / sqrt(se_j^2 + se_s^2)     # standardised difference
# rel <- (m_s - m_j) / m_j
# 
# round(summary(z), 2)
# sum(abs(z) > 2)                # how many exceed 2 MCSE
# binom.test(sum(m_s > m_j), length(m_s))$p.value   # sign test for systematic bias
# plot(m_j, m_s); abline(0, 1, col = "red")
# 
# sum(m_s > m_j)                                     # want ~17, not 34 or 0
# binom.test(sum(m_s > m_j), length(m_s))$p.value    # want non-significant
# summary((m_s - m_j) / sqrt(se_j^2 + se_s^2))       # want |z| mostly < 2
# 
# 
# stan.Nhats <- stan_global_summary[grep("^Corrected_Est", stan_global_summary$variable),] %>%
#   mutate(MCMC = "Stan",
#          year.idx = seq(1, n_year))
# jags.Nhats <- jags.global.summary[grep("^Corrected.Est", jags.global.summary$variable),] %>%
#   mutate(MCMC = "Jags",
#          year.idx = seq(1, n_year))
# 
# Nhats.stan.jags <- rbind(stan.Nhats, jags.Nhats)
# 
# ggplot(Nhats.stan.jags) +
#   geom_point(aes(x = year.idx, y = mean, color = MCMC)) +
#   geom_errorbar(aes(x = year.idx, ymin = q2.5, ymax = q97.5)) +
#   theme(legend.position = "top")
# 
# s_jags <- jm.out$jm$sims.list$Corrected.Est      # or S2
# m_j <- apply(s_jags, 2, mean); se_j <- apply(s_jags, 2, sd)/sqrt(coda::effectiveSize(s_jags))
# sj  <- fit_stan$summary("Corrected_Est")
# m_s <- sj$mean; se_s <- sj$sd / sqrt(sj$ess_bulk)
# 
# z <- (m_s - m_j) / sqrt(se_j^2 + se_s^2)     # standardised difference
# rel <- (m_s - m_j) / m_j
# 
# round(summary(z), 2)
# sum(abs(z) > 2)                # how many exceed 2 MCSE
# #systematic bias
# plot(m_j, m_s); abline(0, 1, col = "red")
# 
# sum(m_s > m_j)                                     # want ~17, not 34 or 0
# binom.test(sum(m_s > m_j), length(m_s))$p.value    # want non-significant
# summary((m_s - m_j) / sqrt(se_j^2 + se_s^2))       # want |z| mostly < 2
# 
# fit_stan$summary("logit_p0")   # expect mean ~1.40, above the 1.39 anchor
# 
# # same comparison on Raw_Est, which has no correction factor
# s_jags <- jm.out$jm$sims.list$Raw.Est      # or S2
# m_j <- apply(s_jags, 2, mean); se_j <- apply(s_jags, 2, sd)/sqrt(coda::effectiveSize(s_jags))
# sj  <- fit_stan$summary("Raw_Est")
# m_s <- sj$mean; se_s <- sj$sd / sqrt(sj$ess_bulk)
# 
# z <- (m_s - m_j) / sqrt(se_j^2 + se_s^2)     # standardised difference
# rel <- (m_s - m_j) / m_j
# 
# round(summary(z), 2)
# sum(abs(z) > 2)                # how many exceed 2 MCSE
# #systematic bias
# plot(m_j, m_s); abline(0, 1, col = "red")
# 
# sum(m_s > m_j)                                     # want ~17, not 34 or 0
# binom.test(sum(m_s > m_j), length(m_s))$p.value    # want non-significant
# summary((m_s - m_j) / sqrt(se_j^2 + se_s^2))       # want |z| mostly < 2
# 
# mean(fit_stan$summary("Raw_Est")$mean) / mean(apply(jm.out$jm$sims.list$Raw.Est, 2, mean))
# sum(m_s > m_j)
# 
# # and check the two correction factors directly
# fit_stan$summary("corr_factor")$mean
# 
# fit_stan$summary("sigma_Obs")$mean
# median(jm.out$jm$sims.list$sd.obs)     # was 0.243
# 
# 
# # implied mean detection in each fit
# mean(plogis(jm.out$jm$sims.list$alpha))           # JAGS, before covariates
# fit_stan$summary(c("BF_Fixed","VS_Fixed","logit_p0","sigma_Obs"))$mean
# 
# # vs JAGS BF.Fixed / VS.Fixed
# summary(stan_data$bf); summary(stan_data$vs)   # raw ranges? vs shifted?
# 
# # --- Get Summaries for Year-Specific Parameters ---
# # This will print the estimates for every year automatically
# richards_summary <- fit_stan$summary(
#   variables = c("Max", "S1", "S2", "K", "P1"),
#   "mean", "median", "quantile~0.025", "quantile~0.975", "rhat"
# )
# print(richards_summary)
# 
# library(bayesplot)
# library(ggplot2)
# 
# # Extract draws for the parameters you want to check
# # (Using a few years of Max and S1 as an example)
# draws_diagnostic <- fit_stan$draws(c("mean_prob", "Max[1]", "Max[10]", "S1[1]"))
# 
# # Generate trace plots
# mcmc_trace(draws_diagnostic) +
#   theme_minimal() +
#   labs(title = "MCMC Chain Trace Plots")
# 
# # Extract all S1 year estimates
# s1_draws <- fit_stan$draws("S1")
# 
# # Plot the intervals across all years sequentially
# mcmc_intervals(s1_draws) +
#   theme_minimal() +
#   labs(title = "Posterior Estimates for S1 Across Years",
#        x = "Parameter Value", y = "Year Index")
# 
# # Compare the fixed effect coefficients
# covariate_draws <- fit_stan$draws(c("BF_Fixed", "VS_Fixed"))
# 
# mcmc_areas(covariate_draws, prob = 0.95) +
#   theme_minimal() +
#   labs(title = "Posterior Distributions of Detection Covariates",
#        x = "Effect Size")
# 
# # 1. Extract the mean_N matrix from the posterior
# # format = "draws_matrix" converts it into a clean math grid
# mean_N_posterior <- fit_stan$draws("mean_N", format = "draws_matrix")
# 
# # 2. Pick a specific year to visualize (e.g., Year 5)
# target_year <- 5
# day_columns <- paste0("mean_N[", 1:n.days, ",", target_year, "]")
# 
# # Extract only the columns belonging to that year
# year_data <- mean_N_posterior[, day_columns]
# 
# # 3. Calculate the mean and 95% Credible Interval for each day
# trajectory <- data.frame(
#   Day = 1:n.days,
#   Mean = colMeans(year_data),
#   Lower = apply(year_data, 2, quantile, probs = 0.025),
#   Upper = apply(year_data, 2, quantile, probs = 0.975)
# )
# 
# # 4. Plot the migration curve!
# ggplot(trajectory, aes(x = Day, y = Mean)) +
#   geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "blue", alpha = 0.15) +
#   geom_line(color = "blue", size = 1) +
#   theme_minimal() +
#   labs(title = paste("Estimated Gray Whale Migration Curve: Year", target_year),
#        x = "Days into Season",
#        y = "Estimated Expected Abundance (mean_N)")



##############################################
# 
# # --- 5. Inspect Results ---
# params.1.stan <- c("S1", "beta_S1", "S2", "P", "phi",
#                    "sigma_proc_P", "Corrected_Est", "Max")
# 
# # --- Get Summaries for Specific Global Parameters ---
# stan_global_summary <- fit_stan$summary(
#   variables = params.1.stan,
#   default_summary_measures(), 
#   default_convergence_measures(),
#   extra_quantiles = ~quantile2(., probs = c(0.025, 0.975)))
# 
# print(stan_global_summary)
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
# stan.S1 <- stan_global_summary[grep("^S1\\[", stan_global_summary$variable),] %>%
#   mutate(MCMC = "Stan",
#          year.idx = seq(1, n_year))
# jags.S1 <- jags.global.summary[grep("^S1\\[", jags.global.summary$variable),] %>%
#   mutate(MCMC = "Jags",
#          year.idx = seq(1, n_year))
# 
# S1.stan.jags <- rbind(stan.S1, jags.S1)
# 
# ggplot(S1.stan.jags) +
#   geom_point(aes(x = year.idx, y = mean, color = MCMC)) +
#   geom_errorbar(aes(x = year.idx, ymin = q2.5, ymax = q97.5)) +
#   theme(legend.position = "top")
# 
# s_jags <- jm.out$jm$sims.list$S1      # or S2
# m_j <- apply(s_jags, 2, mean); se_j <- apply(s_jags, 2, sd)/sqrt(coda::effectiveSize(s_jags))
# sj  <- fit_stan$summary("S1")
# m_s <- sj$mean; se_s <- sj$sd / sqrt(sj$ess_bulk)
# 
# z <- (m_s - m_j) / sqrt(se_j^2 + se_s^2)     # standardised difference
# rel <- (m_s - m_j) / m_j
# 
# round(summary(z), 2)
# sum(abs(z) > 2)                # how many exceed 2 MCSE
# binom.test(sum(m_s > m_j), length(m_s))$p.value   # sign test for systematic bias
# plot(m_j, m_s); abline(0, 1, col = "red")
# 
# sum(m_s > m_j)                                     # want ~17, not 34 or 0
# binom.test(sum(m_s > m_j), length(m_s))$p.value    # want non-significant
# summary((m_s - m_j) / sqrt(se_j^2 + se_s^2))       # want |z| mostly < 2
# 
# 
# stan.Nhats <- stan_global_summary[grep("^Corrected_Est", stan_global_summary$variable),] %>%
#   mutate(MCMC = "Stan",
#          year.idx = seq(1, n_year))
# jags.Nhats <- jags.global.summary[grep("^Corrected.Est", jags.global.summary$variable),] %>%
#   mutate(MCMC = "Jags",
#          year.idx = seq(1, n_year))
# 
# Nhats.stan.jags <- rbind(stan.Nhats, jags.Nhats)
# 
# ggplot(Nhats.stan.jags) +
#   geom_point(aes(x = year.idx, y = mean, color = MCMC)) +
#   geom_errorbar(aes(x = year.idx, ymin = q2.5, ymax = q97.5)) +
#   theme(legend.position = "top")
# 
# s_jags <- jm.out$jm$sims.list$Corrected.Est      # or S2
# m_j <- apply(s_jags, 2, mean); se_j <- apply(s_jags, 2, sd)/sqrt(coda::effectiveSize(s_jags))
# sj  <- fit_stan$summary("Corrected_Est")
# m_s <- sj$mean; se_s <- sj$sd / sqrt(sj$ess_bulk)
# 
# z <- (m_s - m_j) / sqrt(se_j^2 + se_s^2)     # standardised difference
# rel <- (m_s - m_j) / m_j
# 
# round(summary(z), 2)
# sum(abs(z) > 2)                # how many exceed 2 MCSE
# #systematic bias
# plot(m_j, m_s); abline(0, 1, col = "red")
# 
# sum(m_s > m_j)                                     # want ~17, not 34 or 0
# binom.test(sum(m_s > m_j), length(m_s))$p.value    # want non-significant
# summary((m_s - m_j) / sqrt(se_j^2 + se_s^2))       # want |z| mostly < 2
# 
# fit_stan$summary("logit_p0")   # expect mean ~1.40, above the 1.39 anchor
# 
# # same comparison on Raw_Est, which has no correction factor
# s_jags <- jm.out$jm$sims.list$Raw.Est      # or S2
# m_j <- apply(s_jags, 2, mean); se_j <- apply(s_jags, 2, sd)/sqrt(coda::effectiveSize(s_jags))
# sj  <- fit_stan$summary("Raw_Est")
# m_s <- sj$mean; se_s <- sj$sd / sqrt(sj$ess_bulk)
# 
# z <- (m_s - m_j) / sqrt(se_j^2 + se_s^2)     # standardised difference
# rel <- (m_s - m_j) / m_j
# 
# round(summary(z), 2)
# sum(abs(z) > 2)                # how many exceed 2 MCSE
# #systematic bias
# plot(m_j, m_s); abline(0, 1, col = "red")
# 
# sum(m_s > m_j)                                     # want ~17, not 34 or 0
# binom.test(sum(m_s > m_j), length(m_s))$p.value    # want non-significant
# summary((m_s - m_j) / sqrt(se_j^2 + se_s^2))       # want |z| mostly < 2
# 
# mean(fit_stan$summary("Raw_Est")$mean) / mean(apply(jm.out$jm$sims.list$Raw.Est, 2, mean))
# sum(m_s > m_j)
# 
# # and check the two correction factors directly
# fit_stan$summary("corr_factor")$mean
# 
# fit_stan$summary("sigma_Obs")$mean
# median(jm.out$jm$sims.list$sd.obs)     # was 0.243
# 
# 
# # implied mean detection in each fit
# mean(plogis(jm.out$jm$sims.list$alpha))           # JAGS, before covariates
# fit_stan$summary(c("BF_Fixed","VS_Fixed","logit_p0","sigma_Obs"))$mean
# 
# # vs JAGS BF.Fixed / VS.Fixed
# summary(stan_data$bf); summary(stan_data$vs)   # raw ranges? vs shifted?
# 
# # --- Get Summaries for Year-Specific Parameters ---
# # This will print the estimates for every year automatically
# richards_summary <- fit_stan$summary(
#   variables = c("Max", "S1", "S2", "K", "P1"),
#   "mean", "median", "quantile~0.025", "quantile~0.975", "rhat"
# )
# print(richards_summary)
# 
# library(bayesplot)
# library(ggplot2)
# 
# # Extract draws for the parameters you want to check
# # (Using a few years of Max and S1 as an example)
# draws_diagnostic <- fit_stan$draws(c("mean_prob", "Max[1]", "Max[10]", "S1[1]"))
# 
# # Generate trace plots
# mcmc_trace(draws_diagnostic) +
#   theme_minimal() +
#   labs(title = "MCMC Chain Trace Plots")
# 
# # Extract all S1 year estimates
# s1_draws <- fit_stan$draws("S1")
# 
# # Plot the intervals across all years sequentially
# mcmc_intervals(s1_draws) +
#   theme_minimal() +
#   labs(title = "Posterior Estimates for S1 Across Years",
#        x = "Parameter Value", y = "Year Index")
# 
# # Compare the fixed effect coefficients
# covariate_draws <- fit_stan$draws(c("BF_Fixed", "VS_Fixed"))
# 
# mcmc_areas(covariate_draws, prob = 0.95) +
#   theme_minimal() +
#   labs(title = "Posterior Distributions of Detection Covariates",
#        x = "Effect Size")
# 
# # 1. Extract the mean_N matrix from the posterior
# # format = "draws_matrix" converts it into a clean math grid
# mean_N_posterior <- fit_stan$draws("mean_N", format = "draws_matrix")
# 
# # 2. Pick a specific year to visualize (e.g., Year 5)
# target_year <- 5
# day_columns <- paste0("mean_N[", 1:n.days, ",", target_year, "]")
# 
# # Extract only the columns belonging to that year
# year_data <- mean_N_posterior[, day_columns]
# 
# # 3. Calculate the mean and 95% Credible Interval for each day
# trajectory <- data.frame(
#   Day = 1:n.days,
#   Mean = colMeans(year_data),
#   Lower = apply(year_data, 2, quantile, probs = 0.025),
#   Upper = apply(year_data, 2, quantile, probs = 0.975)
# )
# 
# # 4. Plot the migration curve!
# ggplot(trajectory, aes(x = Day, y = Mean)) +
#   geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "blue", alpha = 0.15) +
#   geom_line(color = "blue", size = 1) +
#   theme_minimal() +
#   labs(title = paste("Estimated Gray Whale Migration Curve: Year", target_year),
#        x = "Days into Season",
#        y = "Estimated Expected Abundance (mean_N)")
# 
# 
# # --- Posterior Predictive Simulation Loop ---
# library(bayesplot)
# library(posterior)
# 
# cat("Extracting posterior draws...\n")
# # Extract parameters as flat iteration-by-column matrices
# mean_N_mat <- as.matrix(fit_stan$draws("mean_N", format = "matrix"))
# r_mat      <- as.matrix(fit_stan$draws("r", format = "matrix"))
# p_mat      <- as.matrix(fit_stan$draws("obs_prob", format = "matrix"))
# 
# S <- nrow(p_mat)       # Total MCMC iterations saved
# V <- nrow(flat_df)     # Total unique data observations
# yrep <- matrix(0, nrow = S, ncol = V)
# 
# # Pre-calculate column lookups so R doesn't have to search text inside the loop
# mean_N_cols <- match(paste0("mean_N[", flat_df$day_idx, ",", flat_df$year_idx, "]"), colnames(mean_N_mat))
# r_cols      <- match(paste0("r[", flat_df$day_idx, ",", flat_df$year_idx, "]"), colnames(r_mat))
# 
# cat("Simulating replicated datasets...\n")
# for (s in 1:S) {
#   # Grab the parameters for this specific MCMC step across all observations
#   mu_val  <- mean_N_mat[s, mean_N_cols]
#   phi_val <- r_mat[s, r_cols]
#   p_val   <- p_mat[s, ]
#   
#   # Structural safety floor
#   mu_val[mu_val < 1e-6] <- 1e-6
#   
#   # Step 1: Simulate the true unobserved abundance (N) for this iteration
#   # (NOTE: For your Poisson models, swap this line to: N_sim <- rpois(V, lambda = mu_val))
#   N_sim <- rnbinom(V, size = phi_val, mu = mu_val)
#   
#   # Step 2: Simulate the observation process to get final whale counts
#   yrep[s, ] <- rbinom(V, size = N_sim, prob = p_val)
# }
# 
# # Define your actual raw observed data vector
# y <- flat_df$n
# 
# # We look at the first 50 simulated datasets to keep the plot clean
# ppc_dens_overlay(y, yrep[1:50, ]) +
#   theme_minimal() +
#   labs(title = "Posterior Predictive Overlays",
#        x = "Whale Count Value", y = "Density")
# 
# # 1. Does the model accurately predict the maximum number of whales seen in a single survey?
# plot_max <- ppc_stat(y, yrep, stat = "max") + 
#   labs(title = "Checking Maximum Values")
# 
# # 2. Does the model accurately handle zero-inflation (proportion of zero counts)?
# prop_zero <- function(x) mean(x == 0)
# plot_zero <- ppc_stat(y, yrep, stat = "prop_zero") + 
#   labs(title = "Checking Proportion of Zeros")
# 
# # View them side by side
# library(gridExtra)
# grid.arrange(plot_max, plot_zero, ncol = 2)
# 
# # Because looking at thousands of rows at once is overwhelming, 
# # let's look at a subset of 100 observations to see the fit clearly.
# subset_idx <- 1:100
# 
# ppc_intervals(y[subset_idx], yrep[, subset_idx]) +
#   theme_minimal() +
#   labs(title = "Predictive Intervals vs. Observed Counts (Subset)",
#        x = "Observation Index", y = "Whale Count")