# Running the best model in STAN
# 

rm(list = ls())
library(tidyverse)
library(posterior)

# Need to install cmdnstanr from the Stan package repo:
# # Tell R to look at the Stan repository instead of CRAN
#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)

# This downloads and builds the C++ toolchain Stan needs to run fast
#check_cmdstan_toolchain(fix = TRUE) # Verifies you have a C++ compiler
#install_cmdstan(cores = 4)          # Downloads and installs CmdStan

source("Granite_Canyon_Counts_fcns.R")

# Create data input for JAGS - this will have to be changed later:
# Minimum length of observation periods in minutes
min.dur <- 60 #10 #85 #
Run.date <- Sys.Date()

# These are the ending year of each season - for example, 2022 in the following vector indicates
# for the 2021/2022 season. These data were extracted using Extract_Data_All_v2.Rmd
# Data prior to the 2009/2010 season are in Laake's ERAnalayis package. 
years <- c(2008, 2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024, 2025, 2026)
data.dir <- "RData/V2.1_May2026"
max.day <- 100

#jags.input.list <- AllData2JagsInput_NoBUGS(min.dur, years = years, data.dir, max.day)                       
jm.out <- readRDS("RData/JAGS_Richards_HSSM_M1a2_1968to2026_min60_2026-07-27_NoBUGS.rds")

#jags.input.list$jags.data["N"] <- NULL
# Modify jags data to rearrange days and provide zeros for t = 1 and t = max.day
jags.data <- jm.out$jags.input$jags.data

# --- 1. Flatten Your Existing JAGS Arrays ---
# (This simulates rebuilding your messy arrays into a long-form data frame)
flat_data_list <- list()
counter <- 1

# add 2 (day 1 and 100) to # periods:
periods <- jags.data$periods + 2

for (y in 1:jags.data$n.year) {
  for (s in 1:jags.data$n.station[y]) {
    for (d in 2:(periods[y, s] - 1)) {
      
      flat_data_list[[counter]] <- data.frame(
        n = jags.data$n[d, s, y],
        bf = jags.data$bf[(d - 1), s, y],
        vs = jags.data$vs[(d - 1), s, y],
        obs = jags.data$obs.fixed[d, s, y],
        watch_length = jags.data$watch.length[(d - 1), s, y],
        year_idx = y,
        day_idx = jags.data$day[d, s, y],
        station_idx = s
      )
      counter <- counter + 1
    }
  }
}

flat_df <- do.call(rbind, flat_data_list)
# Ensure data is sorted sequentially by year and day for indexing
flat_df <- flat_df[order(flat_df$year_idx, flat_df$day_idx), ]

# --- 2. Build the Start/End Pointer Index Matrices ---
start_idx <- matrix(0, nrow = jags.data$n.days, ncol = jags.data$n.year)
end_idx <- matrix(0, nrow = jags.data$n.days, ncol = jags.data$n.year)

for (y in 1:jags.data$n.year) {
  for (t in 1:jags.data$n.days) {
    matching_rows <- which(flat_df$year_idx == y & flat_df$day_idx == t)
    if (length(matching_rows) > 0) {
      start_idx[t, y] <- min(matching_rows)
      end_idx[t, y]  <- max(matching_rows)
    } else {
      start_idx[t, y] <- 1
      end_idx[t, y]  <- 0 # Signifies no observations on this day
    }
  }
}

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

# --- Add these lines right before packaging 'stan_data' ---
storage.mode(start_idx) <- "integer"
storage.mode(end_idx)   <- "integer"

# --- 3. Package Everything for Stan ---
stan_data <- list(
  n_year = jags.data$n.year, 
  n_days = jags.data$n.days, 
  n_observer = max(flat_df$obs), 
  N_flat = nrow(flat_df),
  n = flat_df$n, 
  bf = flat_df$bf, 
  vs = flat_df$vs - 1,      # RAW, and vs shifted
  observer_idx = flat_df$obs, 
  watch_length = flat_df$watch_length,
  year_values = jm.out$jags.input$jags.data$year.index, 
  day_idx = flat_df$day_idx, 
  year_idx = flat_df$year_idx,
  
  sd_beta1_P = 2, beta0_Max_mu = 7.6, sd_beta0_Max = 2, sd_beta1_Max = 2,
  sd_BF = 2, sd_VS = 2, sd_sigma_proc_P = 5, sd_sigma_proc_Max = 5,
  alpha_S_mu = 10, alpha_S_sd = sqrt(10), beta_S_shape = 1, beta_S_rate = 1,
  sigma_Obs_max = 1.5, phi_max = 50,
  gamma_fix = 1.0, mean_prob_mu = 0.8, mean_prob_sd = 0.35,
  year_specific_phi = 0,
  anchor_mu = 1.39, anchor_sd = 0.01, #0.1622,    # 0.01 = parity; raise to 0.1622 to propagate
  boundary_N = 0.0001,
  centred_P = 1, centred_Max = 1,
  use_process_error = 0, use_shape_dev = 0, 
  independent_corr = 0,
  n_period = 20, period_idx = rep(1:20, each = 5), sd_sigma_shape = 0.5
)

# N(1.3863, 0.162) on the logit scale implies p with a 95% CI of [0.744, 0.846]
# and median 0.8. This should be used for the distribution of P. 

# note that mean_prob_mu should be tested for 0.7 and 0.9, as sensitivity tests
n_year <- jags.data$n.year
n_observer <- jags.data$n.obs.fixed
init_fn <- function() list(
  beta0_Max = rnorm(1, 7.6, 0.2),  
  beta1_Max = rnorm(1, 0, 0.05),
  #log_Max_raw = rnorm(n_year, 0, 0.1),  
  sigma_proc_Max = runif(1, 0.2, 0.6),
  beta0_P = runif(1, 42, 48),      
  beta1_P = rnorm(1, 0.21, 0.03),
  P_raw       = rnorm(n_year, 45, 2),      # was rnorm(n_year, 0, 0.1)
  log_Max_raw = rnorm(n_year, 7.6, 0.2),    # was rnorm(n_year, 0, 0.1)
  #P_raw = rnorm(n_year, 0, 0.1),   
  #sigma_proc_P = runif(1, 3.5, 5.5),
  S1 = runif(n_year, 2, 4),        
  S2 = runif(n_year, 2, 4),
  mu_S1 = runif(1, 2.5, 3.5),  
  shape_S1 = runif(1, 8, 12),
  mu_S2 = runif(1, 2.5, 3.5),  
  shape_S2 = runif(1, 8, 12),
  #alpha_S1 = runif(1, 8, 12),      beta_S1 = runif(1, 2.5, 4.5),
  #alpha_S2 = runif(1, 8, 12),      beta_S2 = runif(1, 2.5, 4.5),
  alpha = rnorm(n_observer, 1.39, 0.1),  sigma_Obs = runif(1, 0.15, 0.35),
  logit_p0 = rnorm(1, 1.3863, 0.01),
  BF_Fixed = rnorm(1, 0, 0.1),     VS_Fixed = rnorm(1, 0, 0.1),
  phi = runif(1, 4, 6),
  sigma_shape = runif(1, 0.05, 0.2)
)

# --- 4. Compile and Execute ---
#file <- file.path("models/whale_model.stan")
file <- file.path("models/model_Richards_HSSM_M1a2_fixed_1.stan")
out.file <- "Richards_HSSM_M1a2_stan"

# Compile with aggressive C++ optimization flags
mod <- cmdstan_model(file, 
                     cpp_options = list(stan_threads = TRUE, O = 3))

#mod <- cmdstan_model(file)
if (!file.exists(paste0("RData/", out.file, ".rds"))){
  fit_stan <- mod$sample(
    data            = stan_data,
    init            = init_fn,
    chains          = 4,
    parallel_chains = 4,
    threads_per_chain = 2,
    iter_warmup     = 1000,
    iter_sampling   = 1000,
    adapt_delta     = 0.90  
  )
  
  fit_stan$save_object(file = paste0("RData/", out.file, ".rds"))
  
} else {
  fit_stan <- readRDS(paste0("RData/", out.file, ".rds"))
}

# --- 5. Inspect Results ---
params.1.stan <- c("S1", "beta_S1", "S2", "P", "phi",
                   "sigma_proc_P", "Corrected_Est", "Max")

# --- Get Summaries for Specific Global Parameters ---
stan_global_summary <- fit_stan$summary(
  variables = params.1.stan,
  default_summary_measures(), 
  default_convergence_measures(),
  extra_quantiles = ~quantile2(., probs = c(0.025, 0.975)))

print(stan_global_summary)

# From Jags using regular expression syntax:
params.1.jags <- c("^S1\\[", "^S1.beta", "^S2\\[", "^P\\[", "^r",
                   "^sd.proc.P", "^Corrected.Est", "^Max\\[")

jags.params.idx <- lapply(params.1.jags,
                     FUN = function(x) grep(x, jm.out$posterior.summary$variable)) %>% unlist()

summarise_draws(jm.out$jm$samples, 
                default_summary_measures(),
                default_convergence_measures(),
                extra_quantiles = ~quantile2(., probs = c(0.025, 0.975))) %>%
  as.data.frame() -> jags.convergence.summary

jags.global.summary <- jags.convergence.summary[jags.params.idx,]

stan.S1 <- stan_global_summary[grep("^S1\\[", stan_global_summary$variable),] %>%
  mutate(MCMC = "Stan",
         year.idx = seq(1, n_year))
jags.S1 <- jags.global.summary[grep("^S1\\[", jags.global.summary$variable),] %>%
  mutate(MCMC = "Jags",
         year.idx = seq(1, n_year))

S1.stan.jags <- rbind(stan.S1, jags.S1)

ggplot(S1.stan.jags) +
  geom_point(aes(x = year.idx, y = mean, color = MCMC)) +
  geom_errorbar(aes(x = year.idx, ymin = q2.5, ymax = q97.5)) +
  theme(legend.position = "top")

stan.Nhats <- stan_global_summary[grep("^Corrected_Est", stan_global_summary$variable),] %>%
  mutate(MCMC = "Stan",
         year.idx = seq(1, n_year))
jags.Nhats <- jags.global.summary[grep("^Corrected.Est", jags.global.summary$variable),] %>%
  mutate(MCMC = "Jags",
         year.idx = seq(1, n_year))

Nhats.stan.jags <- rbind(stan.Nhats, jags.Nhats)

ggplot(Nhats.stan.jags) +
  geom_point(aes(x = year.idx, y = mean, color = MCMC)) +
  geom_errorbar(aes(x = year.idx, ymin = q2.5, ymax = q97.5)) +
  theme(legend.position = "top")

# --- Get Summaries for Year-Specific Parameters ---
# This will print the estimates for every year automatically
richards_summary <- fit_stan$summary(
  variables = c("Max", "S1", "S2", "K", "P1"),
  "mean", "median", "quantile~0.025", "quantile~0.975", "rhat"
)
print(richards_summary)

library(bayesplot)
library(ggplot2)

# Extract draws for the parameters you want to check
# (Using a few years of Max and S1 as an example)
draws_diagnostic <- fit_stan$draws(c("mean_prob", "Max[1]", "Max[10]", "S1[1]"))

# Generate trace plots
mcmc_trace(draws_diagnostic) +
  theme_minimal() +
  labs(title = "MCMC Chain Trace Plots")

# Extract all S1 year estimates
s1_draws <- fit_stan$draws("S1")

# Plot the intervals across all years sequentially
mcmc_intervals(s1_draws) +
  theme_minimal() +
  labs(title = "Posterior Estimates for S1 Across Years",
       x = "Parameter Value", y = "Year Index")

# Compare the fixed effect coefficients
covariate_draws <- fit_stan$draws(c("BF_Fixed", "VS_Fixed"))

mcmc_areas(covariate_draws, prob = 0.95) +
  theme_minimal() +
  labs(title = "Posterior Distributions of Detection Covariates",
       x = "Effect Size")

# 1. Extract the mean_N matrix from the posterior
# format = "draws_matrix" converts it into a clean math grid
mean_N_posterior <- fit_stan$draws("mean_N", format = "draws_matrix")

# 2. Pick a specific year to visualize (e.g., Year 5)
target_year <- 5
day_columns <- paste0("mean_N[", 1:n.days, ",", target_year, "]")

# Extract only the columns belonging to that year
year_data <- mean_N_posterior[, day_columns]

# 3. Calculate the mean and 95% Credible Interval for each day
trajectory <- data.frame(
  Day = 1:n.days,
  Mean = colMeans(year_data),
  Lower = apply(year_data, 2, quantile, probs = 0.025),
  Upper = apply(year_data, 2, quantile, probs = 0.975)
)

# 4. Plot the migration curve!
ggplot(trajectory, aes(x = Day, y = Mean)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "blue", alpha = 0.15) +
  geom_line(color = "blue", size = 1) +
  theme_minimal() +
  labs(title = paste("Estimated Gray Whale Migration Curve: Year", target_year),
       x = "Days into Season",
       y = "Estimated Expected Abundance (mean_N)")


# --- Posterior Predictive Simulation Loop ---
library(bayesplot)
library(posterior)

cat("Extracting posterior draws...\n")
# Extract parameters as flat iteration-by-column matrices
mean_N_mat <- as.matrix(fit_stan$draws("mean_N", format = "matrix"))
r_mat      <- as.matrix(fit_stan$draws("r", format = "matrix"))
p_mat      <- as.matrix(fit_stan$draws("obs_prob", format = "matrix"))

S <- nrow(p_mat)       # Total MCMC iterations saved
V <- nrow(flat_df)     # Total unique data observations
yrep <- matrix(0, nrow = S, ncol = V)

# Pre-calculate column lookups so R doesn't have to search text inside the loop
mean_N_cols <- match(paste0("mean_N[", flat_df$day_idx, ",", flat_df$year_idx, "]"), colnames(mean_N_mat))
r_cols      <- match(paste0("r[", flat_df$day_idx, ",", flat_df$year_idx, "]"), colnames(r_mat))

cat("Simulating replicated datasets...\n")
for (s in 1:S) {
  # Grab the parameters for this specific MCMC step across all observations
  mu_val  <- mean_N_mat[s, mean_N_cols]
  phi_val <- r_mat[s, r_cols]
  p_val   <- p_mat[s, ]
  
  # Structural safety floor
  mu_val[mu_val < 1e-6] <- 1e-6
  
  # Step 1: Simulate the true unobserved abundance (N) for this iteration
  # (NOTE: For your Poisson models, swap this line to: N_sim <- rpois(V, lambda = mu_val))
  N_sim <- rnbinom(V, size = phi_val, mu = mu_val)
  
  # Step 2: Simulate the observation process to get final whale counts
  yrep[s, ] <- rbinom(V, size = N_sim, prob = p_val)
}

# Define your actual raw observed data vector
y <- flat_df$n

# We look at the first 50 simulated datasets to keep the plot clean
ppc_dens_overlay(y, yrep[1:50, ]) +
  theme_minimal() +
  labs(title = "Posterior Predictive Overlays",
       x = "Whale Count Value", y = "Density")

# 1. Does the model accurately predict the maximum number of whales seen in a single survey?
plot_max <- ppc_stat(y, yrep, stat = "max") + 
  labs(title = "Checking Maximum Values")

# 2. Does the model accurately handle zero-inflation (proportion of zero counts)?
prop_zero <- function(x) mean(x == 0)
plot_zero <- ppc_stat(y, yrep, stat = "prop_zero") + 
  labs(title = "Checking Proportion of Zeros")

# View them side by side
library(gridExtra)
grid.arrange(plot_max, plot_zero, ncol = 2)

# Because looking at thousands of rows at once is overwhelming, 
# let's look at a subset of 100 observations to see the fit clearly.
subset_idx <- 1:100

ppc_intervals(y[subset_idx], yrep[, subset_idx]) +
  theme_minimal() +
  labs(title = "Predictive Intervals vs. Observed Counts (Subset)",
       x = "Observation Index", y = "Whale Count")
