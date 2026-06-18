# Running the best model in STAN
# 

rm(list = ls())
library(tidyverse)

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

jags.input.list <- AllData2JagsInput_NoBUGS(min.dur, years = years, data.dir, max.day)                        
#jags.input.list$jags.data["N"] <- NULL
# Modify jags data to rearrange days and provide zeros for t = 1 and t = max.day
jags.data <- jags.input.list$jags.data

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
        obs = jags.data$obs[d, s, y],
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

# --- Add these lines right before packaging 'stan_data' ---
storage.mode(start_idx) <- "integer"
storage.mode(end_idx)   <- "integer"

# --- 3. Package Everything for Stan ---
stan_data <- list(
  n_year       = jags.data$n.year,
  n_days       = jags.data$n.days,
  n_obs        = max(flat_df$obs),
  N_flat       = nrow(flat_df),
  n            = flat_df$n,
  bf           = flat_df$bf,
  vs           = flat_df$vs,
  obs          = flat_df$obs,
  watch_length = flat_df$watch_length,
  start_idx    = start_idx,
  end_idx      = end_idx,
  day_idx      = flat_df$day_idx,
  year_idx     = flat_df$year_idx,
  station_idx  = flat_df$station_idx,
  year_values  = c(1:jags.data$n.year)
)

# --- 4. Compile and Execute ---
#file <- file.path("models/whale_model.stan")
file <- file.path("models/model_Richards_HSSM_M1a2_reparam.stan")
# Compile with aggressive C++ optimization flags
mod <- cmdstan_model(file, 
                     cpp_options = list(stan_threads = TRUE, O = 3))

#mod <- cmdstan_model(file)

fit_stan <- mod$sample(
  data            = stan_data,
  chains          = 4,
  parallel_chains = 4,
  threads_per_chain = 2,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  adapt_delta     = 0.95  # High adaptation limit targets complex curves cleanly
)

# --- 5. Inspect Results ---
print(fit_stan$summary(c("S1", "S2", "P", "BF_Fixed", "VS_Fixed", "Max")))


# --- Get Summaries for Specific Global Parameters ---
global_summary <- fit_stan$summary(
  variables = c("mean_prob", "BF_Fixed", "VS_Fixed", "sigma_Obs"),
  "mean", "sd", "rhat", "ess_bulk"
)
print(global_summary)

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
