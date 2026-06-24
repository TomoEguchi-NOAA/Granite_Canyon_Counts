
rm(list = ls())

# Posterior predictive check 
library(bayesplot)
library(ggplot2)
library(reshape2)

# Bring in the best model: this should be found after running Jags_Richards_ModelComparison.R
# which uses tjhe results from Jags_Richards_AllData_n_models.R
best.model <- "M1a2"

min.dur <- 60
YEAR <- 2026  # the last season name
run.date <- "2026-06-18"
jags_output <- readRDS(paste0("RData/JAGS_Richards_HSSM_", 
                              best.model, "_1968to", YEAR, "_min", 
                              min.dur, "_", run.date, "_NoBUGS.rds"))

jags_fit <- jags_output$jm

# ==============================================================================
# 1. ASSUMPTIONS & PREREQUISITES
# ==============================================================================
# This script assumes:
#   1. Your observed data array 'n' is available in your R environment with dimensions [d, s, y]
#   2. You included "kappa" and "r" in the 'variable.names' monitored during your JAGS run.
#   3. 'jags_fit' is your output object (assumed here to be from jagsUI or R2jags)
# ==============================================================================

# Extract the posterior simulations
# Dimensions of sims_kappa will be: [MCMC_iterations, max_days, max_stations, n_years]
sims_kappa <- jags_fit$sims.list$kappa  
sims_r     <- jags_fit$sims.list$r      # Dimensions: [MCMC_iterations, n_years]

n_sims <- dim(sims_kappa)[1]
n_days <- dim(sims_kappa)[2]
n_stat <- dim(sims_kappa)[3]
n_year <- dim(sims_kappa)[4]

# Create a flattened vector of your real observed data where observations actually occurred
# This ignores unobserved or hardcoded placeholder matrix slots
obs_indices <- which(!is.na(n) & n >= 0) 
y_obs       <- n[obs_indices]
N_data      <- length(y_obs)

# Initialize a matrix to hold your replicated datasets: [iterations x data_points]
y_rep <- matrix(NA, nrow = n_sims, ncol = N_data)

# Initialize vectors to calculate the Bayesian p-value via discrepancy metrics
T_obs <- rep(NA, n_sims)
T_rep <- rep(NA, n_sims)

# ==============================================================================
# 2. GENERATE REPLICATED DATA & CALCULATE DISCREPANCY
# ==============================================================================
cat("Simulating posterior predictive datasets...\n")

for (i in 1:n_sims) {
  
  # Reconstruct the 3D kappa matrix for this specific MCMC iteration
  kappa_iter <- sims_kappa[i, , , ]
  
  # Reconstruct the expected count matrix matching the observation structure
  # For each year, scale the dispersion parameter 'r' appropriately
  r_iter <- matrix(NA, nrow = n_days, ncol = n_stat)
  r_matrix <- array(NA, dim = c(n_days, n_stat, n_year))
  for(y in 1:n_year){
    r_matrix[,,y] <- sims_r[i, y]
  }
  
  # Flatten them using the exact same structural index filter as the observed data
  kappa_flat <- kappa_iter[obs_indices]
  r_flat     <- r_matrix[obs_indices]
  
  # Generate a random draw from the Negative Binomial distribution for every single data point
  # Protecting against minor numerical approximations near 0 using pmax
  y_rep[i, ] <- rnbinom(n = N_data, size = r_flat, mu = pmax(kappa_flat, 1e-6))
  
  # Freeman-Tukey Discrepancy Metric: Highly robust against zero-inflated data
  # Formula: T = sum( (sqrt(y) - sqrt(E))^2 )
  T_obs[i] <- sum((sqrt(y_obs) - sqrt(kappa_flat))^2, na.rm = TRUE)
  T_rep[i] <- sum((sqrt(y_rep[i, ]) - sqrt(kappa_flat))^2, na.rm = TRUE)
}

# Calculate the Bayesian p-value (proportion of times simulated discrepancy exceeds observed)
b_p_value <- mean(T_rep > T_obs)
cat("Bayesian p-value (Freeman-Tukey):", round(b_p_value, 3), "\n")
cat("Note: Values near 0.5 indicate excellent fit. Avoid values < 0.05 or > 0.95.\n")

# ==============================================================================
# 3. VISUALIZATION
# ==============================================================================

### Plot 1: Discrepancy Scatter Plot (The Classic PPC Visual)
df_ppc <- data.frame(T_obs = T_obs, T_rep = T_rep)
p1 <- ggplot(df_ppc, aes(x = T_obs, y = T_rep)) +
  geom_point(alpha = 0.4, color = "dodgerblue4") +
  geom_ab_line(intercept = 0, slope = 1, color = "firebrick", size = 1, linetype = "dashed") +
  labs(
    title = "Posterior Predictive Check: Discrepancy Plot",
    subtitle = paste("Bayesian p-value =", round(b_p_value, 3)),
    x = "Discrepancy of Observed Data: T(y, theta)",
    y = "Discrepancy of Replicated Data: T(y_rep, theta)"
  ) +
  theme_minimal()

print(p1)

### Plot 2: Density Overlay of Counts (Using 'bayesplot')
# Subset to a random group of 50 simulated iterations to keep the plot legible
set.seed(42)
selected_sims <- sample(1:n_sims, 50)

# Focus the plot on the lower-to-middle count range to see distribution shapes clearly
p2 <- ppc_dens_overlay(y = y_obs, yrep = y_rep[selected_sims, ]) +
  scale_x_continuous(limits = c(0, quantile(y_obs, 0.95))) +
  labs(
    title = "Posterior Predictive Density Overlay",
    subtitle = "Observed count distribution vs. 50 simulated predictive datasets (95% range focused)",
    x = "Whale Counts (n)"
  )

print(p2)