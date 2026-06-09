#Jags_Richards_AllData_n_models.R
# Runs multiple models and saves results in .rds files
# 


rm(list = ls())

library(tidyverse)
library(loo)

source("Granite_Canyon_Counts_fcns.R")
source("Richards_HSSM_model_definition.R")
options(mc.cores = parallel::detectCores())

# Minimum length of observation periods in minutes
min.dur <- 60 #10 #85 #

model.defs <- data.frame(Lkhd = c(rep("NegBin", 4), rep("Poisson", 4)),
                         P = "time",
                         S1 = rep(c("time", "S1"), 4),
                         S2 = rep(c("time", "S2", "S2", "time"), 2))

Run.date <- Sys.Date()

# These are the ending year of each season - for example, 2022 in the following vector indicates
# for the 2021/2022 season. These data were extracted using Extract_Data_All_v2.Rmd
# Data prior to the 2009/2010 season are in Laake's ERAnalayis package. 
years <- c(2008, 2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024, 2025, 2026)
data.dir <- "RData/V2.1_May2026"
max.day <- 100

#5000 samples
MCMC.params <- list(n.samples = 250000,
                    n.thin = 200,
                    n.burnin = 50000,
                    n.chains = 5)

# 225 samples
# MCMC.params <- list(n.samples = 100,
#                     n.thin = 2,
#                     n.burnin = 10,
#                     n.chains = 5)

jags.params <- c("VS.Fixed", "BF.Fixed",
                 "Max", "K", "K1", "K2", "S1", "S2", "P",
                 "mean.prob", "prob", "obs.prob",
                 "mean.N", "Corrected.Est", "N", "obs.N",
                 #"OBS.RF", "sigma.Obs",
                 "Max.alpha", "Max.beta",
                 "S1.alpha", "S2.alpha",
                 "S1.beta", "S2.beta",
                 "beta0.P", "beta1.P", "sd.proc.P",
                 "beta.p",  "sd.obs",
                 "beta0.Max", "beta1.Max", "sd.proc.Max",
                 "Raw.Est", "beta.obs",
                 "alpha", "r",
                 #"P.alpha", "P.beta",
                 #"K.alpha", "K.beta",
                 #"beta.1",
                 #"N.alpha", "N.obs",
                 "log.lkhd")
model.names <- list()
for (k in 1:nrow(model.defs)){
  model.names[[k]] <- Richards_HSSM_model_definition(K = 1,
                                                     S1 = model.defs[k, "S1"],
                                                     S2 = model.defs[k, "S2"],
                                                     P = "time",
                                                     Max = "time",
                                                     lkhd = model.defs[k, "Lkhd"])
  
  model.name.no.dir <- strsplit(model.names[[k]], split = "models/")[[1]][2]
  
  jm.out <- NoBUGS_Richards_fcn(min.dur = min.dur, 
                                years = years, 
                                data.dir = data.dir, 
                                jags.params = jags.params, 
                                MCMC.params = MCMC.params,
                                Run.date = Run.date,
                                obs.n.min = 10,
                                max.day = 100,
                                N.obs = 10,
                                model.name = model.name.no.dir,
                                ext = ".jags")
  
}

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


