#Stan_ModelComparison
#
# Compares Stan models using Richards function using LOOIC, Pareto-K, 
# rank-normalized Rhat (Rhat), and effective sample sizes (ESS).
# 

rm(list = ls())
source("Granite_Canyon_Counts_fcns.R")
library(tidyverse)
library(posterior)
library(ggplot2)
library(bayesplot)
library(loo)
library(rstanarm)

options(mc.cores = parallel::detectCores())

# The following is from Gemini
# To calculate LOOIC and view your Pareto-$k$ statistics, the loo package requires a pointwise log-likelihood matrix of size $\text{MCMC Iterations} \times \text{Observations}$.
# 
# Because we integrated the latent abundance $N_{t,y}$ out of the model, multiple unique observation sessions occurring on the same day and year share that identical hidden population pool. Therefore, to perform a mathematically valid cross-validation check, we must calculate the log-likelihood at the Day level (Cluster-LOO) rather than the raw row level.
# 
# This function reconstructs that marginalized log-likelihood matrix directly from your saved parameters in R, saving you from having to modify or re-compile your Stan code.
compute_marginal_log_lik <- function(fit_obj, flat_df, start_idx, end_idx, n_days, n_year) {
  library(matrixStats) # Highly efficient log-sum-exp operations
  
  cat("Extracting parameter chains from memory...\n")
  mean_N_mat <- as.matrix(fit_obj$draws("mean_N", format = "matrix"))
  p_mat      <- as.matrix(fit_obj$draws("obs_prob", format = "matrix"))
  
  # Check if model is Negative Binomial or Poisson
  has_dispersion <- "r[1,1]" %in% colnames(as.matrix(fit_obj$draws(variables = "r[1,1]", 
                                                                   format = "matrix")))
  if(has_dispersion) {
    r_mat <- as.matrix(fit_obj$draws("r", format = "matrix"))
  }
  
  S <- nrow(mean_N_mat) # Number of MCMC draws
  
  # Identify days that actually contain observation data
  active_days <- which(start_idx <= end_idx, arr.ind = TRUE)
  num_active_days <- nrow(active_days)
  
  log_lik_mat <- matrix(NA, nrow = S, ncol = num_active_days)
  
  cat("Computing pointwise day-level marginalized log-likelihood...\n")
  for (j in 1:num_active_days) {
    t <- active_days[j, "row"]
    y <- active_days[j, "col"]
    start <- start_idx[t, y]
    end <- end_idx[t, y]
    
    y_counts <- flat_df$n[start:end]
    max_n <- max(y_counts)
    K_val <- max_n + 50 # Matches our optimized integration buffer
    
    mu_t_y <- mean_N_mat[, paste0("mean_N[", t, ",", y, "]")]
    p_cols <- paste0("obs_prob[", start:end, "]")
    p_sub  <- p_mat[, p_cols, drop = FALSE]
    
    # Generate matching matrix of observed counts for vectorization
    y_matrix <- matrix(y_counts, nrow = S, ncol = length(y_counts), byrow = TRUE)
    
    if (t == 1 || t == n_days) {
      # Boundary condition: N is fixed at 0
      if(has_dispersion) {
        log_prob_N <- dnbinom(0, size = r_mat[, paste0("r[", t, ",", y, "]")], 
                              mu = pmax(mu_t_y, 1e-6), log = TRUE)
      } else {
        log_prob_N <- dpois(0, lambda = pmax(mu_t_y, 1e-6), log = TRUE)
      }
      log_prob_obs <- rowSums(dbinom(y_matrix, size = 0, prob = p_sub, log = TRUE))
      log_lik_mat[, j] <- log_prob_N + log_prob_obs
    } else {
      # Marginalize across possible integer states of N
      N_range <- max_n:K_val
      lp_matrix <- matrix(NA, nrow = S, ncol = length(N_range))
      
      for (k in seq_along(N_range)) {
        N_val <- N_range[k]
        if(has_dispersion) {
          log_prob_N <- dnbinom(N_val, size = r_mat[, paste0("r[", t, ",", y, "]")], 
                                mu = pmax(mu_val = mu_t_y, 1e-6), log = TRUE)
        } else {
          log_prob_N <- dpois(N_val, lambda = pmax(mu_t_y, 1e-6), log = TRUE)
        }
        log_prob_obs <- rowSums(dbinom(y_matrix, size = N_val, prob = p_sub, log = TRUE))
        lp_matrix[, k] <- log_prob_N + log_prob_obs
      }
      # Vectorized stable log_sum_exp across columns
      log_lik_mat[, j] <- rowLogSumExps(lp_matrix)
    }
  }
  return(log_lik_mat)
}


min.dur <- 60
YEAR <- 2026  # the last season name

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

YEAR <- max(years)

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
        day_idx = jags.data$day[d, s, y]
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
n.year <- jags.data$n.year
n.days <- jags.data$n.days

for (y in 1:n.year) {
  for (t in 1:n.days) {
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


run.date <- "2026-06-11"
# Model name IDs
model.names <- c( "M1a1", "M2a1", "M3a1", "M4a1", 
                  "M1a2", "M2a2", "M3a2", "M4a2")

# model IDs in the manuscript is in the same order as above but the numbers are
# different:
model.ID <- c(1:length(model.names))

#max.Rhat.big <- list()
prop.big.Rhat <- n.params <- n.big.Rhat <- n.bad.Pareto <- prop.bad.Pareto <- LOOIC <- vector(mode = "numeric", length = length(model.names))

#min.ESS <- vector(mode = "numeric", length = length(model.names))
ESS.bulk <- ESS.tail <- new.Rhat <- LOO <- log_lik <- list()
k <- 5
for (k in 1:length(model.names)){
  .out <- readRDS(paste0("RData/Stan_Richards_HSSM_", 
                         model.names[k], "_1968to", YEAR, "_min", 
                         min.dur, "_", run.date, ".rds"))
  
  # Extract summaries for the core structural parameters
  param_summary <- .out$summary(variables = c("Max", "S1", "S2", 
                                             "K", "mean_prob", "BF_Fixed", 
                                             "VS_Fixed"))
  
  # 1. Check Convergence (Look for a maximum Rhat near 1.00)
  new.Rhat[[k]] <- param_summary$rhat
  #max_rhat <- max(param_summary$rhat, na.rm = TRUE)
  #cat(sprintf("Maximum Rhat: %.3f\n", max_rhat))
  
  # 2. Check Sampling Efficiency (Look for Bulk and Tail ESS > 400)
  ESS.bulk[[k]] <- param_summary$ess_bulk
  ESS.tail[[k]] <- param_summary$ess_tail
  #min_bulk_ess <- min(param_summary$ess_bulk, na.rm = TRUE)
  #min_tail_ess <- min(param_summary$ess_tail, na.rm = TRUE)
  
  # Compute log-likelihood and LOOIC:
  # fit_obj, flat_df, start_idx, end_idx, n_days, n_year
  log_lik[[k]] <- compute_marginal_log_lik(fit_obj = .out, 
                                           flat_df = flat_df, 
                                           start_idx = start_idx, 
                                           end_idx = end_idx, 
                                           n_days = n.days, 
                                           n_year = n.year)
  
  # Generate the formal LOO objects
  LOO[[k]] <- loo(log_lik[[k]], cores = 4)

}

model_comparison <- loo_compare(LOO)
print(model_comparison, simplify = FALSE)

# Posterior predictive check:
library(bayesplot)

# 1. Isolate parameters from your chosen best-performing model object
best_model <- readRDS(paste0("RData/Stan_Richards_HSSM_", 
                       model.names[5] , "_1968to", YEAR, "_min", 
                       min.dur, "_", run.date, ".rds"))

mean_N_mat <- as.matrix(best_model$draws("mean_N", format = "matrix"))
p_mat      <- as.matrix(best_model$draws("obs_prob", format = "matrix"))
has_r      <- "r[1,1]" %in% 
  colnames(as.matrix(best_model$draws(variables = "r[1,1]", format = "matrix")))

if(has_r) { r_mat <- as.matrix(best_model$draws("r", format = "matrix")) }

S <- nrow(p_mat)       
V <- nrow(flat_df)     
yrep <- matrix(0, nrow = S, ncol = V)

mean_N_cols <- match(paste0("mean_N[", flat_df$day_idx, ",", flat_df$year_idx, "]"), 
                     colnames(mean_N_mat))
if(has_r) { r_cols <- match(paste0("r[", flat_df$day_idx, ",", flat_df$year_idx, "]"), 
                            colnames(r_mat)) }

# 2. Simulate replicated survey datasets
for (s in 1:S) {
  mu_val  <- mean_N_mat[s, mean_N_cols]
  p_val   <- p_mat[s, ]
  mu_val[mu_val < 1e-6] <- 1e-6
  
  # Conditional check on state distribution type
  if(has_r) {
    N_sim <- rnbinom(V, size = r_mat[s, r_cols], mu = mu_val)
  } else {
    N_sim <- rpois(V, lambda = mu_val)
  }
  
  yrep[s, ] <- rbinom(V, size = N_sim, prob = p_val)
}

# 3. Generate diagnostic charts for the reviewer
y_observed <- flat_df$n

# Chart A: Density Overlays
ppc_dens_overlay(y_observed, yrep[1:50, ]) + 
  theme_minimal() + labs(title = "Best Model: Posterior Predictive Overlay")

# Chart B: Target Statistics (Recreating Zero Counts and Max Counts)
prop_zero <- function(x) mean(x == 0)
ppc_stat(y_observed, yrep, stat = "prop_zero") + 
  labs(title = "Frequency of Zero-Count Surveys")

ppc_stat(y_observed, yrep, stat = "max") + 
  labs(title = "Maximum Whale Counts Observed")










out.list <- list(LOOIC = LOOIC.n,
                 Rhat = new.Rhat,
                 ESS.bulk = ESS.bulk,
                 ESS.tail = ESS.tail)
saveRDS(out.list, 
        file = paste0("RData/Stan_Richards_Convergence_", YEAR, "_", run.date, ".rds"))

out.table <- data.frame(model = model.names,
                        LOOIC = LOOIC,
                        n.params = n.params,
                        p.big.Rhat = prop.big.Rhat,
                        #LOOIC = LOOIC,
                        n.big.Rhat = n.big.Rhat,
                        n.bad.Pareto = n.bad.Pareto,
                        p.bad.Pareto = prop.bad.Pareto,
                        min.ESS.bulk = min.ESS.bulk,
                        min.ESS.tail = min.ESS.tail) %>%  
  
  mutate(dLOOIC = LOOIC - min(LOOIC, na.rm = T)) %>%
  arrange(dLOOIC)  %>%
  select(model, dLOOIC, n.params, n.big.Rhat, 
         p.big.Rhat, n.bad.Pareto, p.bad.Pareto, 
         min.ESS.bulk, min.ESS.tail, LOOIC)

saveRDS(out.table,
        file = paste0("RData/Richards_ModelComparison_", YEAR, "_", run.date, ".rds"))
#d.t.2 <- Sys.time() - t.2 # 1.68 min

