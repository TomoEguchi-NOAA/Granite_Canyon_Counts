# gap analysis
# Fit the best model to artifical datasets with gaps in different parts of the raw data.
# 
# Tests robustness of the model.
# 
library(cmdstanr)
library(tidyverse)
library(ggplot2)

# Function to evaluate the robustness of a model to missing data
calculate_gap_metrics <- function(full_draws, gap_draws) {
  
  # 1. Percent Bias of the Median (The Shift Test)
  med_full <- median(full_draws)
  med_gap  <- median(gap_draws)
  pct_bias <- abs(med_gap - med_full) / med_full * 100
  
  # 2. CI Expansion Ratio (The Uncertainty Explosion Test)
  ci_full <- quantile(full_draws, probs = c(0.025, 0.975))
  ci_gap  <- quantile(gap_draws, probs = c(0.025, 0.975))
  
  width_full <- ci_full[2] - ci_full[1]
  width_gap  <- ci_gap[2] - ci_gap[1]
  
  expansion_ratio <- width_gap / width_full
  
  # 3. Posterior Probability Overlap (Strict Distribution Test)
  # Calculates the percentage of the gapped draws that fall inside 
  # the 50% credible interval (the interquartile range) of the full model.
  ci_50_full <- quantile(full_draws, probs = c(0.25, 0.75))
  overlap_pct <- mean(gap_draws >= ci_50_full[1] & gap_draws <= ci_50_full[2]) * 100
  
  # Compile the results into a clean dataframe
  results <- data.frame(
    Metric = c("Percent Bias (%)", 
               "CI Expansion Ratio", 
               "Draws inside Full 50% CI (%)"),
    Value = c(round(pct_bias, 2), 
              round(expansion_ratio, 2), 
              round(overlap_pct, 2)),
    Passing_Grade = c("< 5 - 10%", 
                      "< 1.2 - 1.5", 
                      "High is better")
  )
  
  return(results)
}


best.model.name <- 
stan.info <- readRDS("RData//Richards_HSSM_M1a2_0gamma_mod5_ZI1_stan_info.rds")
stan.data <- stan.info$stan.data

# Try creating artificial gaps and see how it performs:

# Three gap scenarios are tested with 2 seasons within each scenario = 6 runs
# 1. The peak gap - the middle of the season days 48 - 52 missing from years 1987
# (a high year, year 16) and 2022 (a low year, year 31)
# 2. The shoulder gap - days 33-37 missing with the same two years
# 3. The tail gap - days 63-67 missing with the same two years
gap.data.1.1 <- simulate_targeted_gap(stan_data = stan.data, 
                                      target_year = 16,
                                      target_days = 48:52)

gap.data.2.1 <- simulate_targeted_gap(stan_data = stan.data, 
                                      target_year = 16,
                                      target_days = 33:37)

gap.data.3.1 <- simulate_targeted_gap(stan_data = stan.data, 
                                      target_year = 16,
                                      target_days = 63:67)

gap.data.1.2 <- simulate_targeted_gap(stan_data = stan.data, 
                                      target_year = 31,
                                      target_days = 48:52)

gap.data.2.2 <- simulate_targeted_gap(stan_data = stan.data, 
                                      target_year = 31,
                                      target_days = 33:37)

gap.data.3.2 <- simulate_targeted_gap(stan_data = stan.data, 
                                      target_year = 31,
                                      target_days = 63:67)

# Running the mod4 version, without Zero Inflation:
mod.file <- file.path(paste0("models//model_Richards_HSSM_mod5_ZI.stan"))


mod <- cmdstan_model(mod.file, 
                     cpp_options = list(stan_threads = TRUE, 
                                        O = 3))

gap.data.list <- list(gap.data.1.1, gap.data.1.2, 
                      gap.data.2.1, gap.data.2.2,
                      gap.data.3.1, gap.data.3.2)

# Parameters to monitor for Poisson models:
params.a1.stan <- c("S1", "S2", "P", 
                    "sigma_proc_P", "Corrected_Est", "Max", "log_N_latent",
                    "peak_day_decade", "beta1_P")

# Negative Binomial model has the over dispersion parameter
params.a2.stan <- c(params.a1.stan, "phi")

if (length(grep("a1", best.model.name)) > 0) params <- params.a1.stan
if (length(grep("a2", best.model.name)) > 0) params <- params.a2.stan

stan.global.summary.gap <- diag.summary.gap <- list()
for (k in 1:length(gap.data.list)){
  out.file <- paste0("Richards_HSSM_", best.model.name, "_gap", k, "_stan")
  n_year <- gap.data.list[[k]]$data$n_year
  n_observer <-gap.data.list[[k]]$data$n_observer
  sd <- gap.data.list[[k]]$data
  #mod <- cmdstan_model(file)
  if (!file.exists(paste0("RData//", out.file, ".rds"))){
    fit_stan <- mod$sample(
      data            = gap.data.list[[k]]$data,
      init            = stan_init_fn,
      chains          = 4,
      parallel_chains = 4,
      threads_per_chain = 2,
      iter_warmup     = 1500,
      iter_sampling   = 2000,
      adapt_delta     = 0.90  
    )
    
    fit_stan$save_object(file = paste0("RData//", out.file, ".rds"))
    saveRDS(list(model.file = mod.file,
                 stan.data = gap.data.list[[k]],
                 init_fn = stan_init_fn(),
                 System = Sys.getenv(),
                 Run.date = Sys.Date()),
            file = paste0("RData//", out.file, "_info.rds"))
    
  } else {
    fit_stan <- readRDS(paste0("RData//", out.file, ".rds"))
  }  
  
  # If zero-inflated, may have one or two extra parameters:
  #if (length(grep("ZI1", models[k])) > 0) params <- c(params, "zi_a")
  #if (length(grep("ZI2", models[k])) > 0) params <- c(params, "zi_a", "zi_b")
  
  stan.global.summary[[k]] <- fit_stan$summary(
    variables = params,
    default_summary_measures(), 
    default_convergence_measures(),
    extra_quantiles = ~quantile2(., probs = c(0.025, 0.975)))
  
  diag.summary[[k]] <- fit_stan$diagnostic_summary()
  
}


Corrected.Estimates.gap <- stan.global.summary %>%
  filter(str_detect(variable, "Corrected")) %>%
  mutate(Start_Years = start.years,
         Data = "Modified")

# Compare the summary statistics with the original fit with all data
fit_stan_orig <- readRDS(paste0("RData//Richards_HSSM_", best.model.name, "_stan.rds"))
stan.global.summary.orig <- fit_stan_orig$summary(
  variables = params,
  default_summary_measures(), 
  default_convergence_measures(),
  extra_quantiles = ~quantile2(., probs = c(0.025, 0.975)))

Corrected.Estimates.Orig <- stan.global.summary.orig %>%
  filter(str_detect(variable, "Corrected")) %>%
  mutate(Start_Years = start.years,
         Data = "Original")

Corrected.Estimates <- rbind(Corrected.Estimates.Orig, Corrected.Estimates.Sim)

ggplot(Corrected.Estimates) +
  geom_pointrange(aes(x = Start_Years, y = median, ymin = q2.5, ymax = q97.5, color = Data))


### Check the following chunk - Gemini created this ###
library(rstan)

# Define the year you injected the gap into (e.g., Year 1)
target_year <- 1

# Extract the Corrected_Est matrix from both models
# rows = MCMC iterations, columns = years
draws_full <- extract(fit_full_data)$Corrected_Est
draws_peak <- extract(fit_peak_gap)$Corrected_Est

# Isolate the vector of draws for the specific target year
target_draws_full <- draws_full[, target_year]
target_draws_peak <- draws_peak[, target_year]

# Run the evaluation
metrics_peak_gap <- calculate_gap_metrics(full_draws = target_draws_full, 
                                          gap_draws  = target_draws_peak)

print(metrics_peak_gap)

