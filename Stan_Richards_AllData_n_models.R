# Running the best model in STAN
# 

rm(list = ls())
library(tidyverse)
library(cmdstanr)
library(loo)

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
  
  # Hyper-parameters
  S1_alpha     = 10.0,  # Match your uniform hyperprior midpoint bounds
  S1_beta      = 1.0,
  S2_alpha     = 10.0,
  S2_beta      = 1.0,
  log_K_mu     = 1.5,
  log_K_sigma  = 0.5
)

# --- 4. Compile and Execute ---
models <- c("M1a2", "M2a2", "M3a2", "M4a2",
            "M1a1", "M2a1", "M3a1", "M4a1")

fit_stan <- list()
for (k in 1:length(models)){
  model.file <- paste0("models/model_Richards_HSSM_", 
                       models[k], ".stan")
  
  out.file <- paste0("Stan_Richards_HSSM_", 
                     models[k], "_1968to", YEAR, "_min", 
                     min.dur, "_", Run.date)
  
  if (!file.exists(out.file)){
    # Compile with aggressive C++ optimization flags
    mod <- cmdstan_model(model.file, 
                         cpp_options = list(stan_threads = TRUE, 
                                            O = 3))
    
    tic <- Sys.time()
    fit_ <- mod$sample(
      data            = stan_data,
      chains          = 4,
      parallel_chains = 4,
      threads_per_chain = 2,
      iter_warmup     = 1000,
      iter_sampling   = 1000,
      adapt_delta     = 0.95  # High adaptation limit targets complex curves cleanly
    )
    toc <- Sys.time() - tic
    # This forces R to read the CSVs and compress them into a permanent file
    
    fit_$save_object(file = paste0("RData/", out.file, ".rds"))
    
    info.stan <- list(out.filename = out.file,
                      Run.Time = toc,
                      Run.Date = Sys.Date(),
                      System = Sys.getenv())
     
     saveRDS(info.stan,
             file = paste0("RData/", out.file, ".info"))
    
  }
}

# Output needs to be modified to fit the following code chunk using loo_compare() 

# Compute the LOO object for a model
loo_model1 <- fit_stan$loo(cores = 4)
print(loo_model1)

# Once you run this for multiple models, you can compare them directly:

comparison <- loo_compare(loo_model1, loo_model2, loo_model3)
print(comparison)

