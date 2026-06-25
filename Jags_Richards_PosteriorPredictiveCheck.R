
rm(list = ls())

# Posterior predictive check 
library(bayesplot)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(tidyr)

source("Granite_Canyon_Counts_fcns.R")
# Bring in the best model: this should be found after running Jags_Richards_ModelComparison.R
# which uses tjhe results from Jags_Richards_AllData_n_models.R
best.model <- "M1a2"
data.dir <- "RData/V2.1_May2026"
min.dur <- 60
max.day <- 100

YEAR <- 2026  # the last season name
run.date <- "2026-06-24"  # This is the version with "kappa" monitored. Previous runs that 
# were used to compare models (Jags_Richards_ModelComparison.R) did not include "kappa" in
# the list of monitored parameters
jags_output <- readRDS(paste0("RData/JAGS_Richards_HSSM_", 
                              best.model, "_1968to", YEAR, "_min", 
                              min.dur, "_", run.date, "_NoBUGS.rds"))

jags_fit <- jags_output$jm

# Need all seasons:
years <- c(2008, 2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024, 2025, 2026)
jags.input.list <- AllData2JagsInput_NoBUGS(min.dur, years = years, 
                                            data.dir = data.dir, 
                                            max.day = max.day)                

all.start.years <- c(jags.input.list$jags.input.Laake$all.start.year,
                     jags.input.list$jags.input.new$start.years)

seasons <- paste0(all.start.years, "/", all.start.years+1)
# ==============================================================================
# 1. ASSUMPTIONS & PREREQUISITES
# ==============================================================================
# This script assumes:
#   1. Your observed data array 'n' is available in your R environment with dimensions [d, s, y]
#   2. You included "kappa" and "r" in the 'variable.names' monitored during your JAGS run.
#   3. 'jags_fit' is your output object (assumed here to be from jagsUI or R2jags)
# ==============================================================================

# ==============================================================================
# Pointwise PPC with Explicit 1-Day Dimension Offset Correction
# ==============================================================================

n <- jags_output$jags.input$jags.data$n

# ==============================================================================
# DIRECT POINTWISE POSTERIOR PREDICTIVE CHECK (PPC) & APPENDIX EXPORTER
# ==============================================================================


#' Execute Direct Pointwise PPC with Offset Correction
#'
#' @param jags_fit The output object from your JAGS run (containing sims.list)
#' @param n The original observed 3D array of counts with dimensions [d, s, y]
#' @param day The number of days since 12/1 - stored in jags.data [d, s, y]
run_direct_ppc <- function(jags_fit, n, day) {
  
  # --- 1. EXTRACT DATA & DIMENSIONS ---
  sims_kappa <- jags_fit$sims.list$kappa  # [MCMC_iter, n_days_kappa, n_stations, n_years]
  sims_r     <- jags_fit$sims.list$r      # [MCMC_iter] (Global dispersion)
  
  n_sims        <- dim(sims_kappa)[1]
  n_days_kappa  <- dim(sims_kappa)[2]
  n_stat        <- dim(n)[2]
  n_year        <- dim(n)[3]
  n_days_n      <- dim(n)[1] # Includes Day 1 (which is anchored to zero in data)
  
  cat(sprintf("MCMC Iterations: %d\n", n_sims))
  cat(sprintf("Observed data 'n' has %d days (Includes Day 1).\n", n_days_n))
  cat(sprintf("Monitored 'kappa' has %d days (Starts at Day 2).\n", n_days_kappa))
  
  # --- 2. MULTI-LAYER COORDINATE FILTERING ---
  # Build an explicit coordinate map to link Day 'd' in data to 'd - 1' in kappa.
  # This cleanly resolves the 1-day loop index mismatch.
  matched_list <- list()
  idx <- 1
  
  for (y in 1:n_year) {
    for (s in 1:n_stat) {
      for (d in 2:n_days_n) {
        
        d_kappa <- d - 1
        
        # Pull coordinates where a real count exists AND JAGS evaluated the node
        if (!is.na(n[d, s, y]) && n[d, s, y] >= 0 && !is.na(sims_kappa[1, d_kappa, s, y])) {
          matched_list[[idx]] <- data.frame(
            Day         = day[d, s, y],
            Period      = d,
            D_Kappa     = d_kappa,
            Station     = s,
            Year        = y,
            Obs_n       = n[d, s, y]
          )
          idx <- idx + 1
        }
      }
    }
  }
  
  df_matched <- do.call(rbind, matched_list)
  N_data     <- nrow(df_matched)
  
  y_obs <- df_matched$Obs_n
  
  if (N_data == 0) {
    stop("No overlapping data points found between observed 'n' and modeled 'kappa'. Check array structures.")
  }
  
  #cat(sprintf("Successfully mapped %d perfect observation-to-model coordinate pairs.\n", N_data))
  
  # --- 3. EXTRACT POSTERIOR KAPPA VALUES EFFICIENTLY ---
  #cat("Extracting posterior trajectories across the coordinate map...\n")
  kappa_matrix <- matrix(NA, nrow = n_sims, ncol = N_data)
  for (j in 1:N_data) {
    kappa_matrix[, j] <- sims_kappa[, df_matched$D_Kappa[j], 
                                    df_matched$Station[j], 
                                    df_matched$Year[j]]
  }
  
  # --- 4. GENERATE REPLICATED SAMPLES (POINTWISE) ---
  cat("Generating pointwise simulated datasets...\n")
  y_rep <- matrix(NA, nrow = n_sims, ncol = N_data)
  T_obs <- rep(NA, n_sims)
  T_rep <- rep(NA, n_sims)
  
  for (i in 1:n_sims) {
    r_val <- sims_r[i]
    
    # Simulate a negative binomial observation for each matched node
    y_rep[i, ] <- rnbinom(n = N_data, 
                          size = r_val, 
                          mu = pmax(kappa_matrix[i, ], 1e-6))
    
    # Freeman-Tukey discrepancy metric on aligned arrays
    T_obs[i] <- sum((sqrt(df_matched$Obs_n) - sqrt(kappa_matrix[i, ]))^2)
    T_rep[i] <- sum((sqrt(y_rep[i, ]) - sqrt(kappa_matrix[i, ]))^2)
  }
  
  b_p_value <- mean(T_rep > T_obs)
  
  # --- 5. POINTWISE STATISTICAL SUMMARY (DIRECT COMPARISON) ---
  cat("Calculating pointwise summary intervals...\n")
  pointwise_summary <- df_matched %>%
    mutate(
      Rep_Median = apply(y_rep, 2, median),
      Rep_L95    = apply(y_rep, 2, quantile, probs = 0.025),
      Rep_U95    = apply(y_rep, 2, quantile, probs = 0.975),
      Rep_L50    = apply(y_rep, 2, quantile, probs = 0.25),
      Rep_U50    = apply(y_rep, 2, quantile, probs = 0.75),
      In_95_CI   = (Obs_n >= Rep_L95) & (Obs_n <= Rep_U95),
      In_50_CI   = (Obs_n >= Rep_L50) & (Obs_n <= Rep_U50)
    )
  
  coverage_95 <- mean(pointwise_summary$In_95_CI) * 100
  coverage_50 <- mean(pointwise_summary$In_50_CI) * 100
  
  cat("\n================ DIRECT COMPARISON REPORT ================\n")
  cat(sprintf("95%% Predictive Interval Coverage: %.2f%%\n", coverage_95))
  cat(sprintf("50%% Predictive Interval Coverage: %.2f%%\n", coverage_50))
  #cat(sprintf("Adjusted Bayesian p-value (Freeman-Tukey): %.3f\n", b_p_value))
  cat("========================================================\n\n")
  
  # --- 6. VISUALIZATION: SINGLE TIMESERIES GENERATOR ---
  plot_ppc_timeseries <- function(target_year, target_station = 1, year_label = NULL) {
    year_data <- pointwise_summary %>%
      filter(Year == target_year, Station == target_station)
    
    if (nrow(year_data) == 0) {
      return(NULL)
    }
    
    display_title <- if(!is.null(year_label)) {
      sprintf("Pointwise PPC: Season %s (Station %d)", 
              as.character(year_label), target_station)
    } else {
      sprintf("Pointwise PPC: Season Index %d (Station %d)", target_year, target_station)
    }
    
    p <- ggplot(year_data, aes(x = Period)) +
      # 95% Predictive Envelope (Light Blue)
      geom_ribbon(aes(ymin = Rep_L95,
                      ymax = Rep_U95, 
                      fill = "95% Pred Interval"), alpha = 0.15) +
      # 50% Predictive Envelope (Medium Blue)
      geom_ribbon(aes(ymin = Rep_L50, ymax = Rep_U50, fill = "50% Pred Interval"), alpha = 0.3) +
      # Modeled Median Path
      geom_line(aes(y = Rep_Median, color = "Median Projection"), linewidth = 1) +
      # Real Observations
      geom_point(aes(y = Obs_n, color = "Observed Counts"), size = 2) +
      scale_fill_manual(values = c("95% Pred Interval" = "skyblue3", 
                                   "50% Pred Interval" = "dodgerblue4")) +
      scale_color_manual(values = c("Median Projection" = "blue", 
                                    "Observed Counts" = "firebrick")) +
      labs(
        title = display_title,
        subtitle = "Observed grey whale counts vs. Posterior Predictive intervals",
        x = "Sequential Observation Period",
        y = "Whale Count / Day",
        fill = "Uncertainty Envelopes",
        color = "Lines"
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(color = "grey35", size = 10)
      )
    
    return(p)
  }
  
  # Peak alignment check
  peak_check <- pointwise_summary %>%
    group_by(Year, Station) %>%
    summarise(
      Obs_Peak_Day = Day[which.max(Obs_n)],
      Rep_Peak_Day = Day[which.max(Rep_Median)],
      Peak_Difference = Obs_Peak_Day - Rep_Peak_Day,
      .groups = "drop"
    )
  
  return(list(
    summary = pointwise_summary,
    peak_alignment = peak_check,
    plot_year = plot_ppc_timeseries,
    y_rep = y_rep,
    y_obs = y_obs,
    n_years = n_year,
    n_stations = n_stat,
    df_matched = df_matched,
    pointwise_summary = pointwise_summary
  ))
}

#' Export All PPC Plots to a Single Multi-Page PDF for Appendix
#'
#' @param ppc_results List output from run_direct_ppc()
#' @param filename String filepath for the generated PDF
#' @param year_labels Vector of strings representing true year/season names
export_appendix_pdf <- function(ppc_results, 
                                filename = "appendix_ppc_plots.pdf", year_labels = NULL) {
  
  cat(sprintf("Opening PDF Graphic Device: %s\n", filename))
  
  # Open landscape letter-sized PDF
  pdf(file = filename, width = 11, height = 8.5)
  
  pages_printed <- 0
  
  for (y in 1:ppc_results$n_years) {
    y_lbl <- if(!is.null(year_labels) && length(year_labels) >= y) year_labels[y] else NULL
    
    for (s in 1:ppc_results$n_stations) {
      
      p <- ppc_results$plot_year(target_year = y, target_station = s, year_label = y_lbl)
      
      if (!is.null(p)) {
        print(p)
        pages_printed <- pages_printed + 1
        cat(sprintf("  Printed Page %d: Season %s, Station %d\n", 
                    pages_printed, ifelse(is.null(y_lbl), as.character(y), y_lbl), s))
      }
    }
  }
  
  dev.off()
  cat(sprintf("\nSUCCESS: Saved %d appendix pages to '%s'.\n", pages_printed, filename))
}

# 1. Execute the function to analyze data and align dimensions
ppc_results <- run_direct_ppc(jags_fit = jags_fit, n = n, day = jags_output$jags.input$jags.data$day)

# 3. Compile every season into a single, landscape-oriented PDF document
export_appendix_pdf(
  ppc_results = ppc_results, 
  filename = "Appendix_A_PPC_Fits.pdf", 
  year_labels = seasons
)

ppc_results$pointwise_summary %>%
  mutate(Season = seasons[Year]) %>%
  filter(Day < 100) -> pointwise_summary

p.ppc <- ggplot(data = pointwise_summary,
                aes(x = Period)) +
  geom_ribbon(aes(ymin = Rep_L95,
                  ymax = Rep_U95, 
                  fill = "95% Pred Interval"), alpha = 0.15) +
  # 50% Predictive Envelope (Medium Blue)
  geom_ribbon(aes(ymin = Rep_L50, 
                  ymax = Rep_U50, 
                  fill = "50% Pred Interval"), alpha = 0.3) +
  # Modeled Median Path
  geom_line(aes(y = Rep_Median, color = "Median Projection"), linewidth = 1) +
  # Real Observations
  geom_point(aes(y = Obs_n, color = "Observed Counts"), size = 0.5) +
  scale_fill_manual(values = c("95% Pred Interval" = "skyblue3", 
                               "50% Pred Interval" = "dodgerblue4")) +
  scale_color_manual(values = c("Median Projection" = "blue", 
                                "Observed Counts" = "firebrick")) +
  facet_wrap(~Season, scales = "free") +
  labs( x = "Sequential Observation Period",
    y = "Whale Count",
    fill = "Uncertainty Envelopes",
    color = "Lines"
  ) +
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(color = "grey35", size = 10))

# 
# 
# 
# run_offset_corrected_ppc <- function(jags_fit, n) {
#   
#   # --- 1. EXTRACT RAW DIMENSIONS ---
#   sims_kappa <- jags_fit$sims.list$kappa  # [MCMC_iter, n_days_kappa, n_stat, n_year]
#   sims_r     <- jags_fit$sims.list$r      # [MCMC_iter]
#   
#   n_sims        <- dim(sims_kappa)[1]
#   n_days_kappa  <- dim(sims_kappa)[2]
#   n_stat        <- dim(n)[2]
#   n_year        <- dim(n)[3]
#   n_days_n      <- dim(n)[1] # Includes Day 1
#   
#   cat(sprintf("Observed data 'n' has %d days (Includes Day 1).\n", n_days_n))
#   cat(sprintf("Monitored 'kappa' has %d days (Starts at Day 2).\n", n_days_kappa))
#   
#   # --- 2. BUILD AN EXPLICIT COORDINATE MAP ---
#   # We loop through years, stations, and days (starting at Day 2)
#   # to safely link coordinates together without index-shifting bugs.
#   matched_list <- list()
#   idx <- 1
#   
#   for (y in 1:n_year) {
#     for (s in 1:n_stat) {
#       for (d in 2:n_days_n) {
#         
#         # Core Adjustment: Day 'd' in data maps to row 'd - 1' in kappa
#         d_kappa <- d - 1
#         
#         # Only pull points where a real observation occurred 
#         # AND JAGS evaluated a tracking window
#         if (!is.na(n[d, s, y]) && n[d, s, y] >= 0 && !is.na(sims_kappa[1, d_kappa, s, y])) {
#           matched_list[[idx]] <- data.frame(
#             Day         = d,
#             D_Kappa     = d_kappa,
#             Station     = s,
#             Year        = y,
#             Obs_n       = n[d, s, y]
#           )
#           idx <- idx + 1
#         }
#       }
#     }
#   }
#   
#   df_matched <- do.call(rbind, matched_list)
#   N_data     <- nrow(df_matched)
#   cat(sprintf("Successfully mapped %d perfect observation-to-model coordinate pairs.\n", N_data))
#   
#   # --- 3. EXTRACT POSTERIOR KAPPA VALUES EFFICIENTLY ---
#   cat("Extracting posterior trajectories across the coordinate map...\n")
#   kappa_matrix <- matrix(NA, nrow = n_sims, ncol = N_data)
#   for (j in 1:N_data) {
#     kappa_matrix[, j] <- sims_kappa[, df_matched$D_Kappa[j], df_matched$Station[j], df_matched$Year[j]]
#   }
#   
#   # --- 4. GENERATES REPLICATED DATASETS & DISCREPANCIES ---
#   cat("Simulating predictive counts and calculating Freeman-Tukey metrics...\n")
#   y_rep <- matrix(NA, nrow = n_sims, ncol = N_data)
#   T_obs <- rep(NA, n_sims)
#   T_rep <- rep(NA, n_sims)
#   
#   for (i in 1:n_sims) {
#     r_val <- sims_r[i]
#     
#     # Generate negative binomial predictive counts
#     y_rep[i, ] <- rnbinom(n = N_data, size = r_val, mu = pmax(kappa_matrix[i, ], 1e-6))
#     
#     # Calculate Freeman-Tukey discrepancy metrics on aligned arrays
#     T_obs[i] <- sum((sqrt(df_matched$Obs_n) - sqrt(kappa_matrix[i, ]))^2)
#     T_rep[i] <- sum((sqrt(y_rep[i, ]) - sqrt(kappa_matrix[i, ]))^2)
#   }
#   
#   b_p_value <- mean(T_rep > T_obs)
#   cat(sprintf("\n>>> Adjusted Bayesian p-value: %.3f <<<\n\n", b_p_value))
#   
#   # --- 5. COMPUTE SUMMARY STATISTICS FOR PLOTTING ---
#   pointwise_summary <- df_matched %>%
#     mutate(
#       Rep_Median = apply(y_rep, 2, median),
#       Rep_L95    = apply(y_rep, 2, quantile, probs = 0.025),
#       Rep_U95    = apply(y_rep, 2, quantile, probs = 0.975),
#       Rep_L50    = apply(y_rep, 2, quantile, probs = 0.25),
#       Rep_U50    = apply(y_rep, 2, quantile, probs = 0.75),
#       In_95_CI   = (Obs_n >= Rep_L95) & (Obs_n <= Rep_U95)
#     )
#   
#   cat(sprintf("95%% Predictive Interval Coverage: %.2f%%\n", mean(pointwise_summary$In_95_CI) * 100))
#   
#   # --- 6. IN-LINE TIMESERIES PLOTTING FUNCTION ---
#   plot_year_results <- function(target_year, target_station = 1) {
#     year_data <- pointwise_summary %>%
#       filter(Year == target_year, Station == target_station)
#     
#     ggplot(year_data, aes(x = Day)) +
#       geom_ribbon(aes(ymin = Rep_L95, ymax = Rep_U95, fill = "95% Pred Interval"), alpha = 0.15) +
#       geom_ribbon(aes(ymin = Rep_L50, ymax = Rep_U50, fill = "50% Pred Interval"), alpha = 0.3) +
#       geom_line(aes(y = Rep_Median, color = "Model Median"), size = 1) +
#       geom_point(aes(y = Obs_n, color = "Observed Data"), size = 2) +
#       scale_fill_manual(values = c("95% Pred Interval" = "skyblue3", "50% Pred Interval" = "dodgerblue4")) +
#       scale_color_manual(values = c("Model Median" = "blue", "Observed Data" = "firebrick")) +
#       labs(
#         title = sprintf("Corrected Pointwise PPC: Year %d (Station %d)", target_year, target_station),
#         x = "Julian Day", y = "Whale Count", fill = "Uncertainty Envelopes", color = "Lines"
#       ) +
#       theme_minimal() +
#       theme(legend.position = "bottom")
#   }
#   
#   return(list(
#     p_value = b_p_value,
#     summary = pointwise_summary,
#     plot_year = plot_year_results
#   ))
# }

# # Run the corrected PPC script
# ppc_output <- run_offset_corrected_ppc(jags_fit = jags_fit, n = n)
# 
# # Plot a year to view your red data points sitting cleanly inside the blue prediction interval
# ppc_output$plot_year(target_year = 1, target_station = 1)
