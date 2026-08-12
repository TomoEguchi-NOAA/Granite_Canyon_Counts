# =====================================================================
# Posterior predictive checks for the Richards HSSM gray whale model
#
# Works directly from an existing cmdstanr fit -- no refit needed.
# kappa is reconstructed from saved draws of log_N_latent (a transformed
# parameter, saved by default), alpha, BF_Fixed, VS_Fixed and phi:
#
#   log kappa[i] = log_N_latent[day[i], year[i]]
#                + log(watch_length[i])
#                + log inv_logit(alpha[obs[i]] + BF*bf[i] + VS*vs[i])
#
# Usage:
#   res <- ppc_richards(fit_stan, stan.data$stan.data, n_draws = 500)
#   res$pvalues          # table of Bayesian p-values
#   res$plots$resid_day  # etc.
# =====================================================================

library(posterior)
library(ggplot2)
library(dplyr)
library(tidyr)

# sd: Stan data
# ps: output of ppc_setup

# ---------------------------------------------------------------------
# 1. Reconstruct kappa and draw replicate datasets
# ---------------------------------------------------------------------
ppc_setup <- function(fit, sd, n_draws = 1500, seed = 1) {
  set.seed(seed)
  N   <- sd$N_flat
  dm  <- fit$draws(format = "draws_matrix")
  idx <- sort(sample(nrow(dm), min(n_draws, nrow(dm))))

  # Pull only the log_N_latent cells actually used by an observation.
  cols_N <- sprintf("log_N_latent[%d,%d]", sd$day_idx, sd$year_idx)
  cols_a <- sprintf("alpha[%d]", sd$observer_idx)
  stopifnot(all(cols_N %in% colnames(dm)), all(unique(cols_a) %in% colnames(dm)))

  logN  <- dm[idx, cols_N, drop = FALSE]           # n_draws x N
  alpha <- dm[idx, cols_a, drop = FALSE]           # n_draws x N
  BF    <- as.numeric(dm[idx, "BF_Fixed"])
  VS    <- as.numeric(dm[idx, "VS_Fixed"])

  eta <- alpha +
    outer(BF, sd$bf) +
    outer(VS, sd$vs)
  log_p <- -log1p(exp(-eta))                       # log inv_logit, stable
  log_kappa <- logN + log_p +
    matrix(log(sd$watch_length), nrow = length(idx), ncol = N, byrow = TRUE)
  kappa <- exp(log_kappa)

  is_nb <- "phi[1]" %in% colnames(dm)
  phi   <- if (is_nb) as.numeric(dm[idx, "phi[1]"]) else NA_real_

  zi_mode <- sd$zi_mode
  if (zi_mode > 0) {
    za <- as.numeric(dm[idx, "zi_a[1]"])
    zb <- if (zi_mode == 2) as.numeric(dm[idx, "zi_b[1]"]) else rep(0, length(idx))
  } else {
    za = zb = NA
  }
  
  # Replicate datasets
  y_rep <- matrix(NA_integer_, nrow = length(idx), ncol = N)
  for (d in seq_along(idx)) {
    y <- if (is_nb) rnbinom(N, mu = kappa[d, ], size = phi[d])
    else       rpois(N, lambda = kappa[d, ])
    if (zi_mode > 0) {
      lk   <- log(kappa[d, ])
      lpi  <- za[d] + zb[d] * (lk - mean(lk))
      y[runif(N) < plogis(lpi)] <- 0L        # structural zeros
    }
    y_rep[d, ] <- y
  }
  
  list(y = sd$n, y_rep = y_rep, kappa = kappa, phi = phi, is_nb = is_nb,
       draw_idx = idx, zi_mode = zi_mode, za = za, zb = zb, sd = sd)
}

# ---------------------------------------------------------------------
# 2. Randomised quantile (Dunn-Smyth) residuals
#    Exactly standard normal under a correct model, even for counts.
#    These are the right tool for detecting residual STRUCTURE.
# ---------------------------------------------------------------------
ds_residuals <- function(ps, draw = 1) {
  y <- ps$y; 
  k <- ps$kappa[draw, ]
  if (ps$is_nb) {
    lo <- pnbinom(y - 1, mu = k, size = ps$phi[draw])
    hi <- pnbinom(y,     mu = k, size = ps$phi[draw])
  } else {
    lo <- ppois(y - 1, lambda = k); hi <- ppois(y, lambda = k)
  }
  if (ps$zi_mode > 0) {
    pi_ <- plogis(ps$za[draw] + ps$zb[draw] * (log(k) - mean(log(k))))
    lo  <- ifelse(y == 0, 0, pi_ + (1 - pi_) * lo)   # mass below y
    hi  <- pi_ + (1 - pi_) * hi
  }
  qnorm(pmin(pmax(runif(length(y), lo, hi), 1e-10), 1 - 1e-10))
}

# The following was before ZI models were introduced
# ds_residuals <- function(ps, draw = 1) {
#   y <- ps$y; k <- ps$kappa[draw, ]
#   if (ps$is_nb) {
#     lo <- pnbinom(y - 1, mu = k, size = ps$phi[draw])
#     hi <- pnbinom(y,     mu = k, size = ps$phi[draw])
#   } else {
#     lo <- ppois(y - 1, lambda = k)
#     hi <- ppois(y,     lambda = k)
#   }
#   qnorm(pmin(pmax(runif(length(y), lo, hi), 1e-10), 1 - 1e-10))
# }

# ---------------------------------------------------------------------
# 3. Discrepancy statistics + Bayesian p-values
#    p near 0 or 1 indicates the model fails to reproduce that feature.
# ---------------------------------------------------------------------
ppc_pvalues <- function(ps) {
  stats <- list(
    prop_zero  = function(v) mean(v == 0),
    maximum    = function(v) max(v),
    sd_counts  = function(v) sd(v),
    mean_count = function(v) mean(v),
    q95        = function(v) as.numeric(quantile(v, 0.95)),
    n_over_100 = function(v) sum(v > 100)
  )
  out <- lapply(names(stats), function(nm) {
    f <- stats[[nm]]
    T_obs <- f(ps$y)
    T_rep <- apply(ps$y_rep, 1, f)
    data.frame(statistic = nm, observed = T_obs,
               rep_median = median(T_rep),
               rep_lo = quantile(T_rep, .025), rep_hi = quantile(T_rep, .975),
               p_value = mean(T_rep >= T_obs))
  })
  res <- do.call(rbind, out); rownames(res) <- NULL
  res$flag <- ifelse(res$p_value < 0.025 | res$p_value > 0.975, "**", "")
  res
}

# ---------------------------------------------------------------------
# 4. THE comment-#2 check: within-season residual autocorrelation.
#    A curve too rigid to track the migration leaves runs of same-signed
#    residuals along the season. Aggregate to daily means first, since
#    multiple stations/watch periods share a day.
# ---------------------------------------------------------------------
ppc_autocorr <- function(ps, n_rep = 1000) {
  sd <- ps$sd
  acf1 <- function(r) {
    df <- data.frame(y = sd$year_idx, d = sd$day_idx, r = r) |>
      group_by(y, d) |> 
      summarise(r = mean(r), .groups = "drop") |>
      arrange(y, d)
    
    # lag-1 correlation within each season, averaged over seasons
    vals <- df |> 
      group_by(y) |>
      summarise(a = if (n() > 5) cor(r[-n()], r[-1]) else NA_real_,
                .groups = "drop")
    mean(vals$a, na.rm = TRUE)
  }
  
  obs <- mean(replicate(20, acf1(ds_residuals(ps, draw = 1))))
  
  rep <- sapply(seq_len(min(n_rep, nrow(ps$y_rep))), function(d) {
    yr <- ps$y_rep[d, ]
    k  <- ps$kappa[d, ]
    lo <- if (ps$is_nb) pnbinom(yr - 1, mu = k, size = ps$phi[d]) else ppois(yr - 1, k)
    hi <- if (ps$is_nb) pnbinom(yr,     mu = k, size = ps$phi[d]) else ppois(yr,     k)
    acf1(qnorm(pmin(pmax(runif(length(yr), lo, hi), 1e-10), 1 - 1e-10)))
  })
  
  list(observed = obs, rep_median = median(rep),
       rep_interval = quantile(rep, c(.025, .975)),
       p_value = mean(rep >= obs))
}

# ---------------------------------------------------------------------
# 5. Plots
# ---------------------------------------------------------------------
ppc_plots <- function(ps) {
  sd <- ps$sd
  p_lo <- apply(ps$y_rep, 2, quantile, .025)
  p_hi <- apply(ps$y_rep, 2, quantile, .975)
  p_md <- apply(ps$y_rep, 2, median)
  cover <- mean(sd$n >= p_lo & sd$n <= p_hi)

  df <- data.frame(obs = sd$n, med = p_md, lo = p_lo, hi = p_hi,
                   day = sd$day_idx, year = factor(sd$year_idx),
                   resid = ds_residuals(ps, 1))

  list(
    coverage = cover,

    # (a) observed vs predicted with 95% predictive band
    obs_pred = ggplot(df, aes(med, obs)) +
      geom_linerange(aes(xmin = lo, xmax = hi), colour = "grey75", alpha = .4) +
      geom_point(size = .6, alpha = .5) +
      geom_abline(slope = 1, intercept = 0, colour = "red") +
      labs(x = "Posterior predictive median", y = "Observed count",
           subtitle = sprintf("95%% predictive interval coverage = %.1f%% (nominal 95%%)",
                              100 * cover)) + theme_bw(),

    # (b) DS residuals vs day of season -- systematic curve misfit
    resid_day = ggplot(df, aes(day, resid)) +
      geom_hline(yintercept = 0, colour = "red") +
      geom_point(size = .5, alpha = .3) +
      geom_smooth(method = "loess", se = TRUE, span = .4) +
      labs(x = "Day of season", y = "Randomised quantile residual",
           subtitle = "Flat loess = Richards curve adequate; systematic curvature = too rigid") +
      theme_bw(),

    # (c) residuals vs fitted -- mean/variance adequacy
    resid_fit = ggplot(df, aes(med, resid)) +
      geom_hline(yintercept = 0, colour = "red") +
      geom_point(size = .5, alpha = .3) +
      geom_smooth(method = "loess", se = TRUE) +
      scale_x_continuous(trans = "log1p") +
      labs(x = "Predicted count (log1p scale)", y = "Residual") + theme_bw(),

    # (d) density overlay
    dens = ggplot() +
      geom_density(data = data.frame(v = as.vector(ps$y_rep[1:min(50, nrow(ps$y_rep)), ]),
                                     g = "replicate"),
                   aes(v, group = g), colour = "grey70", alpha = .3) +
      geom_density(data = data.frame(v = sd$n), aes(v), colour = "red", linewidth = 1) +
      scale_x_continuous(trans = "log1p") +
      labs(x = "Count (log1p scale)", y = "Density",
           subtitle = "Red = observed, grey = replicates") + theme_bw(),

    # (e) per-season residual QQ, to find seasons that fit poorly
    resid_season = ggplot(df, aes(sample = resid)) +
      stat_qq(size = .4) + stat_qq_line(colour = "red") +
      facet_wrap(~ year, scales = "free") +
      labs(subtitle = "Per-season QQ of randomised quantile residuals") +
      theme_bw(base_size = 7)
  )
}

# ---------------------------------------------------------------------
# 6. Wrapper
# ---------------------------------------------------------------------
ppc_richards <- function(fit, sd, n_draws = 1500, seed = 1) {
  ps <- ppc_setup(fit, sd, n_draws, seed)
  list(setup    = ps,
       pvalues  = ppc_pvalues(ps),
       autocorr = ppc_autocorr(ps),
       plots    = ppc_plots(ps))
}

# The following functions are from Gemini
# Assuming 'stan_data' is the list you currently pass to the rstan::sampling() function
# --- How to use this ---
# 1. Generate the gap data
# gap_simulation <- simulate_weather_gaps(my_stan_data, gap_length = 5)
#
# 2. Fit your top models to this new data
# fit_gapped <- sampling(my_stan_model, data = gap_simulation$data, ...)
#
# 3. Compare the 'Raw_Est' from the gapped fit to the 'Raw_Est' from the full data fit.
# The model whose Raw_Est shifts the least is the most robust interpolator.
# 
# Assuming 'stan_data' is the list you currently pass to the rstan::sampling() function

simulate_weather_gaps <- function(stan_data, gap_length = 4) {
  # Make a copy of the original data to modify
  gapped_data <- stan_data
  
  # Track which flat indices we are removing
  rows_to_remove <- integer(0)
  
  # Loop through each year to inject a continuous gap
  for (y in 1:stan_data$n_year) {
    
    # Find all flat indices corresponding to observations in year 'y'
    year_rows <- which(stan_data$year_idx == y)
    
    if (length(year_rows) > 0) {
      # Get the unique days actually observed this year
      days_in_year <- sort(unique(stan_data$day_idx[year_rows]))
      
      # Ensure we don't pick a start day that pushes the gap beyond our observed data
      valid_start_days <- days_in_year[days_in_year <= (max(days_in_year) - gap_length)]
      
      if (length(valid_start_days) > 0) {
        # Randomly select a starting day for the gap
        gap_start <- sample(valid_start_days, 1)
        
        # Define the contiguous block of days to remove
        gap_days <- gap_start:(gap_start + gap_length - 1)
        
        # Identify the exact row indices in the flat arrays for these days
        gap_indices <- year_rows[stan_data$day_idx[year_rows] %in% gap_days]
        rows_to_remove <- c(rows_to_remove, gap_indices)
      }
    }
  }
  
  # Remove these rows from all flat arrays expected by your Stan model
  gapped_data$n            <- stan_data$n[-rows_to_remove]
  gapped_data$bf           <- stan_data$bf[-rows_to_remove]
  gapped_data$vs           <- stan_data$vs[-rows_to_remove]
  gapped_data$observer_idx <- stan_data$observer_idx[-rows_to_remove]
  gapped_data$watch_length <- stan_data$watch_length[-rows_to_remove]
  gapped_data$day_idx      <- stan_data$day_idx[-rows_to_remove]
  gapped_data$year_idx     <- stan_data$year_idx[-rows_to_remove]
  
  # Update the total flat observation count
  gapped_data$N_flat <- stan_data$N_flat - length(rows_to_remove)
  
  # Return both the new data list and the removed indices (for validation later)
  return(list(
    data = gapped_data, 
    removed_indices = rows_to_remove
  ))
}

# --- How to use this ---
# 1. Generate the gap data
# gap_simulation <- simulate_weather_gaps(my_stan_data, gap_length = 5)
#
# 2. Fit your top models to this new data
# fit_gapped <- sampling(my_stan_model, data = gap_simulation$data, ...)
#
# 3. Compare the 'Raw_Est' from the gapped fit to the 'Raw_Est' from the full data fit.
# The model whose Raw_Est shifts the least is the most robust interpolator.

# A function to remove specific days from a specific year in your flattened Stan data
simulate_targeted_gap <- function(stan_data, target_year, target_days) {
  
  # Make a copy of the original data to modify
  gapped_data <- stan_data
  
  # Find the exact row indices in the flat arrays where:
  # 1. The year matches the target_year
  # 2. The day is within our target_days vector
  rows_to_remove <- which(stan_data$year_idx == target_year & 
                            stan_data$day_idx %in% target_days)
  
  # Safety check: Did we actually find any rows to remove?
  if (length(rows_to_remove) == 0) {
    warning("No observations found for the specified year and days. Returning original data.")
    return(list(data = stan_data, removed_indices = integer(0)))
  }
  
  # Remove these rows from all flat arrays expected by your Stan model
  gapped_data$n            <- stan_data$n[-rows_to_remove]
  gapped_data$bf           <- stan_data$bf[-rows_to_remove]
  gapped_data$vs           <- stan_data$vs[-rows_to_remove]
  gapped_data$observer_idx <- stan_data$observer_idx[-rows_to_remove]
  gapped_data$watch_length <- stan_data$watch_length[-rows_to_remove]
  gapped_data$day_idx      <- stan_data$day_idx[-rows_to_remove]
  gapped_data$year_idx     <- stan_data$year_idx[-rows_to_remove]
  
  # Update the total flat observation count
  gapped_data$N_flat <- stan_data$N_flat - length(rows_to_remove)
  
  # Print a helpful summary to the console
  message(sprintf("Removed %d observation periods from Year %d (Days: %s)", 
                  length(rows_to_remove), 
                  target_year, 
                  paste(target_days, collapse = ", ")))
  
  # Return both the new data list and the removed indices
  return(list(
    data = gapped_data, 
    removed_indices = rows_to_remove
  ))
}


library(ggplot2)
library(dplyr)

# Assuming 'y_obs' is your original vector of counts (stan_data$n)
# Assuming 'y_rep' is a matrix of posterior predictive draws for those counts
# (You would need to generate y_rep in your generated quantities block, 
# e.g., y_rep[i] = neg_binomial_2_log_rng(log_kappa[i], phi[1]))

plot_rootogram <- function(y_obs, y_rep) {
  
  # Calculate observed frequencies for each count (0, 1, 2, ...)
  max_count <- max(y_obs, y_rep)
  obs_freq <- tabulate(y_obs + 1, nbins = max_count + 1)
  
  # Calculate expected frequencies from the posterior draws
  # Apply tabulate to each draw, then take the mean across all draws
  rep_freq_matrix <- apply(y_rep, 1, function(row) tabulate(row + 1, nbins = max_count + 1))
  exp_freq <- rowMeans(rep_freq_matrix)
  
  # Create a dataframe for plotting
  df <- data.frame(
    count = 0:max_count,
    obs_freq = obs_freq,
    exp_freq = exp_freq
  ) %>%
    # Filter out extremely long tails with zero observations and zero expectations
    filter(obs_freq > 0 | exp_freq > 0.1)
  
  # Calculate square roots for the rootogram
  df$sqrt_obs <- sqrt(df$obs_freq)
  df$sqrt_exp <- sqrt(df$exp_freq)
  
  # For a "hanging" rootogram, the bars hang from the expected curve
  df$bar_bottom <- df$sqrt_exp - df$sqrt_obs
  
  ggplot(df, aes(x = count)) +
    geom_rect(aes(xmin = count - 0.4, xmax = count + 0.4, 
                  ymin = bar_bottom, ymax = sqrt_exp), 
              fill = "lightgray", color = "black") +
    geom_line(aes(y = sqrt_exp), color = "red", linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      title = "Hanging Rootogram",
      x = "Whale Count",
      y = "Sqrt(Frequency)"
    ) +
    theme_minimal()
}

# plot_rootogram(stan_data$n, posterior_predictive_matrix)

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

