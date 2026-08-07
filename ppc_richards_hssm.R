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

# ---------------------------------------------------------------------
# 1. Reconstruct kappa and draw replicate datasets
# ---------------------------------------------------------------------
ppc_setup <- function(fit, sd, n_draws = 500, seed = 1) {
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
       draw_idx = idx, sd = sd)
}

# ---------------------------------------------------------------------
# 2. Randomised quantile (Dunn-Smyth) residuals
#    Exactly standard normal under a correct model, even for counts.
#    These are the right tool for detecting residual STRUCTURE.
# ---------------------------------------------------------------------
ds_residuals <- function(ps, draw = 1) {
  y <- ps$y; k <- ps$kappa[draw, ]
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
# 4. THE reviewer-#2 check: within-season residual autocorrelation.
#    A curve too rigid to track the migration leaves runs of same-signed
#    residuals along the season. Aggregate to daily means first, since
#    multiple stations/watch periods share a day.
# ---------------------------------------------------------------------
ppc_autocorr <- function(ps, n_rep = 200) {
  sd <- ps$sd
  acf1 <- function(r) {
    df <- data.frame(y = sd$year_idx, d = sd$day_idx, r = r) |>
      group_by(y, d) |> summarise(r = mean(r), .groups = "drop") |>
      arrange(y, d)
    # lag-1 correlation within each season, averaged over seasons
    vals <- df |> group_by(y) |>
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
ppc_richards <- function(fit, sd, n_draws = 500, seed = 1) {
  ps <- ppc_setup(fit, sd, n_draws, seed)
  list(setup    = ps,
       pvalues  = ppc_pvalues(ps),
       autocorr = ppc_autocorr(ps),
       plots    = ppc_plots(ps))
}
