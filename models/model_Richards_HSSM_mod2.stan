// =====================================================================
// Richards hierarchical state-space model, ENP gray whale abundance,
// Granite Canyon CA.  Stan port of model_Richards_HSSM_M1a2.jags
// M1 = S1 and S2 both season-specific;  a2 = Negative Binomial
//
// Verified line-by-line against the JAGS source. In PARITY MODE
// (use_process_error = 0, use_shape_dev = 0, anchor_sd = 1e-6) this is
// algebraically identical to the JAGS model, with two exceptions noted
// at the bottom of this header.
//
// JAGS -> Stan notes:
//  * dnorm(mu, tau) is PRECISION; Stan normal(mu, sigma) is SD.
//    sd = 1/sqrt(tau).  All conversions passed as data, see R list below.
//  * dnegbin(p, r) with p = r/(r+kappa) has mean kappa, var kappa+kappa^2/r,
//    which is exactly neg_binomial_2(kappa, r). So phi == r directly.
//  * The Richards curve is computed in log space with log1p_exp() rather
//    than by forming exp() and clamping. Numerically identical, but
//    differentiable everywhere, which HMC requires and Gibbs does not.
//
// TWO DELIBERATE DEPARTURES FROM JAGS (both switchable):
//  1. anchor_sd lets uncertainty in Durban's baseline p propagate into N.
//     JAGS pins it at exactly logit(0.80) = 1.39 with zero uncertainty.
//     Set anchor_sd = 0.01 (1e-6) for parity; set it to the SE of logit(0.80)
//     to propagate properly.
//  2. corr.factor: JAGS draws ONE value per iteration shared across all
//     years (perfectly correlated). Reproduced here. Set
//     independent_corr = 1 for an independent draw per year.
// =====================================================================

data {
  int<lower=1> n_year;
  int<lower=1> n_days;
  int<lower=1> n_observer;          // n.obs.fixed
  int<lower=1> N_flat;              // total watch periods, d >= 2 only

  array[N_flat] int<lower=0> n;     // observed counts
  vector[N_flat] bf;                // RAW Beaufort  (JAGS bf.1)
  vector[N_flat] vs;                // RAW visibility MINUS 1 (JAGS vs.1 - 1)
  array[N_flat] int<lower=1> observer_idx;   // obs.fixed
  vector<lower=0>[N_flat] watch_length;      // in days
  vector[n_year] year_values;                // year.index, centred

  array[N_flat] int<lower=1> day_idx;        // JAGS day[d,s,y], must be >= 2
  array[N_flat] int<lower=1> year_idx;

  // ---- Prior scales, all on the SD scale ------------------------------
  real<lower=0> sd_beta1_P;         // 2.0     <- dnorm(0, 0.25)
  real beta0_Max_mu;                // 7.6
  real<lower=0> sd_beta0_Max;       // 2.0     <- dnorm(7.6, 0.25)
  real<lower=0> sd_beta1_Max;       // 2.0     <- dnorm(0, 0.25)
  real<lower=0> sd_BF;              // 2.0     <- dnorm(0, 0.25)
  real<lower=0> sd_VS;              // 2.0     <- dnorm(0, 0.25)
  real<lower=0> sd_sigma_proc_P;    // 5.0     <- dnorm(0, 1/5^2)T(0,)
  real<lower=0> sd_sigma_proc_Max;  // 5.0     <- dnorm(0, 1/5^2)T(0,)

  real alpha_S_mu;                  // 10.0    <- dnorm(10, 0.1)T(0,)
  real<lower=0> alpha_S_sd;         // sqrt(10) = 3.1623
  real<lower=0> beta_S_shape;       // 1.0     <- dgamma(1, 1)
  real<lower=0> beta_S_rate;        // 1.0

  // Direct prior for a CONSTANT S1/S2 (used only when S*_by_season = 0).
  // Retaining the gamma hyperpriors over a single draw implies a prior with
  // no finite mean (beta ~ Gamma(1,1) => E[1/beta] diverges), so constant-S
  // models get a weakly informative gamma instead. Gamma(10, 1): median ~9.7,
  // 90% interval ~[5.4, 15.7], comfortably covering the posterior region.
  real<lower=0> S_const_shape;      // e.g. 10.0
  real<lower=0> S_const_rate;       // e.g. 1.0

  real<lower=0> sigma_Obs_max;      // 1.5     <- dunif(0, 1.5)
  real<lower=0> phi_max;            // 50.0    <- dunif(0, 50)

  // ---- Detection anchor -------------------------------------------------
  real anchor_mu;                   // 1.39 = logit(0.80)
  real<lower=0> anchor_sd;          // 1e-6 for JAGS parity; else SE

  real<lower=0> boundary_N;         // 0.0001, JAGS mean.N[1,y] and [n.days,y]

  // ---- Switches (all 0 = parity mode) ----------------------------------
  int<lower=0, upper=1> use_process_error;
  int<lower=0, upper=1> use_shape_dev;
  int<lower=0, upper=1> independent_corr;
  // Centred vs non-centred season effects. Use CENTRED (1) when the
  // process SDs are well away from zero and the per-season data are
  // informative -- which is this model. Non-centred (0) only helps when
  // the SD can approach zero.
  int<lower=0, upper=1> centred_P;
  int<lower=0, upper=1> centred_Max;

  // ---- Structural sensitivity switches ---------------------------------
  // use_trend_*  = 0 drops the linear-in-year term entirely (beta1 removed
  //                from the model, not merely initialised at zero).
  // use_pooling_Max = 0 removes hierarchical shrinkage on the season levels
  //                altogether: each log_Max[y] gets an independent diffuse
  //                prior, as in the calf-production model. Use this to test
  //                whether shrinkage toward the trend inflates low seasons.
  //                When 0, beta0_Max / beta1_Max / sigma_proc_Max are not
  //                identified and simply sample their priors -- ignore them.
  // ---- Plateau (P1 != P2) ----------------------------------------------
  // delta = plateau width in days. The descent limb centres at P + delta/2
  // and the ascent limb at P - delta/2, so delta = 0 reproduces the single-P
  // model exactly and P remains the plateau MIDPOINT (peak timing, and its
  // trend, keep their interpretation). delta >= 0 by construction, so there
  // is no P1/P2 ridge and no label switching.
  // ---- MODEL STRUCTURE (Table 1) ---------------------------------------
  //   S1_by_season  S2_by_season  likelihood_NB      model
  //        1             1              0/1          M1a1 / M1a2
  //        0             0              0/1          M2a1 / M2a2
  //        1             0              0/1          M3a1 / M3a2
  //        0             1              0/1          M4a1 / M4a2
  int<lower=0, upper=1> S1_by_season;
  int<lower=0, upper=1> S2_by_season;
  int<lower=0, upper=1> likelihood_NB;   // 0 = Poisson, 1 = Negative Binomial

  int<lower=0, upper=1> use_plateau;
  int<lower=0, upper=1> plateau_by_year;   // 0 = one shared delta
  real<lower=0> sd_delta;                  // half-normal scale, e.g. 5 days

  int<lower=0, upper=1> use_trend_P;
  int<lower=0, upper=1> use_trend_Max;
  int<lower=0, upper=1> use_pooling_Max;
  real<lower=0> sd_unpooled_Max;   // e.g. 5.0, on the log scale

  // ---- Shared periodic deviation, only if use_shape_dev = 1 ------------
  int<lower=1> n_period;
  array[n_days] int<lower=1> period_idx;
  real<lower=0> sd_sigma_shape;
}

transformed data {
  array[N_flat] int flat_idx;                       // column-major into matrix
  vector[N_flat] log_watch_length = log(watch_length);
  real log_gamma_term = log(2 * exp(1.0) - 1);      // JAGS hard-codes exp(1)
  real inv_exp_gamma  = 1.0 / exp(1.0);
  real log_boundary_N = log(boundary_N);

  for (i in 1:N_flat)
    flat_idx[i] = (year_idx[i] - 1) * n_days + day_idx[i];
}

parameters {
  vector<lower=0>[n_year] S1;       // untruncated above, as in JAGS
  vector<lower=0>[n_year] S2;
  // Gamma hyperparameters sampled as (mean, shape) rather than (shape, rate):
  // alpha and beta are identified mainly through their ratio and trace a long
  // curved ridge. The original JAGS prior is preserved exactly via the
  // Jacobian in the model block, so this is a reparameterisation only.
  // Hyperparameters exist only when there is a group of seasons to inform them.
  array[S1_by_season] real<lower=0> mu_S1;
  array[S1_by_season] real<lower=0> shape_S1;
  array[S2_by_season] real<lower=0> mu_S2;
  array[S2_by_season] real<lower=0> shape_S2;

  vector[n_observer] alpha;                    // cell means, NO global intercept
  real<lower=0, upper=sigma_Obs_max> sigma_Obs;   // dunif(0, sigma_Obs_max)
  real logit_p0;                               // anchor, ~fixed when anchor_sd tiny

  real BF_Fixed;
  real VS_Fixed;

  real<lower=30, upper=60> beta0_P;            // dunif(30, 60)
  array[use_trend_P] real beta1_P;             // absent when use_trend_P = 0
  real<lower=0> sigma_proc_P;
  vector[n_year] P_raw;

  real beta0_Max;
  array[use_trend_Max] real beta1_Max;         // absent when use_trend_Max = 0
  real<lower=0> sigma_proc_Max;
  vector[n_year] log_Max_raw;

  // Dispersion exists only under the Negative Binomial likelihood.
  vector<lower=0, upper=phi_max>[likelihood_NB] phi;   // dunif(0, phi_max)

  matrix[use_process_error ? n_days : 0,
         use_process_error ? n_year : 0] log_N_raw;
  array[use_process_error] real<lower=0> sigma_process;

  // Plateau width: length 0 (off), 1 (shared), or n_year (season-specific)
  vector<lower=0>[use_plateau ? (plateau_by_year ? n_year : 1) : 0] delta;

  sum_to_zero_vector[n_period] shape_raw;
  real<lower=0> sigma_shape;
}

transformed parameters {
  // Reported on the original JAGS scale for comparison. Zero-sized, hence
  // absent from output, when the corresponding S is constant.
  array[S1_by_season] real alpha_S1 = rep_array(0.0, S1_by_season);
  array[S1_by_season] real beta_S1  = rep_array(0.0, S1_by_season);
  array[S2_by_season] real alpha_S2 = rep_array(0.0, S2_by_season);
  array[S2_by_season] real beta_S2  = rep_array(0.0, S2_by_season);
  vector[n_year] P;
  vector[n_year] log_Max;
  vector[n_year] S1_y = S1_by_season ? S1 : rep_vector(S1[1], n_year);
  vector[n_year] S2_y = S2_by_season ? S2 : rep_vector(S2[1], n_year);
  vector[n_year] delta_y = rep_vector(0.0, n_year);
  real b1_P   = 0.0;
  real b1_Max = 0.0;
  vector[n_year] mu_P;
  vector[n_year] mu_Max;
  matrix[n_days, n_year] log_N_latent;

  if (S1_by_season) { alpha_S1[1] = shape_S1[1]; beta_S1[1] = shape_S1[1] / mu_S1[1]; }
  if (S2_by_season) { alpha_S2[1] = shape_S2[1]; beta_S2[1] = shape_S2[1] / mu_S2[1]; }

  if (use_plateau)
    delta_y = plateau_by_year ? delta : rep_vector(delta[1], n_year);
  if (use_trend_P)   b1_P   = beta1_P[1];
  if (use_trend_Max) b1_Max = beta1_Max[1];
  mu_P   = beta0_P   + b1_P   * year_values;
  mu_Max = beta0_Max + b1_Max * year_values;

  // In centred mode the raw vector IS the season effect; in non-centred
  // mode it is a standard normal that gets shifted and scaled.
  P       = centred_P   ? P_raw       : mu_P   + sigma_proc_P   * P_raw;
  log_Max = centred_Max ? log_Max_raw : mu_Max + sigma_proc_Max * log_Max_raw;

  {
    vector[n_period] shape_dev = use_shape_dev ? sigma_shape * shape_raw
                                               : rep_vector(0.0, n_period);
    for (y in 1:n_year) {
      real inv_S1 = 1.0 / S1_y[y];
      real inv_S2 = 1.0 / S2_y[y];

      // JAGS fixes the first and last day; no observations fall on them.
      log_N_latent[1, y]      = log_boundary_N;
      log_N_latent[n_days, y] = log_boundary_N;

      for (t in 2:(n_days - 1)) {
        // Descent limb centred later, ascent limb earlier, by delta/2 each.
        real p_minus_t1 = P[y] + 0.5 * delta_y[y] - t;   // C1, descent (S1)
        real p_minus_t2 = P[y] - 0.5 * delta_y[y] - t;   // C2, ascent  (S2)
        // C1 = [1 + (2e - 1) exp(-(P-t)/S1)]^(-1/e)
        // C2 = [1 + (2e - 1) exp(+(P-t)/S2)]^(-1/e)
        real log_C1 = -inv_exp_gamma * log1p_exp(log_gamma_term - inv_S1 * p_minus_t1);
        real log_C2 = -inv_exp_gamma * log1p_exp(log_gamma_term + inv_S2 * p_minus_t2);

        log_N_latent[t, y] = log_Max[y] + log_C1 + log_C2 + shape_dev[period_idx[t]];

        if (use_process_error)
          log_N_latent[t, y] += sigma_process[1] * log_N_raw[t, y];
      }
    }
  }
}

model {
  // ---- Richards trend priors -------------------------------------------
  // beta0_P is dunif(30,60) via its declared bounds.
  if (use_trend_P) beta1_P ~ normal(0, sd_beta1_P);
  sigma_proc_P   ~ normal(0, sd_sigma_proc_P);      // half-normal via <lower=0>
  if (centred_P)  P_raw ~ normal(mu_P, sigma_proc_P);
  else            P_raw ~ std_normal();

  beta0_Max      ~ normal(beta0_Max_mu, sd_beta0_Max);
  if (use_trend_Max) beta1_Max ~ normal(0, sd_beta1_Max);
  sigma_proc_Max ~ normal(0, sd_sigma_proc_Max);
  if (!use_pooling_Max)
    log_Max_raw ~ normal(beta0_Max_mu, sd_unpooled_Max);   // independent levels
  else if (centred_Max)
    log_Max_raw ~ normal(mu_Max, sigma_proc_Max);
  else
    log_Max_raw ~ std_normal();

  // ---- Shape parameters, hyperparameters estimated ----------------------
  if (S1_by_season) {
    S1          ~ gamma(alpha_S1[1], beta_S1[1]);
    alpha_S1[1] ~ normal(alpha_S_mu, alpha_S_sd);   // dnorm(10, 0.1)T(0,)
    beta_S1[1]  ~ gamma(beta_S_shape, beta_S_rate); // dgamma(1, 1)
    target += log(shape_S1[1]) - 2 * log(mu_S1[1]); // Jacobian, (alpha,beta)->(shape,mu)
  } else {
    S1 ~ gamma(S_const_shape, S_const_rate);        // direct, weakly informative
  }
  if (S2_by_season) {
    S2          ~ gamma(alpha_S2[1], beta_S2[1]);
    alpha_S2[1] ~ normal(alpha_S_mu, alpha_S_sd);
    beta_S2[1]  ~ gamma(beta_S_shape, beta_S_rate);
    target += log(shape_S2[1]) - 2 * log(mu_S2[1]);
  } else {
    S2 ~ gamma(S_const_shape, S_const_rate);
  }

  // ---- Detection --------------------------------------------------------
  logit_p0 ~ normal(anchor_mu, anchor_sd);
  alpha    ~ normal(logit_p0, sigma_Obs);           // JAGS: dnorm(1.39, tau.obs)
  // sigma_Obs is dunif(0, sigma_Obs_max) via its declared bounds.
  BF_Fixed ~ normal(0, sd_BF);
  VS_Fixed ~ normal(0, sd_VS);

  // phi is dunif(0, phi_max) via its declared bounds.

  if (use_process_error) {
    sigma_process[1]     ~ normal(0, 1);
    to_vector(log_N_raw) ~ std_normal();
  }
  if (use_plateau) delta ~ normal(0, sd_delta);   // half-normal via <lower=0>
  shape_raw   ~ std_normal();
  sigma_shape ~ normal(0, sd_sigma_shape);

  // ---- Observation model -------------------------------------------------
  {
    vector[N_flat] eta = alpha[observer_idx] + BF_Fixed * bf + VS_Fixed * vs;
    vector[N_flat] log_kappa = to_vector(log_N_latent)[flat_idx]
                             + log_watch_length
                             + log_inv_logit(eta);
    if (likelihood_NB) n ~ neg_binomial_2_log(log_kappa, phi[1]);
    else               n ~ poisson_log(log_kappa);
  }
}

generated quantities {
  vector[n_year] Raw_Est;
  vector[n_year] Corrected_Est;

  // ---- Realised curve shape, per season ---------------------------------
  // P is the MIDPOINT of the two limb centres, not the mode: whenever
  // S1 != S2 the actual peak sits several days away from P. These
  // quantities describe the fitted curve itself, so peak_day_slope is the
  // trend in observed peak timing and is the number to report alongside
  // (or instead of) beta1_P.
  vector[n_year] peak_day;      // mode, sub-day via parabolic interpolation
  vector[n_year] peak_height;   // N at the mode
  vector[n_year] peak_width;    // days with N above half of the seasonal peak
  real peak_day_slope;          // days per year, OLS through peak_day
  real peak_day_decade;         // days per decade
  vector[N_flat] log_lik;
  vector[n_year] Max = exp(log_Max);
  real p0 = inv_logit(logit_p0);
  vector[n_observer] obs_p = inv_logit(alpha);
  real corr_factor = normal_rng(1.0875, 0.03625);   // one draw, shared

  {
    vector[N_flat] eta = alpha[observer_idx] + BF_Fixed * bf + VS_Fixed * vs;
    vector[N_flat] log_kappa = to_vector(log_N_latent)[flat_idx]
                             + log_watch_length
                             + log_inv_logit(eta);
    for (i in 1:N_flat)
      log_lik[i] = likelihood_NB
                 ? neg_binomial_2_log_lpmf(n[i] | log_kappa[i], phi[1])
                 : poisson_log_lpmf(n[i] | log_kappa[i]);
  }

  for (y in 1:n_year) {
    int arg_t = 2;
    real best = negative_infinity();
    real half;
    int w = 0;

    for (t in 2:(n_days - 1))
      if (log_N_latent[t, y] > best) { best = log_N_latent[t, y]; arg_t = t; }

    // Parabolic refinement on the log curve, for sub-day resolution.
    if (arg_t > 2 && arg_t < n_days - 1) {
      real fm = log_N_latent[arg_t - 1, y];
      real f0 = log_N_latent[arg_t,     y];
      real fp = log_N_latent[arg_t + 1, y];
      real den = fm - 2 * f0 + fp;
      peak_day[y] = arg_t + (den < -1e-10 ? 0.5 * (fm - fp) / den : 0.0);
    } else {
      peak_day[y] = arg_t;
    }
    peak_height[y] = exp(best);

    half = 0.5 * peak_height[y];
    for (t in 2:(n_days - 1))
      if (exp(log_N_latent[t, y]) > half) w += 1;
    peak_width[y] = w;

    Raw_Est[y] = sum(exp(col(log_N_latent, y)));
    Corrected_Est[y] = Raw_Est[y]
                     * (independent_corr ? normal_rng(1.0875, 0.03625) : corr_factor);
  }

  {
    real ybar = mean(year_values);
    real pbar = mean(peak_day);
    real num = 0;
    real den = 0;
    for (y in 1:n_year) {
      num += (year_values[y] - ybar) * (peak_day[y] - pbar);
      den += square(year_values[y] - ybar);
    }
    peak_day_slope  = num / den;
    peak_day_decade = 10 * peak_day_slope;
  }
}
