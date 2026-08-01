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
//     Set anchor_sd = 1e-6 for parity; set it to the SE of logit(0.80)
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
  real<lower=0> alpha_S1;  real<lower=0> beta_S1;
  real<lower=0> alpha_S2;  real<lower=0> beta_S2;
  
  real<lower=0> mu_S1;      // mean of S1
  real<lower=0> shape_S1;   // shape

  real<lower=0> mu_S2;      // mean of S2
  real<lower=0> shape_S2;   // shape


  vector[n_observer] alpha;                    // cell means, NO global intercept
  real<lower=0, upper=sigma_Obs_max> sigma_Obs;   // dunif(0, sigma_Obs_max)
  real logit_p0;                               // anchor, ~fixed when anchor_sd tiny

  real BF_Fixed;
  real VS_Fixed;

  real<lower=30, upper=60> beta0_P;            // dunif(30, 60)
  real beta1_P;
  real<lower=0> sigma_proc_P;
  vector[n_year] P_raw;

  real beta0_Max;
  real beta1_Max;
  real<lower=0> sigma_proc_Max;
  vector[n_year] log_Max_raw;

  real<lower=0, upper=phi_max> phi;            // dunif(0, phi_max)

  matrix[use_process_error ? n_days : 0,
         use_process_error ? n_year : 0] log_N_raw;
  array[use_process_error] real<lower=0> sigma_process;

  sum_to_zero_vector[n_period] shape_raw;
  real<lower=0> sigma_shape;
}

transformed parameters {
  vector[n_year] P;
  vector[n_year] log_Max;
  matrix[n_days, n_year] log_N_latent;

  for (y in 1:n_year) {
    P[y]       = beta0_P   + beta1_P   * year_values[y] + sigma_proc_P   * P_raw[y];
    log_Max[y] = beta0_Max + beta1_Max * year_values[y] + sigma_proc_Max * log_Max_raw[y];
  }

  {
    vector[n_period] shape_dev = use_shape_dev ? sigma_shape * shape_raw
                                               : rep_vector(0.0, n_period);
    for (y in 1:n_year) {
      real inv_S1 = 1.0 / S1[y];
      real inv_S2 = 1.0 / S2[y];

      // JAGS fixes the first and last day; no observations fall on them.
      log_N_latent[1, y]      = log_boundary_N;
      log_N_latent[n_days, y] = log_boundary_N;

      for (t in 2:(n_days - 1)) {
        real p_minus_t = P[y] - t;
        // C1 = [1 + (2e - 1) exp(-(P-t)/S1)]^(-1/e)
        // C2 = [1 + (2e - 1) exp(+(P-t)/S2)]^(-1/e)
        real log_C1 = -inv_exp_gamma * log1p_exp(log_gamma_term - inv_S1 * p_minus_t);
        real log_C2 = -inv_exp_gamma * log1p_exp(log_gamma_term + inv_S2 * p_minus_t);

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
  beta1_P        ~ normal(0, sd_beta1_P);
  sigma_proc_P   ~ normal(0, sd_sigma_proc_P);      // half-normal via <lower=0>
  P_raw          ~ std_normal();

  beta0_Max      ~ normal(beta0_Max_mu, sd_beta0_Max);
  beta1_Max      ~ normal(0, sd_beta1_Max);
  sigma_proc_Max ~ normal(0, sd_sigma_proc_Max);
  log_Max_raw    ~ std_normal();

  // ---- Shape parameters, hyperparameters estimated ----------------------
  //S1       ~ gamma(alpha_S1, beta_S1);
  S1 ~ gamma(shape_S1, shape_S1 / mu_S1);
  //S2       ~ gamma(alpha_S2, beta_S2);
  S2 ~ gamma(shape_S2, shape_S2 / mu_S2);

  alpha_S1 ~ normal(alpha_S_mu, alpha_S_sd);        // half-normal via <lower=0>
  alpha_S2 ~ normal(alpha_S_mu, alpha_S_sd);
  beta_S1  ~ gamma(beta_S_shape, beta_S_rate);
  beta_S2  ~ gamma(beta_S_shape, beta_S_rate);

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
  shape_raw   ~ std_normal();
  sigma_shape ~ normal(0, sd_sigma_shape);

  // ---- Observation model -------------------------------------------------
  {
    vector[N_flat] eta = alpha[observer_idx] + BF_Fixed * bf + VS_Fixed * vs;
    vector[N_flat] log_kappa = to_vector(log_N_latent)[flat_idx]
                             + log_watch_length
                             + log_inv_logit(eta);
    n ~ neg_binomial_2_log(log_kappa, phi);
  }
}

generated quantities {
  vector[n_year] Raw_Est;
  vector[n_year] Corrected_Est;
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
      log_lik[i] = neg_binomial_2_log_lpmf(n[i] | log_kappa[i], phi);
  }

  for (y in 1:n_year) {
    Raw_Est[y] = sum(exp(col(log_N_latent, y)));
    Corrected_Est[y] = Raw_Est[y]
                     * (independent_corr ? normal_rng(1.0875, 0.03625) : corr_factor);
  }
}
