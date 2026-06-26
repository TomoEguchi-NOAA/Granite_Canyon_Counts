// A STAN model for Hierarchical State Space Modeling using 
// the Richards function to estimate gray whale abundance from 
// count data at Granite Canyon, CA.  

// The current model naming convention is: 
// M1: P = 'time', Max = 'time', S1 = 'time', S2 = 'time, K = 1 
// M2: P = 'time', Max = 'time', S1 = 'S1', S2 = 'S2', K = 1 
// M3: P = 'time', Max = 'time', S1 = 'time', S2 = 'S2', K = 1 
// M4: P = 'time', Max = 'time', S1 = 'S1', S2 = time, K = 1 

// Poisson likelihood (Poisson): a1 
// Negative binomial likelihood (NegBin): a2 

data {
  int<lower=1> n_year;
  int<lower=1> n_days;
  int<lower=1> n_obs;        
  int<lower=1> N_flat;      
  
  array[N_flat] int n;      
  array[N_flat] real bf;    
  array[N_flat] real vs;    
  array[N_flat] int obs;    
  vector[N_flat] watch_length;
  
  // Mapping indices for the observation loop
  array[N_flat] int day_idx;
  array[N_flat] int year_idx;
  
  real<lower=0> S1_alpha; real<lower=0> S1_beta;
  real<lower=0> S2_alpha; real<lower=0> S2_beta;
  real log_K_mu;          real<lower=0> log_K_sigma;
}

parameters {
  vector<lower=0>[n_year] Max;
  //vector<lower=0, upper=20>[n_year] S1;
  vector<lower=0, upper=20>[n_year] S2;
  vector<lower=1e-3>[n_year] K;
  vector<lower=30, upper=60>[n_year] P;
  
  vector[n_obs] OBS_RF;     
  //vector<lower=0>[n_year] r; // Dispersion parameter per year
  
  real mean_prob;
  real BF_Fixed;
  real VS_Fixed;
  real<lower=0,upper=20> S1;
  //real<lower=0,upper=20> S2;
  
  real<lower=0> sigma_Obs;
  real<lower=0> sigma_process;
  
  // Latent abundance process (continuous approximation)
  matrix[n_days, n_year] log_N_raw;
}

transformed parameters {
  vector[N_flat] obs_prob;
  matrix[n_days, n_year] mean_N;
  matrix[n_days, n_year] N_latent; // The "True" Abundance
  matrix[n_days, n_year] log_N_latent;
  
  // 1. Calculate Detection Probability
  for (i in 1:N_flat) {
    obs_prob[i] = inv_logit(mean_prob + (BF_Fixed * bf[i]) + (VS_Fixed * vs[i]) + 
                            OBS_RF[obs[i]]); 
  }

  // 2. Richards Migration Curve
  for (y in 1:n_year) {
    for (t in 1:n_days) {
      if (t == 1 || t == n_days) {
        mean_N[t, y] = 1e-6; 
      } else {
        real base1 = 1 + (2 * exp(K[y]) - 1) * exp((1 / (-S1)) * (P[y] - t));
        real base2 = 1 + (2 * exp(K[y]) - 1) * exp((1 / S2[y]) * (P[y] - t));
        mean_N[t, y] = Max[y] * (pow(max([base1, 1e-9]), -1/exp(K[y])) * pow(max([base2, 1e-9]), -1/exp(K[y])));
      }
      
      // Boundary conditions for Latent State
      if (t == 1 || t == n_days) {
		N_latent[t, y] = 1e-6;
		log_N_latent[t,y] = -11.5;
	  } else {
	    log_N_latent[t, y] = log(mean_N[t, y] + 1e-6) + (sigma_process * log_N_raw[t, y]);
        N_latent[t, y] = exp(log_N_latent[t, y]);
	  }
    }
  }
}

model {
  // Priors
  Max ~ normal(2000, 800);
  S1 ~ gamma(S1_alpha, S1_beta);
  S2 ~ gamma(S2_alpha, S2_beta);
  K ~ lognormal(log_K_mu, log_K_sigma);
  P ~ uniform(30, 60);
  mean_prob ~ uniform(-1.5, 3);
  BF_Fixed ~ normal(0, 30);
  VS_Fixed ~ normal(0, 30);
  sigma_Obs ~ uniform(0, 3);
  OBS_RF ~ normal(0, sigma_Obs);
  sigma_process ~ normal(0, 1);
  //r ~ uniform(0, 50);

  // --- PROCESS MODEL ---
  for (y in 1:n_year) {
    for (t in 2:(n_days - 1)) {
      log_N_latent[t, y] ~ normal(log(mean_N[t, y] + 1e-6), sigma_process);
    }
  }

  // --- OBSERVATION MODEL (Negative Binomial) ---
  for (i in 1:N_flat) {
    int d = day_idx[i];
    int y = year_idx[i];
    
    // N_latent[d,y] is the true abundance; obs_prob[i] is the detection rate
    // We use the same Negative Binomial parameterization as your JAGS model:
    // n ~ NegBin(kappa, r) where kappa = N * p
    // 1. Calculate the number of whales that could be seen given the effort
    real whales_available = N_latent[d, y] * watch_length[i];
    
    // 2. Expected mean (kappa)
    real kappa = whales_available * obs_prob[i];
    
    // 3. Negative Binomial Likelihood
    n[i] ~ poisson(kappa);
  }
}

generated quantities {
  vector[n_year] Raw_Est;
  vector[n_year] Corrected_Est;
  
  // Define correction factor parameters
  real mean_corr = 1.0875;
  real sd_corr = 0.03625;
  
  for (y in 1:n_year) {
    // Sum the latent abundance for the year
    Raw_Est[y] = sum(N_latent[1:n_days, y]);
    
    // Sample the correction factor distribution for this iteration
    real corr_factor = normal_rng(mean_corr, sd_corr);
    
    // Apply correction
    Corrected_Est[y] = Raw_Est[y] * corr_factor;
  }
}
