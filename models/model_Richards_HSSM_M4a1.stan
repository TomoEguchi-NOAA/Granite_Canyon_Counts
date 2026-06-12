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
  array[N_flat] real bf;     // FIXED: Now correctly declared as an array of reals
  array[N_flat] real vs;     
  array[N_flat] int obs;     
  vector[N_flat] watch_length;
  
  array[n_days, n_year] int start_idx;
  array[n_days, n_year] int end_idx;
  
  real<lower=0> S1_alpha; real<lower=0> S1_beta;
  real<lower=0> S2_alpha; real<lower=0> S2_beta;
  real log_K_mu;          real<lower=0> log_K_sigma;
}

parameters {
  vector<lower=0>[n_year] Max;
  vector<lower=0, upper=20>[n_year] S2;
  vector<lower=1e-3>[n_year] K;
  vector<lower=30, upper=60>[n_year] P;
  
  real S1; 
  real mean_prob;
  real BF_Fixed;
  real VS_Fixed;
  real<lower=0> sigma_Obs;
  vector[n_obs] OBS_RF;     
  
  matrix<lower=1e-3, upper=50>[n_days, n_year] r;
}

transformed parameters {
  vector[N_flat] obs_prob;
  matrix[n_days, n_year] mean_N;

  for (i in 1:N_flat) {
    obs_prob[i] = inv_logit(mean_prob + 
                            (BF_Fixed * bf[i]) + 
                            (VS_Fixed * vs[i]) + 
                            OBS_RF[obs[i]] + 
                            log(watch_length[i] + 1e-9)); 
  }

  for (y in 1:n_year) {
    for (t in 1:n_days) {
      if (t == 1 || t == n_days) {
        mean_N[t, y] = 1e-6; 
      } else {
        real base1 = 1 + (2 * exp(K[y]) - 1) * exp((1 / (-S1)) * (P[y] - t));
        real base2 = 1 + (2 * exp(K[y]) - 1) * exp((1 / S2[y]) * (P[y] - t));
        
        real safe_base1 = max([base1, 1e-9]);
        real safe_base2 = max([base2, 1e-9]);
        
        real M1 = pow(safe_base1, -1 / exp(K[y]));
        real M2 = pow(safe_base2, -1 / exp(K[y]));
        
        mean_N[t, y] = Max[y] * (M1 * M2);
      }
    }
  }
}

model {
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

// --- Optimized Likelihood (Vectorized Marginalization Loop) ---
  for (y in 1:n_year) {
    for (t in 1:n_days) {
      int start = start_idx[t, y];
      int end = end_idx[t, y];
      
      if (start <= end) { // Check if we have observations for this specific day
        
        // Find maximum observed count for the day
        int max_n = 0;
        for (i in start:end) {
          if (n[i] > max_n) max_n = n[i];
        }
        
        // --- OPTIMIZATION: Tighter integration buffer ---
        // Changed from +120 to +50. Adjust this if your detection prob is tiny.
        int K_val = max_n + 50; 
        vector[K_val - max_n + 1] lp;
        
        if (t == 1 || t == n_days) {
          real log_prob_N = poisson_lpmf(0 | mean_N[t, y]);
          
          // --- OPTIMIZATION: Vectorized Binomial ---
          real log_prob_obs = binomial_lpmf(n[start:end] | 0, obs_prob[start:end]);
          
          target += log_prob_N + log_prob_obs;
        } else {
          for (N_val in max_n:K_val) {
            real log_prob_N = poisson_lpmf(N_val | mean_N[t, y]);
            
            // --- OPTIMIZATION: Vectorized Binomial ---
            // Eliminates the inner 'for (i in start:end)' loop entirely!
            real log_prob_obs = binomial_lpmf(n[start:end] | N_val, obs_prob[start:end]);
            
            lp[N_val - max_n + 1] = log_prob_N + log_prob_obs;
          }
          target += log_sum_exp(lp); 
        }
      }
    }
  }
  
}
