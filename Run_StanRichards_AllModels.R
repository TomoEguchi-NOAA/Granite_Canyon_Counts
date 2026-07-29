# Running the best model in STAN
# 
# Translation of the JAGS Richards model to Stan is possible because 
# if N ~ Poisson(A) and n | N ~ Binomial(N, p), then n ~ Poisson(Ap) exactly;
# 
# The manuscript's whole argument for single-observer estimation rests on anchoring baseline detection at 0.80 via an informative prior — but this code doesn't do that. With N and p both free and only single-observer counts, the ridge is essentially flat, and the sampler will wander along it forever. This alone could account for the failure independently of the clamp. Corrected to a normal prior centred at logit(0.80) = 1.3863, with the SD passed as data - can do sensitivity analysis. SD of the detection probability was set at 0.1622 as the default.
# 
# 
rm(list = ls())
library(tidyverse)
library(posterior)

# Need to install cmdnstanr from the Stan package repo:
# # Tell R to look at the Stan repository instead of CRAN
#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)

# This downloads and builds the C++ toolchain Stan needs to run fast
#check_cmdstan_toolchain(fix = TRUE) # Verifies you have a C++ compiler
#install_cmdstan(cores = 4)          # Downloads and installs CmdStan

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


jags.input <- NoBUGS_Jags_input(min.dur = min.dur, 
                                years = years, 
                                data.dir = data.dir, 
                                max.day = max.day, 
                                obs.n.min = 10, N.obs = 10)

jags.data <- jags.input$jags.data

# // ---- MODEL STRUCTURE (Table 1) ---------------------------------------
#   //   S1_by_season  S2_by_season  likelihood_NB      model
# //        1             1              0/1          M1a1 / M1a2
# //        0             0              0/1          M2a1 / M2a2
# //        1             0              0/1          M3a1 / M3a2
# //        0             1              0/1          M4a1 / M4a2

S1.vec <- c(1, 0)
S2.vec <- c(1, 0)
lkhd.NB.vec <- c(1,0)

# --- Add these lines right before packaging 'stan_data' ---
# storage.mode(start_idx) <- "integer"
# storage.mode(end_idx)   <- "integer"

# The following plot can be useful: Set it aside for the manuscript. Predicted-versus-observed counts across all watch periods and all 34 seasons, correlation 0.88, with the scatter fanning out proportionally to the mean exactly as a negative binomial should — that is a direct, visual answer to reviewer #2. It shows the Richards curve tracking the observed counts across the full range without any spline flexibility, and the mean-variance relationship supporting the NB choice.

# It's also most of the posterior predictive check your Methods promise and never deliver. Regenerate it from posterior draws rather than plug-in medians, add the 95% predictive band, and colour by season — that plus the within-season residual autocorrelation check gives you the whole GOF section.
# 
# reconstruct kappa at the JAGS posterior medians, using the STAN data list
# Pm <- apply(jm.out$jm$sims.list$P, 2, median)
# Mm <- apply(jm.out$jm$sims.list$Max, 2, median)
# S1m <- apply(jm.out$jm$sims.list$S1, 2, median)
# S2m <- apply(jm.out$jm$sims.list$S2, 2, median)
# g <- exp(1)
# pmt <- Pm[stan_data$year_idx] - stan_data$day_idx
# C1 <- (1 + (2*g-1)*exp(-pmt/S1m[stan_data$year_idx]))^(-1/g)
# C2 <- (1 + (2*g-1)*exp( pmt/S2m[stan_data$year_idx]))^(-1/g)
# pred <- Mm[stan_data$year_idx]*C1*C2 * stan_data$watch_length * 0.8
# plot(pred, stan_data$n); abline(0,1,col="red")
# cor(pred, stan_data$n)

# Create an inits function

# all models were fit into one file with various switches:
file <- file.path("models/model_Richards_HSSM_mod2.stan")
#model <- list()
#m <- 1

for (S1 in 1:length(S1.vec)){
  S1_by_season <- ifelse(S1.vec[S1] == 1, 1, 0)
  for (S2 in 1:length(S2.vec)){
    S2_by_season <- ifelse(S2.vec[S2] == 1, 1, 0)
    for (lkhd in 1:length(lkhd.NB.vec)){
      likelihood_NB <- ifelse(lkhd.NB.vec[lkhd] == 1, 1, 0)
      model.lkhd <- ifelse(lkhd.NB.vec[lkhd] == 1, "a2", "a1")      
      if (S1_by_season == 1 & S2_by_season == 1){
        model.M <- "M1"
      } else if (S1_by_season == 1 & S2_by_season == 0){
        model.M <- "M3"
      } else if (S1_by_season == 0 & S2_by_season == 1){
        model.M <- "M4"
      } else {
        model.M <- "M2"
      }
      
      model <- paste0(model.M, model.lkhd)
      out.file <- paste0("Richards_HSSM_", model, "_stan")
      #m <- m + 1
      # Compile with aggressive C++ optimization flags
      mod <- cmdstan_model(file, 
                           cpp_options = list(stan_threads = TRUE, 
                                              O = 3))
      
      stan.data <- create.stan.data(jags.data = jags.data,
                                    S1_by_season = S1_by_season,
                                    S2_by_season = S2_by_season,
                                    likelihood_NB = likelihood_NB)
      n_year <- stan.data$jags.data$n.year
      n_observer <- stan.data$jags.data$n.obs.fixed
      
      init_fn <- function() {
        out <- list(
          beta0_Max = rnorm(1, 7.6, 0.2),  
          #beta1_Max = rnorm(1, 0, 0.05),
          #log_Max_raw = rnorm(n_year, 0, 0.1),  
          sigma_proc_Max = runif(1, 0.2, 0.6),
          beta0_P = runif(1, 42, 48),      
          #beta1_P = rnorm(1, 0.21, 0.03),
          P_raw = rnorm(n_year, 45, 2),      # was rnorm(n_year, 0, 0.1)
          log_Max_raw = rnorm(n_year, 7.6, 0.2),  # was rnorm(n_year, 0, 0.1)
          #P_raw = rnorm(n_year, 0, 0.1),   
          #sigma_proc_P = runif(1, 3.5, 5.5),
          S1 = runif(n_year, 2, 4),        
          S2 = runif(n_year, 2, 4),
          mu_S1 = runif(1, 2.5, 3.5),  
          shape_S1 = runif(1, 8, 12),
          mu_S2 = runif(1, 2.5, 3.5),  
          shape_S2 = runif(1, 8, 12),
          #alpha_S1 = runif(1, 8, 12),      beta_S1 = runif(1, 2.5, 4.5),
          #alpha_S2 = runif(1, 8, 12),      beta_S2 = runif(1, 2.5, 4.5),
          alpha = rnorm(n_observer, 1.39, 0.1),  
          sigma_Obs = runif(1, 0.15, 0.35),
          logit_p0 = rnorm(1, 1.3863, 0.01),
          BF_Fixed = rnorm(1, 0, 0.1),     
          VS_Fixed = rnorm(1, 0, 0.1),
          sigma_shape = runif(1, 0.05, 0.2)
        )
        # only supply inits for parameters that actually exist in this configuration
        
        if (stan.data$stan.data$use_trend_P)   out$beta1_P   <- array(rnorm(1, 0.21, 0.03), dim = 1)
        if (stan.data$stan.data$use_trend_Max) out$beta1_Max <- array(rnorm(1, 0,    0.05), dim = 1)
        if (stan.data$stan.data$use_plateau) {
          k <- if (stan.data$stan.data$plateau_by_year) n_year else 1
          out$delta <- array(runif(k, 0.5, 2), dim = k)
        }
        if (stan.data$stan.data$use_process_error) {
          out$log_N_raw     <- matrix(rnorm(n_days * n_year, 0, 0.1), n_days, n_year)
          out$sigma_process <- array(runif(1, 0.1, 0.3), dim = 1)
        }
        
        if (stan.data$stan.data$S1_by_season) {
          out$mu_S1 <- array(runif(1,2.5,3.5),1)
          out$shape_S1 <- array(runif(1,8,12),1)
        }
        
        if (stan.data$stan.data$S2_by_season) {
          out$mu_S2 <- array(runif(1,2.5,3.5),1)
          out$shape_S2 <- array(runif(1,8,12),1)
        }
        
        if (stan.data$stan.data$likelihood_NB == 1) out$phi = runif(1, 4, 6)
        return(out)
      }
      
      #mod <- cmdstan_model(file)
      if (!file.exists(paste0("RData/", out.file, ".rds"))){
        fit_stan <- mod$sample(
          data            = stan.data$stan.data,
          init            = init_fn,
          chains          = 4,
          parallel_chains = 4,
          threads_per_chain = 2,
          iter_warmup     = 1000,
          iter_sampling   = 1000,
          adapt_delta     = 0.90  
        )
        
        fit_stan$save_object(file = paste0("RData/", out.file, ".rds"))
        saveRDS(list(stan.data = stan.data$stan.data,
                     jags.data = stan.data$jags.data,
                     init_fn = init_fn()),
                file = paste0("RData/", out.file, "info.rds"))
        
      } else {
        fit_stan <- readRDS(paste0("RData/", out.file, ".rds"))
      }
    }
  }
}


