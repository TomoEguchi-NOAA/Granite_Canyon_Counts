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

# To run all models, keep the following three lines:
S1.vec <- c(1, 0)
S2.vec <- c(1, 0)
lkhd.NB.vec <- c(1,0)

# To run a specific model, set the vectors accordingly:
S1.vec <- c(1)
S2.vec <- c(1)
lkhd.NB.vec <- c(1)

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
file <- file.path("models/model_Richards_HSSM_mod3.stan")
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
      out.file <- paste0("Richards_HSSM_", model, "_1gamma_stan")
      #m <- m + 1
      # Compile with aggressive C++ optimization flags
      mod <- cmdstan_model(file, 
                           cpp_options = list(stan_threads = TRUE, 
                                              O = 3))
      
      stan.data <- create.stan.data(jags.data = jags.data,
                                    S1_by_season = S1_by_season,
                                    S2_by_season = S2_by_season,
                                    likelihood_NB = likelihood_NB,
                                    estimate_gamma = 1,
                                    separate_gamma = 0)
      
      n_year <- stan.data$jags.data$n.year
      n_observer <- stan.data$jags.data$n.obs.fixed
      
      #mod <- cmdstan_model(file)
      if (!file.exists(paste0("RData/", out.file, ".rds"))){
        fit_stan <- mod$sample(
          data            = stan.data$stan.data,
          init            = stan_init_fn,
          chains          = 4,
          parallel_chains = 4,
          threads_per_chain = 2,
          iter_warmup     = 1500,
          iter_sampling   = 2000,
          adapt_delta     = 0.90  
        )
        
        fit_stan$save_object(file = paste0("RData/", out.file, ".rds"))
        saveRDS(list(stan.data = stan.data$stan.data,
                     jags.data = stan.data$jags.data,
                     init_fn = stan_init_fn()),
                file = paste0("RData/", out.file, "_info.rds"))
        
      } else {
        fit_stan <- readRDS(paste0("RData/", out.file, ".rds"))
      }
    }
  }
}


