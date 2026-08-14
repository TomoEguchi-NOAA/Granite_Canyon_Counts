# Comparing all models in STAN
# 
#  Collects Stan runs and compares results. Do posterior predictive checks on
#  all models. 
#  
# To summarize sensitivity analyses, use Summary_Stan_Sensitivity.R
# 
# This script was converted from Stan_Richards_AllData.R
# 

rm(list = ls())
library(tidyverse)
library(posterior)
library(bayesplot)

source("Granite_Canyon_Counts_fcns.R")
source("ppc_richards_hssm.R")

min.dur <- 60 #10 #85 #
Run.date <- Sys.Date()

# These are the ending year of each season - for example, 2022 in the following vector indicates
# for the 2021/2022 season. These data were extracted using Extract_Data_All_v2.Rmd
# Data prior to the 2009/2010 season are in Laake's ERAnalayis package. 
years <- c(2008, 2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024, 2025, 2026)
data.dir <- "RData//V2.1_May2026"
max.day <- 100

jags.input <- NoBUGS_Jags_input(min.dur = min.dur, 
                                years = years, 
                                data.dir = data.dir, 
                                max.day = max.day, 
                                obs.n.min = 10, N.obs = 10)

jags.data <- jags.input$jags.data
start.years <- jags.data$start.years

#stan.data <- create.stan.data(jags.data = jags.data)

# // ---- MODEL STRUCTURE (Table 1) ---------------------------------------
#   //   S1_by_season  S2_by_season  likelihood_NB      model
# //        1             1              0/1          M1a1 / M1a2
# //        0             0              0/1          M2a1 / M2a2
# //        1             0              0/1          M3a1 / M3a2
# //        0             1              0/1          M4a1 / M4a2

# ZI0 and ZI1 indicate without or with constant probability zero-inflation model

# Model comparison is done with no gamma parameter. Adding one gamma is used as
# a refinement to the best model. Zero inflated models also are use as another
# refinement given the excessive zeros found in the best model + a gamma parameter
models <- c("M1a1_0gamma_mod5_ZI0", "M2a1_0gamma_mod5_ZI0", 
            "M3a1_0gamma_mod5_ZI0", "M4a1_0gamma_mod5_ZI0",
            "M1a2_0gamma_mod5_ZI0", "M2a2_0gamma_mod5_ZI0", 
            "M3a2_0gamma_mod5_ZI0", "M4a2_0gamma_mod5_ZI0",
            "M1a1_0gamma_mod5_ZI1", "M2a1_0gamma_mod5_ZI1", 
            "M3a1_0gamma_mod5_ZI1", "M4a1_0gamma_mod5_ZI1",
            "M1a2_0gamma_mod5_ZI1", "M2a2_0gamma_mod5_ZI1", 
            "M3a2_0gamma_mod5_ZI1", "M4a2_0gamma_mod5_ZI1")
#            "M1a2_1gamma_mod5_ZI1", "M1a2_1gamma_mod5_ZI2")

# Parameters to monitor for Poisson models:
params.a1.stan <- c("S1", "S2", "P", 
                    "sigma_proc_P", "Corrected_Est", 
                    "Max", "log_N_latent",
                    "peak_day_decade", "beta1_P")

# Negative Binomial model has the over dispersion parameter
params.a2.stan <- c(params.a1.stan, "phi")

diag.summary <- LOO.out <- stan.global.summary <- ppc.res <- list()
#zi_mode <- 1
k <- 3

if (!file.exists("RData//Stan_Richards_Results.rds")){
  for (k in 1:length(models)) {
    out.file <- paste0("Richards_HSSM_", models[k], "_stan")
    
    #mod <- cmdstan_model(file)
    fit_stan <- readRDS(paste0("RData//", out.file, ".rds"))
    
    if (length(grep("a1", models[k])) > 0) params <- params.a1.stan
    if (length(grep("a2", models[k])) > 0) params <- params.a2.stan
    
    # If zero-inflated, may have one or two extra parameters:
    if (length(grep("ZI1", models[k])) > 0) params <- c(params, "zi_a")
    if (length(grep("ZI2", models[k])) > 0) params <- c(params, "zi_a", "zi_b")
    
    # The default zi_mode = 1. So, mod4 runs have to set zi_mode = 0
    #if (length(grep("ZI0", models[k])) == 0) zi_mode = 0
    
    stan.info <- readRDS(paste0("RData//", out.file, "_info.rds"))
    stan.data <- stan.info$stan.data
    #stan.data$zi_mode <- zi_mode
    
    LOO.out[[k]] <- fit_stan$loo()
    stan.global.summary[[k]] <- fit_stan$summary(
      variables = params,
      default_summary_measures(), 
      default_convergence_measures(),
      extra_quantiles = ~quantile2(., probs = c(0.025, 0.975)))
    
    diag.summary[[k]] <- fit_stan$diagnostic_summary()
    
    # --- Posterior Predictive Simulation Loop ---
    # res$autocorr is the key one. Within-season lag-1 autocorrelation of daily-averaged 
    # residuals, compared against the same statistic computed on replicate datasets. 
    # If the Richards curve were too rigid to track the migration, residuals would 
    # come in same-signed runs along the season and the observed autocorrelation 
    # would exceed the replicate distribution. If the observed value sits inside 
    # the replicate interval, it is a direct evidence the curve is flexible enough.
    
    # plots$resid_day is the same thing visually. A flat loess through residuals 
    # against day-of-season says the curve captures the shape; systematic 
    # curvature — a dip at the peak, say — would say it doesn't.
    
    # pvalues covers proportion of zeros, maximum, SD, mean, and the upper tail. 
    # The flag column marks anything outside [0.025, 0.975].
    
    # plots$coverage gives the empirical coverage of 95% predictive intervals. 
    # This is the cleanest single number for the Results — if it lands near 95%, 
    # the observation model is calibrated. It also detects Poisson vs NB 
    # automatically by looking for phi[1], so run it unchanged across all 
    # eight models. Worth doing: comparing coverage and the zeros p-value between 
    # M1a1 and M1a2 is a much more direct justification for the negative binomial 
    # than ΔLOOIC is.
    
    # plots$resid_season faceted QQ plots will show whether any particular season 
    # fits badly. Worth checking against the ones with sparse effort, and 
    # especially 2025/2026.
    
    # I increase the n_draws to 1500. Results seem to change if not enough samples
    # are drawn. 1500 seems to be in a good place.
    ppc.res[[k]] <- ppc_richards(fit_stan, stan.data, n_draws = 1500)
    
    # res$pvalues
    # res$autocorr
    # res$plots$resid_day
  }
  
  models.out <- list(ppc = ppc.res,
                     LOO = LOO.out,
                     diag = diag.summary,
                     summary = stan.global.summary)
  
  saveRDS(models.out, file = "RData//Stan_Richards_Results.rds")
  
} else {
  models.out <- readRDS(file = "RData//Stan_Richards_Results.rds")
}

diag.summary <- models.out$diag
LOO.out <- models.out$LOO
ppc.res <- models.out$ppc
stan.global.summary <- models.out$summary

rm(list = "models.out")

### find the best model:

diag.n_divergent.df <- do.call(rbind,
                               lapply(diag.summary, 
                                      FUN = function(x) x$num_divergent)) %>%
  as.data.frame()

# if any of the entries in the above dataframe is >>0, uncomment the following chunk:
# diag.n_divergent.df$Model <- as.vector(models)
# 
# diag.n_max_treedepth.df <- do.call(rbind,
#                                lapply(diag.summary, FUN = function(x) x$num_max_treedepth)) %>%
#   as.data.frame()
# diag.n_max_treedepth.df$Model <- as.vector(models)
# 
# diag.ebfmi.df <- do.call(rbind,
#                          lapply(diag.summary, FUN = function(x) x$ebfmi)) %>%
#   as.data.frame()
# diag.ebfmi.df$Model <- as.vector(models)
##################################################################################

# --- Model comparison ---
LOOIC.df <- do.call(rbind, lapply(LOO.out, FUN = function(x) x$estimates["looic",])) %>%
  as.data.frame()
LOOIC.df$Model <- as.vector(models) 
LOOIC.df %>%  select(Model, Estimate, SE) %>%
  mutate(dLOOIC = Estimate- min(Estimate),
         rownames_to_column(.,var = "ID")) %>% #-> tmp
  arrange(by = dLOOIC) -> LOOIC.df

if (!file.exists("RData//LOOIC.rds"))
  saveRDS(LOOIC.df, file = "RData//LOOIC.rds")

Top.4 <- LOOIC.df[1:4, "ID"] %>% 
  as.numeric()

# Excluding NB models:
loo::loo_compare(LOO.out[Top.4]) %>%
  data.frame() %>%
  mutate(Model = LOOIC.df[1:4, "Model"])-> loo.comp.df

if (!file.exists("RData//LOO_compare.rds"))
  saveRDS(loo.comp.df, "RData//LOO_compare.rds")

best.model <- Top.4[1]    # running all models but excluding NB models
pk <- LOO.out[[best.model]]$diagnostics$pareto_k

best.model.name <- models[best.model]

# Bring in the data for the best fit model - some switches in the data list 
# are different among models:
# fit.info <- readRDS(paste0("RData//Richards_HSSM_", best.model.name, "_stan_info.rds"))
# stan.data <- fit.info$stan.data

which(pk > 0.7)
if (length(which(pk>0.7)) > 0){
  sd <- stan.data
  data.frame(day = sd$day_idx, 
             year = sd$year_idx, 
             n = sd$n)[which(pk > 0.7), ]
  
}

