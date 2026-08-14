# Sensitivity analysis of Stan models as per Claude's suggestions
# with some modification to them

# Stan_Richards_Model_Comparison.R
# LOOIC (elpd) indicated that M1a2 with zero inflation was the 
# best model, with M1a2 without zero inflation in the second
# dLOOIC = 29.2. These two models are modified to see if they
# can be better. Here is the list of tweaking:
# 
#  1. estimate_gamma = 1 - provide a distribution to the gamma parameter
#  2. Two gamma parameters 
#  3. Second parameter for ZI part (intercept + slope * (log_kappa - lk_bar))
#  4. use_trend_P = 0 - is there a trend in P? (sen4)
#  5. use_pooling_Max = 0 — removes hierarchical shrinkage on the season levels
#                           altogether: each log_Max[y] gets an independent diffuse
#                           prior. Use this to test whether shrinkage toward the 
#                           trend inflates low seasons. When 0, beta0_Max / 
#                           beta1_Max / sigma_proc_Max are not identified and 
#                           simply sample their priors (sen1)  
#  6. use_shape_dev = 1 Periodic deviation from the curve is shared
#  7. use_shape_dev = 1 AND n_period = 10. Periodic grouping is changed from 5 to 10 days
#  8. gamma_prior_mu = 0.0
#  9. gamma_prior_sd = 2.0
# 10. anchor_mu = qlogis(0.70) (sen2)
# 11. anchor_mu = qlogis(0.90) — detection sensitivity (sen3)



rm(list = ls())
library(tidyverse)
library(posterior)
library(cmdstanr)

source("Granite_Canyon_Counts_fcns.R")

# Sensitivity analyses are done only on two models (M1a2_0gamma_mod5_ZI1 and M1a2_0gamma_mod5_ZI0)
# The default values for S1, S2, and likelihood in create.stan.data are set
# to run the M1a2 model as of 2026-07-29


# // ---- MODEL STRUCTURE (Table 1) ---------------------------------------
#   //   S1_by_season  S2_by_season  likelihood_NB      model
# //        1             1              0/1          M1a1 / M1a2
# //        0             0              0/1          M2a1 / M2a2
# //        1             0              0/1          M3a1 / M3a2
# //        0             1              0/1          M4a1 / M4a2

run.sensitivity.stan <- function(sd, mod, out.file){
  
  if (!file.exists(paste0("RData//", out.file, ".rds"))){
    fit_stan <- mod$sample(
      data            = sd,
      init            = stan_init_fn,
      chains          = 4,
      parallel_chains = 4,
      threads_per_chain = 2,
      iter_warmup     = 1500,
      iter_sampling   = 2000,
      adapt_delta     = 0.90  
    )
    
    fit_stan$save_object(file = paste0("RData//", out.file, ".rds"))
    saveRDS(list(model.file = file,
                 stan.data = sd,
                 init_fn = stan_init_fn(),
                 System = Sys.getenv(),
                 Run.date = Sys.Date()),
            file = paste0("RData//", out.file, "_info.rds"))
    
  }
}


extract.results <- function(models, params){
  # Bring in all the results, compute LOOOIC
  LOO.out <- diag.summary <- stan.summary <- list()
  for (k in 1:length(models)){
    fit_stan <- readRDS(
      paste0("RData//Richards_HSSM_",
             models[k], "_stan.rds"))
    
    LOO.out[[k]] <- fit_stan$loo()
    if (length(grep("P0", models[k]))){
      params <- params[params != "beta1_P"]
    }
    
    if (length(grep("ZI1", models[k])) > 0) {
      params.1 <- c(params, "zi_a")
    } else if (length(grep("ZI2", models[k])) > 0){
      params.1 <- c(params, "zi_a", "zi_b")
    } else if (length(grep("ZI0", models[k])) > 0){
      params.1 <- params
    }
    
    stan.summary[[k]] <- fit_stan$summary(
      variables = params.1,
      default_summary_measures(), 
      default_convergence_measures(),
      extra_quantiles = ~quantile2(., probs = c(0.025, 0.975)))
    
    diag.summary[[k]] <- fit_stan$diagnostic_summary()
    
  }
  
  LOOIC.df <- do.call(
    rbind, 
    lapply(LOO.out, 
           FUN = function(x) x$estimates["looic",])) %>%
    as.data.frame() %>%
    mutate(Model = as.vector(models)) %>%  
    select(Model, Estimate, SE) %>%
    mutate(dLOOIC = Estimate- min(Estimate),
           rownames_to_column(.,var = "ID")) %>% #-> tmp
    arrange(by = dLOOIC)
  
  out <- list(LOO.out = LOO.out,
              stan.summary = stan.summary,
              diag.summary = diag.summary,
              LOOIC.df = LOOIC.df)
  return(out)
}

mod.file <- file.path("models//model_Richards_HSSM_mod5_ZI.stan")
params <- c("S1", "S2", "P", "sigma_shape",
            "sigma_proc_P", "Corrected_Est", 
            "Max", "log_N_latent",
            "peak_day_decade", "beta1_P", "phi")

## NEED TO DEFINE THESE FIRST. THEY ARE USED IN THE INTITS FUNCTION 2026-08-13
## sd IS STAN DATA. GENERIC STAN INPUT DATA ARE NEEDED
n_year <- sd$n_year
n_observer <- sd$n_observer

# First two models from the initial run:
models <- c("ZI1",
            "ZI0")
mod <- cmdstan_model(mod.file, 
                     cpp_options = list(stan_threads = TRUE, 
                                        O = 3))

# Sensitivity 1: Add a gamma parameter
for (k in 1:length(models)){
  info <- readRDS(paste0("RData//Richards_HSSM_M1a2_0gamma_mod5_",
                         models[k], "_stan_info.rds"))
  
  out.file <- paste0("Richards_HSSM_M1a2_1gamma_mod5_",
                     models[k], "_stan")
  
  sd <- info$stan.data
  sd$estimate_gamma <- 1
  run.sensitivity.stan(sd, mod, out.file)

}

### Compare results among the four
sen.1.models <- c("M1a2_0gamma_mod5_ZI1",
                  "M1a2_0gamma_mod5_ZI0",
                  "M1a2_1gamma_mod5_ZI1",
                  "M1a2_1gamma_mod5_ZI0")

sen.1.out <- extract.results(sen.1.models, params)

saveRDS(sen.1.out, 
        file = "RData//sensitivity_1_out.rds")

# This model selection process ended up with 1gamma models
# are better than 0gamma models. ZI0 is LOOIC(ZI1) + 38.6
# The best modesl are M1a2_1gamma_mod5_ZI1 and M1a2_1gamma_mod5_ZI0

# Sensitivity 2: Add a second gamma parameter

### Add the second gamma parameter to the top two models
models <- c("ZI1",
            "ZI0")
mod <- cmdstan_model(mod.file, 
                     cpp_options = list(stan_threads = TRUE, 
                                        O = 3))
for (k in 1:length(models)){
  info <- readRDS(paste0("RData//Richards_HSSM_M1a2_1gamma_mod5_",
                         models[k], "_stan_info.rds"))
  
  out.file <- paste0("Richards_HSSM_M1a2_2gamma_mod5_",
                     models[k], "_stan")
  
  sd <- info$stan.data
  sd$estimate_gamma <- 1
  sd$separate_gamma <- 1
  run.sensitivity.stan(sd, mod, out.file)
  
}


### Compare results among the four
sen.2.models <- c("M1a2_1gamma_mod5_ZI1",
                  "M1a2_1gamma_mod5_ZI0",
                  "M1a2_2gamma_mod5_ZI1",
                  "M1a2_2gamma_mod5_ZI0")

sen.2.out <- extract.results(sen.2.models, params)
saveRDS(sen.2.out, 
        file = "RData//sensitivity_2_out.rds")

## This selection process ended up with 1gamma. 2gammas didn't
## do well. 
## 2gamma_ZI0 had divergence issues and low ESS
## 2gamma_ZI1 was the best but divergent issues. Also low ESS.
## So, keep 1gamma_ZI1 as the best for now. 
## ZI1 is a good idea.

# Sensitivity 3: Add the second term to Zero inflation:
models <- c("1gamma")
mod <- cmdstan_model(mod.file, 
                     cpp_options = list(stan_threads = TRUE, 
                                        O = 3))
for (k in 1:length(models)){
  info <- readRDS(paste0("RData//Richards_HSSM_M1a2_",
                         models[k], "_mod5_ZI1_stan_info.rds"))
  
  out.file <- paste0("Richards_HSSM_M1a2_", models[k], 
                     "_mod5_ZI2_stan")
  
  sd <- info$stan.data
  sd$estimate_gamma <- 1
  #sd$separate_gamma <- 0
  sd$zi_mode <- 2
  
  run.sensitivity.stan(sd, mod, out.file)
  
}


### Compare results 
sen.3.models <- c("M1a2_1gamma_mod5_ZI1",
                  "M1a2_1gamma_mod5_ZI2")

sen.3.out <- extract.results(sen.3.models, params)
saveRDS(sen.3.out, 
        file = "RData//sensitivity_3_out.rds")

# Look at the 2 ZI parameters:
# stan.summary[[2]] %>%
#   filter(variable %in% c("zi_a[1]", "zi_b[1]"))

# They converged well and estimates seem pretty good. So, 1gamma_ZI2 is best now
# The best model is M1a2_1gamma_mod5_ZI2. dLOOIC = 132.9

## Sensitivity 4: use_trend_P = 0 - is there a trend in P? (sen4)
info <- readRDS("RData//Richards_HSSM_M1a2_1gamma_mod5_ZI2_stan_info.rds")

out.file <- paste0("Richards_HSSM_M1a2_1gamma_mod5_ZI2_P0_stan")

mod <- cmdstan_model(mod.file, 
                     cpp_options = list(stan_threads = TRUE, 
                                        O = 3))

sd <- info$stan.data
sd$estimate_gamma <- 1
sd$zi_mode <- 2
sd$use_trend_P <- 0

run.sensitivity.stan(sd, mod, out.file)

### Compare results 
sen.4.models <- c("M1a2_1gamma_mod5_ZI2",
                  "M1a2_1gamma_mod5_ZI2_P0")

sen.4.out <- extract.results(sen.4.models, params)
saveRDS(sen.4.out, 
        file = "RData//sensitivity_4_out.rds")

## Sensitivity 5: use_pooling_Max = 0
info <- readRDS("RData//Richards_HSSM_M1a2_1gamma_mod5_ZI2_stan_info.rds")

out.file <- paste0("Richards_HSSM_M1a2_1gamma_mod5_ZI2_PoolMax0_stan")

mod <- cmdstan_model(mod.file, 
                     cpp_options = list(stan_threads = TRUE, 
                                        O = 3))

sd <- info$stan.data
sd$estimate_gamma <- 1
sd$zi_mode <- 2
sd$use_trend_P <- 1
sd$use_pooling_Max <- 0

run.sensitivity.stan(sd, mod, out.file)

### Compare results 
sen.5.models <- c("M1a2_1gamma_mod5_ZI2",
                  "M1a2_1gamma_mod5_ZI2_PoolMax0")

sen.5.out <- extract.results(sen.5.models, params)
saveRDS(sen.5.out, 
        file = "RData//sensitivity_5_out.rds")

## Sensitivity 6: use_shape_dev = 1
info <- readRDS("RData//Richards_HSSM_M1a2_1gamma_mod5_ZI2_stan_info.rds")
out.file <- paste0("Richards_HSSM_M1a2_1gamma_mod5_ZI2_ShapeDev1_stan")

mod <- cmdstan_model(mod.file, 
                     cpp_options = list(stan_threads = TRUE, 
                                        O = 3))

sd <- info$stan.data
sd$estimate_gamma <- 1
sd$zi_mode <- 2
sd$use_shape_dev <- 1

run.sensitivity.stan(sd, mod, out.file)

### Compare results 
sen.6.models <- c("M1a2_1gamma_mod5_ZI2",
                  "M1a2_1gamma_mod5_ZI2_ShapeDev1")

sen.6.out <- extract.results(sen.6.models, params)
saveRDS(sen.6.out, 
        file = "RData//sensitivity_6_out.rds")

### THE FOLLOWING IS NOT WORKING 2026-08-13
### n_years is not found in the inits function. It worked above... 

## Sensitiviy 7: use_shape_dev = 1 AND n_period = 10. 
## Periodic grouping is changed from 5 to 10 days
info <- readRDS("RData//Richards_HSSM_M1a2_1gamma_mod5_ZI2_stan_info.rds")
out.file <- paste0("Richards_HSSM_M1a2_1gamma_mod5_ZI2_ShapeDev1_Period10_stan")

mod <- cmdstan_model(mod.file, 
                     cpp_options = list(stan_threads = TRUE, 
                                        O = 3))

sd <- info$stan.data
sd$estimate_gamma <- 1
sd$zi_mode <- 2
sd$use_shape_dev <- 1
sd$n_period <- 10
n_year <- sd$n_year
n_days <- sd$n_days
n_observer <- sd$n_observer
run.sensitivity.stan(sd, mod, out.file)

### Compare results 
sen.7.models <- c("M1a2_1gamma_mod5_ZI2",
                  "M1a2_1gamma_mod5_ZI2_ShapeDev1_Period10")

sen.7.out <- extract.results(sen.7.models, params)
saveRDS(sen.7.out, 
        file = "RData//sensitivity_7_out.rds")

## Sensitivity 8: gamma_prior_mu = 0

info <- readRDS("RData//Richards_HSSM_M1a2_1gamma_mod5_ZI2_stan_info.rds")
out.file <- paste0("Richards_HSSM_M1a2_1gamma_mod5_ZI2_gammaMu0_stan")

mod <- cmdstan_model(mod.file, 
                     cpp_options = list(stan_threads = TRUE, 
                                        O = 3))

sd <- info$stan.data
sd$estimate_gamma <- 1
sd$gamma_prior_mu <- 0
sd$zi_mode <- 2
sd$use_shape_dev <- 0

run.sensitivity.stan(sd, mod, out.file)

### Compare results 
sen.8.models <- c("M1a2_1gamma_mod5_ZI2",
                  "M1a2_1gamma_mod5_ZI2_gammaMu0")

sen.8.out <- extract.results(sen.8.models, params)
saveRDS(sen.8.out, 
        file = "RData//sensitivity_8_out.rds")

## Sensitiviy 9: #  9. gamma_prior_sd = 2.0
info <- readRDS("RData//Richards_HSSM_M1a2_1gamma_mod5_ZI2_stan_info.rds")
out.file <- paste0("Richards_HSSM_M1a2_1gamma_mod5_ZI2_gammaSD2_stan")

mod <- cmdstan_model(mod.file, 
                     cpp_options = list(stan_threads = TRUE, 
                                        O = 3))

sd <- info$stan.data
sd$estimate_gamma <- 1
sd$gamma_prior_mu <- 1
sd$gamma_prior_sd <- 2
sd$zi_mode <- 2
sd$use_shape_dev <- 0

run.sensitivity.stan(sd, mod, out.file)

### Compare results 
sen.9.models <- c("M1a2_1gamma_mod5_ZI2",
                  "M1a2_1gamma_mod5_ZI2_gammaSD2")

sen.9.out <- extract.results(sen.9.models, params)
saveRDS(sen.9.out, 
        file = "RData//sensitivity_9_out.rds")

# 
# sensitivity.table <- data.frame(ID = paste0("sen", seq(0, 10)),
#                                 Sensitivity = c("Parity",
#                                                 "use_pooling_Max = 0",
#                                                 "anchor_mu = qlogis(0.7)",
#                                                 "anchor_mu = qlogis(0.9)",
#                                                 "use_trend_P = 0",
#                                                 "gamma_prior_mu = 0",
#                                                 "gamma_prios_sd = 2.0",
#                                                 "Base",
#                                                 "anchor_mu = qlogis(0.825)",
#                                                 "use_shape_dev = 1",
#                                                 "n_period = 10"),
#                                 ID.2 = c("Parity", "A", "B", "C", "D", "E", "F", "Base", "G", "H", "I"))
# 
# if (sensitivity == "sen1"){
#   stan.data <- create.stan.data(jags.data = jags.data, 
#                                 use_pooling_Max = 0)
# } else if (sensitivity == "sen2"){
#   stan.data <- create.stan.data(jags.data = jags.data, 
#                                 anchor_mu = qlogis(0.7))
# } else if (sensitivity == "sen3"){
#   stan.data <- create.stan.data(jags.data = jags.data, 
#                                 anchor_mu = qlogis(0.9))
# } else if (sensitivity == "sen4"){
#   stan.data <- create.stan.data(jags.data = jags.data, 
#                                 use_trend_P = 0)
# } else if (sensitivity == "sen0"){
#   stan.data <- create.stan.data(jags.data = jags.data, 
#                                 anchor_sd = 0.01, estimate_gamma = 0)
# } else if (sensitivity == "sen5"){
#   stan.data <- create.stan.data(jags.data = jags.data, 
#                                 gamma_prior_mu = 0)
# } else if (sensitivity == "sen6"){
#   stan.data <- create.stan.data(jags.data = jags.data, 
#                                 gamma_prior_sd = 2.0)
# } else if (sensitivity == "sen7"){
#   stan.data <- create.stan.data(jags.data = jags.data)
# } else if (sensitivity == "sen8"){
#   stan.data <- create.stan.data(jags.data = jags.data, 
#                                 anchor_mu = qlogis(0.825))
# } else if (sensitivity == "sen9"){
#   stan.data <- create.stan.data(jags.data = jags.data, 
#                                 use_shape_dev = 1)
# } else if (sensitivity == "sen10"){
#   stan.data <- create.stan.data(jags.data = jags.data, 
#                                 use_shape_dev = 1,
#                                 n_period = 10)
# }
# 
# # Create an inits function
# n_year <- stan.data$jags.data$n.year
# n_observer <- stan.data$jags.data$n.obs.fixed
# 
# file <- file.path("models//model_Richards_HSSM_mod4.stan")
# out.file <- paste0("Richards_HSSM_", model, "_mod4_", sensitivity, "_stan")
# 
# # --- 5. Inspect Results ---
# params.1.stan <- c("S1", "beta_S1", "S2", "P", "phi",
#                    "sigma_proc_P", "Corrected_Est", "Max", "gamma_free[1]")
# 
# # Jags output:
# jm.out <- readRDS("RData/JAGS_Richards_HSSM_M1a2_1968to2026_min60_2026-07-27_NoBUGS.rds")
# 
# # From Jags using regular expression syntax:
# params.1.jags <- c("^S1\\[", "^S1.beta", "^S2\\[", "^P\\[", "^r",
#                    "^sd.proc.P", "^Corrected.Est", "^Max\\[")
# 
# jags.params.idx <- lapply(params.1.jags,
#                           FUN = function(x) grep(x, jm.out$posterior.summary$variable)) %>% unlist()
# 
# summarise_draws(jm.out$jm$samples, 
#                 default_summary_measures(),
#                 default_convergence_measures(),
#                 extra_quantiles = ~quantile2(., probs = c(0.025, 0.975))) %>%
#   as.data.frame() -> jags.convergence.summary
# 
# jags.global.summary <- jags.convergence.summary[jags.params.idx,]
# 
# jags.Nhats <- jags.global.summary[grep("^Corrected.Est", jags.global.summary$variable),] %>%
#   mutate(MCMC = "Jags",
#          year.idx = seq(1, n_year))
# jags.Nhats$start.year <- jags.data$start.years
# 
# # Use samples to compute SE of Nhats
# Nhats_jags <- jm.out$jm$sims.list$Corrected.Est      # or S2
# mean.Nhats_jags <- apply(Nhats_jags, 2, mean) 
# se.Nhats_jags <- apply(Nhats_jags, 2, sd)/sqrt(coda::effectiveSize(Nhats_jags))


# 
# # --- Get Summaries for Specific Global Parameters ---
# stan_global_summary <- fit_stan$summary(
#   variables = params.1.stan,
#   default_summary_measures(), 
#   default_convergence_measures(),
#   extra_quantiles = ~quantile2(., probs = c(0.025, 0.975)))
# 
# stan.Nhats <- stan_global_summary[grep("^Corrected_Est", stan_global_summary$variable),] %>%
#   mutate(MCMC = "Stan",
#          year.idx = seq(1, n_year)) 
# stan.Nhats$start.year <- jags.data$start.years
# 
# Nhats.stan.jags <- rbind(stan.Nhats, jags.Nhats)
# 
# ggplot(Nhats.stan.jags) +
#   geom_point(aes(x = year.idx, y = mean, color = MCMC),
#              alpha = 0.6) +
#   geom_errorbar(aes(x = year.idx, ymin = q2.5, ymax = q97.5)) +
#   theme(legend.position = "top")
# 
# Nhats_stan_summary  <- fit_stan$summary("Corrected_Est")
# mean.Nhats_stan <- Nhats_stan_summary$mean
# se.Nhats_stan <- Nhats_stan_summary$sd / sqrt(Nhats_stan_summary$ess_bulk)
# 
# z <- (mean.Nhats_stan - mean.Nhats_jags) / sqrt(se.Nhats_jags^2 + se.Nhats_stan^2)     # standardised difference
# rel <- (mean.Nhats_stan - mean.Nhats_jags) / mean.Nhats_jags
# 
# round(summary(z), 2)
# sum(abs(z) > 2)                # how many exceed 2 MCSE
# 
# #systematic bias
# mean.Nhats.df <- data.frame(Jags = mean.Nhats_jags,
#                             Stan = mean.Nhats_stan,
#                             start.year = jags.data$start.years)
# 
# title.txt <- paste0(sensitivity, "(", 
#                     sensitivity.table[grep(sensitivity, sensitivity.table$ID), "ID.2"],  "): ",
#                     sensitivity.table[grep(sensitivity, sensitivity.table$ID), "Sensitivity"])
# 
# ggplot(mean.Nhats.df) +
#   geom_point(aes(x = Jags, y = Stan, color = start.year)) +
#   geom_abline(slope = 1, color = "red") +
#   labs(title = title.txt) +
#   xlab("Jags Mean Corrected Est") +
#   ylab("Stan Mean Corrected Est") 
# 
# sum(mean.Nhats_stan > mean.Nhats_jags)                                     # want ~17, not 34 or 0
# binom.test(sum(mean.Nhats_stan > mean.Nhats_jags), length(mean.Nhats_stan))$p.value    # want non-significant
# summary((mean.Nhats_stan - mean.Nhats_jags) / sqrt(se.Nhats_jags^2 + se.Nhats_stan^2))       # want |z| mostly < 2

# If there were divergence - see below:
# dv <- fit_stan$sampler_diagnostics(format = "df")$divergent__
# dr <- fit_stan$draws(format = "df")
# colMeans(dr[dv == 1, c("sigma_proc_P","sigma_proc_Max","sigma_Obs","phi")])
# colMeans(dr[dv == 0, c("sigma_proc_P","sigma_proc_Max","sigma_Obs","phi")])

# colMeans(dr[dv == 1, c("delta[1]","sigma_proc_P","sigma_proc_Max","sigma_Obs","phi")])
# colMeans(dr[dv == 0, c("delta[1]","sigma_proc_P","sigma_proc_Max","sigma_Obs","phi")])
# 
# fit_stan$summary(c("peak_day_decade", "peak_day_slope"))
# fit_stan$summary("peak_day")
# fit_stan$summary("peak_width")
# # 
# pd <- fit_stan$draws("peak_day", format = "matrix")
# P  <- fit_stan$draws("P",        format = "matrix")
# off <- pd - P                                     # mode minus midpoint, per season
# sl <- apply(off, 1, \(o) coef(lm(o ~ stan.data$stan.data$year_values))[2])
# quantile(sl * 10, c(0.025, 0.5, 0.975))           # days per decade of drift
