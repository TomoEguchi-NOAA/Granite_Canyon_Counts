# Running all models in STAN
# 
#  Collects Stan runs and compares results. Do posterior predictive checks on
#  all models. 
#  
#  Run all models in Run_StanRichards_AllModels.R. It will save results
#  in .rds files. 
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

# models <- c("M1a1_1gamma_mod4", "M2a1_1gamma_mod4", "M3a1_1gamma_mod4", "M4a1_1gamma_mod4",
#             "M1a2_1gamma_mod4", "M2a2_1gamma_mod4", "M3a2_1gamma_mod4", "M4a2_1gamma_mod4",
#             "M1a2_1gamma_mod5_ZI1", "M1a2_1gamma_mod5_ZI2")

# It's been found that the Negative-Binomial models (a2s) are better than the 
# Poisson models (a1s). Also, among the Negative Binomial models, "M1" (S1 and 
# S2 are season-specific) was the best. Also having a gamma parameter, shared
# by the both sides of the curve helps. So, The Zero-Inflated models (ZIs) are
# compared to the best non-ZI model (M1a2_1gamma_mod4). "mod4" refers to the 4th
# modification to the original model. mod5 incorporated a switch to add the Zero-
# Inflated part, the rest is identical to mod4. 

# sen9 refers to the sensitivity run to check Periodic deviation from the curve 
# use_shape_dev = 1 in the mod4 model.
models <- c("M1a2_1gamma_mod4", "M1a2_1gamma_mod4_sen9",
            "M1a2_1gamma_mod5_ZI1", "M1a2_1gamma_mod5_ZI2")

params.a1.stan <- c("S1", "S2", "P", 
                    "sigma_proc_P", "Corrected_Est", "Max", "log_N_latent",
                    "gamma_free", "peak_day_decade", "beta1_P")

params.a2.stan <- c("S1", "S2", "P", "phi",
                    "sigma_proc_P", "Corrected_Est", "Max", "log_N_latent",
                    "gamma_free", "peak_day_decade", "beta1_P")

diag.summary <- LOO.out <- stan.global.summary <- ppc.res <- list()
zi_mode <- 1
k <- 3

for (k in 1:length(models)) {
  out.file <- paste0("Richards_HSSM_", models[k], "_stan")
  
  #mod <- cmdstan_model(file)
  fit_stan <- readRDS(paste0("RData//", out.file, ".rds"))
  
  if (length(grep("a1", models[k])) > 0) params <- params.a1.stan
  if (length(grep("a2", models[k])) > 0) params <- params.a2.stan
  if (length(grep("ZI1", models[k])) > 0) params <- c(params, "zi_a")
  if (length(grep("ZI2", models[k])) > 0) params <- c(params, "zi_a", "zi_b")

  # The default zi_mode = 1. So, mod4 runs have to set zi_mode = 0
  if (length(grep("ZI", models[k])) == 0) zi_mode = 0
  
  stan.data <- create.stan.data(jags.data = jags.data, zi_mode = zi_mode)
  
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
  
  # plots$resid_day is the same thing visually, and it's the figure to put in the 
  # paper. A flat loess through residuals against day-of-season says the curve 
  # captures the shape; systematic curvature — a dip at the peak, say — 
  # would say it doesn't.
  
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
 
   ppc.res[[k]] <- ppc_richards(fit_stan, stan.data$stan.data, n_draws = 1500)
  
  # res$pvalues
  # res$autocorr
  # res$plots$resid_day
}

### find the best model:

#### For one model:
# stan.global.summary[[1]] %>% filter(variable == "sigma_proc_P")
# 
# sd_y <- sd(jags.data$year.index)
# b1 <- as.vector(fit_stan$draws("beta1_P[1]", format = "matrix"))
# sP <- as.vector(fit_stan$draws("sigma_proc_P", format = "matrix"))
# R  <- (abs(b1)*sd_y)^2 / ((abs(b1)*sd_y)^2 + sP^2)
# quantile(R, c(0.025, 0.5, 0.975))
#   
# ppc.res[[1]]$autocorr
# lapply(ppc.res, function(x) x$autocorr$observed)
# 
# ppc.res[[1]]$plots$resid_day
####

diag.n_divergent.df <- do.call(rbind,
                               lapply(diag.summary, 
                                      FUN = function(x) x$num_divergent)) %>%
  as.data.frame()

# if any of the entries in the above dataframe is >0, uncomment the following chunk:
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
  mutate(dLOOIC = Estimate- min(Estimate)) %>% #-> tmp
  arrange(by = dLOOIC) -> LOOIC.df

# ZI2 is the best according to LOOIC - but other measures need to be considered
coverage.df <- do.call(rbind, lapply(ppc.res, FUN = function(x) x$plots$coverage)) %>%
  as.data.frame() %>%
  rename(Coverage = V1)
coverage.df$Model <- as.vector(models)
# They are about the same

# all a2 models are good.
LOOIC.df %>%
  left_join(coverage.df, by = "Model") -> LOOIC.coverage.df

# if there are many models, select top three and provide the indices here:
# loo::loo_compare(LOO.out[[5]], LOO.out[[9]], LOO.out[[10]]) %>%
#   data.frame() -> loo.comp.df

loo::loo_compare(LOO.out[[1]], LOO.out[[2]], LOO.out[[3]],
                 LOO.out[[4]]) %>%
  data.frame() -> loo.comp.df

#best.model <- 10  # when running all models
#best.model <- 2   # when running just ZI models
best.model <- 4    # Two ZI models and non-ZI M1a2_1gamma
pk <- LOO.out[[best.model]]$diagnostics$pareto_k
which(pk > 0.7)
if (length(which(pk>0.7)) > 0){
  sd <- stan.data$stan.data
  data.frame(day = sd$day_idx, 
             year = sd$year_idx, 
             n = sd$n)[which(pk > 0.7), ]
  
}

pvalues.df <- do.call(rbind, lapply(ppc.res, FUN = function(x) x$pvalues)) %>%
  as.data.frame()
pvalues.df$Model <- rep(models, each = 6)

# compare Corrected Estimates among the models:
Corrected.Est <- list()
for (k in 1:length(models)){
  Corrected.Est[[k]] <- stan.global.summary[[k]] %>%
    filter(str_detect(variable, "Corrected_Est")) %>%
    mutate(model = models[k],
           start.years = start.years)
  
}

Corrected.Est.df <- do.call(rbind, Corrected.Est)
ggplot(Corrected.Est.df) +
  geom_pointrange(aes(x = start.years,
                      y = median,
                      ymin = q2.5,
                      ymax = q97.5,
                      color = model)) +
  theme(legend.position = "top")

#### Claude's suggestion about PPS results:
res <- ppc.res[[best.model]]
ps <- res$setup
sd <- ps$sd
zobs <- sd$n == 0
zrep <- colMeans(ps$y_rep == 0)

## sanity check:
## 0.138 is the proportion of zeros (observed = 0.138 in ppc.res[[3]]$pvalues)
stan.global.summary[[best.model]] %>%
  filter(str_detect(variable, "zi_a")) -> za.summary

stan.global.summary[[best.model]] %>%
  filter(str_detect(variable, "zi_b")) -> zb.summary

stan.global.summary[[best.model]] %>%
  filter(str_detect(variable, "phi")) -> phi.summary

lk  <- log(ps$kappa[1, ])
pi_ <- plogis(za.summary$mean + zb.summary$mean * (lk - mean(lk)))
mean(pi_ + (1 - pi_) * dnbinom(0, mu = ps$kappa[1, ], size = phi.summary$mean))   # vs 0.138
##

## Check auto correlation:
# 
ppc.autocorr <- lapply(ppc.res, FUN = function(x){
  return(data.frame(observed = x$autocorr$observed,
                    rep_median = x$autocorr$rep_median,
                    rep_low = x$autocorr$rep_interval[1],
                    rep_high = x$autocorr$rep_interval[2],
                    p_value = x$autocorr$p_value))
}) %>% do.call(rbind, .) %>%
  mutate(model = models) 
  
ppc.res[[3]]$autocorr



chk <- function(g, lab) {
  data.frame(group = lab, bin = levels(cut(g, 5)),
             obs = tapply(zobs, cut(g, 5), mean),
             pred = tapply(zrep, cut(g, 5), mean))
}
rbind(chk(sd$day_idx, "day of season"),
      chk(sd$watch_length, "watch length"),
      chk(sd$bf, "Beaufort"))

#k <- 10
out.file <- paste0("Richards_HSSM_", models[best.model], "_stan")
#mod <- cmdstan_model(file)
fit_stan <- readRDS(paste0("RData//", out.file, ".rds"))

lN <- fit_stan$draws("log_N_latent", format = "draws_matrix")
tail_frac <- sapply(1:34, function(y) {
  cols_all  <- sprintf("log_N_latent[%d,%d]", 1:100, y)
  cols_late <- sprintf("log_N_latent[%d,%d]", 77:100, y)
  mean(rowSums(exp(lN[, cols_late])) / rowSums(exp(lN[, cols_all])))
})
round(quantile(tail_frac, c(0, .5, 1)), 3)

names(tail_frac) <- jags.data$start.years          # or whatever your season labels are
sort(round(tail_frac, 3), decreasing = TRUE)[1:8]

y25 <- which(jags.data$start.years == 2025)
range(sd$day_idx[sd$year_idx == y25])              # does effort span the season?
length(unique(sd$day_idx[sd$year_idx == y25]))     # how many distinct days surveyed
sum(sd$watch_length[sd$year_idx == y25])           # total effort vs other seasons
fit_stan$summary("P")$mean[y25]                    # fitted peak day
fit_stan$summary("S1")$mean[y25]; fit_stan$summary("S2")$mean[y25]
# Finding that a large part of abundance (21%) comes from the unobserved days (days 1 - 36).

obs_frac <- sapply(1:34, function(y) {
  dr <- range(sd$day_idx[sd$year_idx == y])
  cols <- sprintf("log_N_latent[%d,%d]", 1:100, y)
  cin  <- sprintf("log_N_latent[%d,%d]", dr[1]:dr[2], y)
  mean(rowSums(exp(lN[, cin])) / rowSums(exp(lN[, cols])))
})
names(obs_frac) <- jags.data$start.years
sort(round(obs_frac, 3))[1:10]
Laake.Run.Date <- "2026-02-23"
Laake.abundance <- read.csv(file = paste0("Data//all_estimates_Laake_2026", 
                                              "_", Laake.Run.Date, ".csv")) %>%
  mutate(LCL = CL.low,
         UCL = CL.high) %>% 
  na.omit()
Eguchi.abundance <- fit_stan$summary("Corrected_Est")
NEguchi <- Eguchi.abundance$median
NLaake <- Laake.abundance$Nhat
plot(obs_frac, (NEguchi - NLaake) / NLaake)   # the relationship that matters
# There is no relationship

asc <- sapply(1:34, function(y) {
  d0 <- min(sd$day_idx[sd$year_idx == y])
  cols <- sprintf("log_N_latent[%d,%d]", 1:100, y)
  cpre <- sprintf("log_N_latent[%d,%d]", 1:(d0-1), y)
  mean(rowSums(exp(lN[, cpre])) / rowSums(exp(lN[, cols])))
})
plot(asc, (NEguchi - NLaake) / NLaake)


s2 <- fit_stan$draws("S2", format = "draws_matrix")
data.frame(season = jags.data$start.years,
           S2 = colMeans(s2), sd = apply(s2, 2, sd))[c(y25, order(-tail_frac)[1:5]), ]
fit_stan$summary("mu_S2")   # population mean it shrinks toward

raw90 <- sapply(1:34, function(y) {
  cols <- sprintf("log_N_latent[%d,%d]", 1:100, y)
  c90  <- sprintf("log_N_latent[%d,%d]", 1:90,  y)
  mean(rowSums(exp(lN[, c90])) / rowSums(exp(lN[, cols])))
})
summary(1 - raw90)                    # % of each estimate coming from days 91-100
cor(1 - raw90, (NEguchi - NLaake) / NLaake)

#jm.out <- readRDS()
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
# params.1.jags <- c("^S1\\[", "^S1.beta", "^S2\\[", "^P\\[", "^r",
#                    "^sd.proc.P", "^Corrected.Est", "^Max\\[")
# 
# jags.params.idx <- lapply(params.1.jags,
#                      FUN = function(x) grep(x, jm.out$posterior.summary$variable)) %>% unlist()
# 
# summarise_draws(jm.out$jm$samples, 
#                 default_summary_measures(),
#                 default_convergence_measures(),
#                 extra_quantiles = ~quantile2(., probs = c(0.025, 0.975))) %>%
#   as.data.frame() -> jags.convergence.summary
# 
# jags.global.summary <- jags.convergence.summary[jags.params.idx,]
# 
# stan.S1 <- stan_global_summary[grep("^S1\\[", stan_global_summary$variable),] %>%
#   mutate(MCMC = "Stan",
#          year.idx = seq(1, n_year))
# jags.S1 <- jags.global.summary[grep("^S1\\[", jags.global.summary$variable),] %>%
#   mutate(MCMC = "Jags",
#          year.idx = seq(1, n_year))
# 
# S1.stan.jags <- rbind(stan.S1, jags.S1)
# 
# ggplot(S1.stan.jags) +
#   geom_point(aes(x = year.idx, y = mean, color = MCMC)) +
#   geom_errorbar(aes(x = year.idx, ymin = q2.5, ymax = q97.5)) +
#   theme(legend.position = "top")
# 
# s_jags <- jm.out$jm$sims.list$S1      # or S2
# m_j <- apply(s_jags, 2, mean); se_j <- apply(s_jags, 2, sd)/sqrt(coda::effectiveSize(s_jags))
# sj  <- fit_stan$summary("S1")
# m_s <- sj$mean; se_s <- sj$sd / sqrt(sj$ess_bulk)
# 
# z <- (m_s - m_j) / sqrt(se_j^2 + se_s^2)     # standardised difference
# rel <- (m_s - m_j) / m_j
# 
# round(summary(z), 2)
# sum(abs(z) > 2)                # how many exceed 2 MCSE
# binom.test(sum(m_s > m_j), length(m_s))$p.value   # sign test for systematic bias
# plot(m_j, m_s); abline(0, 1, col = "red")
# 
# sum(m_s > m_j)                                     # want ~17, not 34 or 0
# binom.test(sum(m_s > m_j), length(m_s))$p.value    # want non-significant
# summary((m_s - m_j) / sqrt(se_j^2 + se_s^2))       # want |z| mostly < 2
# 
# 
# stan.Nhats <- stan_global_summary[grep("^Corrected_Est", stan_global_summary$variable),] %>%
#   mutate(MCMC = "Stan",
#          year.idx = seq(1, n_year))
# jags.Nhats <- jags.global.summary[grep("^Corrected.Est", jags.global.summary$variable),] %>%
#   mutate(MCMC = "Jags",
#          year.idx = seq(1, n_year))
# 
# Nhats.stan.jags <- rbind(stan.Nhats, jags.Nhats)
# 
# ggplot(Nhats.stan.jags) +
#   geom_point(aes(x = year.idx, y = mean, color = MCMC)) +
#   geom_errorbar(aes(x = year.idx, ymin = q2.5, ymax = q97.5)) +
#   theme(legend.position = "top")
# 
# s_jags <- jm.out$jm$sims.list$Corrected.Est      # or S2
# m_j <- apply(s_jags, 2, mean); se_j <- apply(s_jags, 2, sd)/sqrt(coda::effectiveSize(s_jags))
# sj  <- fit_stan$summary("Corrected_Est")
# m_s <- sj$mean; se_s <- sj$sd / sqrt(sj$ess_bulk)
# 
# z <- (m_s - m_j) / sqrt(se_j^2 + se_s^2)     # standardised difference
# rel <- (m_s - m_j) / m_j
# 
# round(summary(z), 2)
# sum(abs(z) > 2)                # how many exceed 2 MCSE
# #systematic bias
# plot(m_j, m_s); abline(0, 1, col = "red")
# 
# sum(m_s > m_j)                                     # want ~17, not 34 or 0
# binom.test(sum(m_s > m_j), length(m_s))$p.value    # want non-significant
# summary((m_s - m_j) / sqrt(se_j^2 + se_s^2))       # want |z| mostly < 2
# 
# fit_stan$summary("logit_p0")   # expect mean ~1.40, above the 1.39 anchor
# 
# # same comparison on Raw_Est, which has no correction factor
# s_jags <- jm.out$jm$sims.list$Raw.Est      # or S2
# m_j <- apply(s_jags, 2, mean); se_j <- apply(s_jags, 2, sd)/sqrt(coda::effectiveSize(s_jags))
# sj  <- fit_stan$summary("Raw_Est")
# m_s <- sj$mean; se_s <- sj$sd / sqrt(sj$ess_bulk)
# 
# z <- (m_s - m_j) / sqrt(se_j^2 + se_s^2)     # standardised difference
# rel <- (m_s - m_j) / m_j
# 
# round(summary(z), 2)
# sum(abs(z) > 2)                # how many exceed 2 MCSE
# #systematic bias
# plot(m_j, m_s); abline(0, 1, col = "red")
# 
# sum(m_s > m_j)                                     # want ~17, not 34 or 0
# binom.test(sum(m_s > m_j), length(m_s))$p.value    # want non-significant
# summary((m_s - m_j) / sqrt(se_j^2 + se_s^2))       # want |z| mostly < 2
# 
# mean(fit_stan$summary("Raw_Est")$mean) / mean(apply(jm.out$jm$sims.list$Raw.Est, 2, mean))
# sum(m_s > m_j)
# 
# # and check the two correction factors directly
# fit_stan$summary("corr_factor")$mean
# 
# fit_stan$summary("sigma_Obs")$mean
# median(jm.out$jm$sims.list$sd.obs)     # was 0.243
# 
# 
# # implied mean detection in each fit
# mean(plogis(jm.out$jm$sims.list$alpha))           # JAGS, before covariates
# fit_stan$summary(c("BF_Fixed","VS_Fixed","logit_p0","sigma_Obs"))$mean
# 
# # vs JAGS BF.Fixed / VS.Fixed
# summary(stan_data$bf); summary(stan_data$vs)   # raw ranges? vs shifted?
# 
# # --- Get Summaries for Year-Specific Parameters ---
# # This will print the estimates for every year automatically
# richards_summary <- fit_stan$summary(
#   variables = c("Max", "S1", "S2", "K", "P1"),
#   "mean", "median", "quantile~0.025", "quantile~0.975", "rhat"
# )
# print(richards_summary)
# 
# library(bayesplot)
# library(ggplot2)
# 
# # Extract draws for the parameters you want to check
# # (Using a few years of Max and S1 as an example)
# draws_diagnostic <- fit_stan$draws(c("mean_prob", "Max[1]", "Max[10]", "S1[1]"))
# 
# # Generate trace plots
# mcmc_trace(draws_diagnostic) +
#   theme_minimal() +
#   labs(title = "MCMC Chain Trace Plots")
# 
# # Extract all S1 year estimates
# s1_draws <- fit_stan$draws("S1")
# 
# # Plot the intervals across all years sequentially
# mcmc_intervals(s1_draws) +
#   theme_minimal() +
#   labs(title = "Posterior Estimates for S1 Across Years",
#        x = "Parameter Value", y = "Year Index")
# 
# # Compare the fixed effect coefficients
# covariate_draws <- fit_stan$draws(c("BF_Fixed", "VS_Fixed"))
# 
# mcmc_areas(covariate_draws, prob = 0.95) +
#   theme_minimal() +
#   labs(title = "Posterior Distributions of Detection Covariates",
#        x = "Effect Size")
# 
# # 1. Extract the mean_N matrix from the posterior
# # format = "draws_matrix" converts it into a clean math grid
# mean_N_posterior <- fit_stan$draws("mean_N", format = "draws_matrix")
# 
# # 2. Pick a specific year to visualize (e.g., Year 5)
# target_year <- 5
# day_columns <- paste0("mean_N[", 1:n.days, ",", target_year, "]")
# 
# # Extract only the columns belonging to that year
# year_data <- mean_N_posterior[, day_columns]
# 
# # 3. Calculate the mean and 95% Credible Interval for each day
# trajectory <- data.frame(
#   Day = 1:n.days,
#   Mean = colMeans(year_data),
#   Lower = apply(year_data, 2, quantile, probs = 0.025),
#   Upper = apply(year_data, 2, quantile, probs = 0.975)
# )
# 
# # 4. Plot the migration curve!
# ggplot(trajectory, aes(x = Day, y = Mean)) +
#   geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "blue", alpha = 0.15) +
#   geom_line(color = "blue", size = 1) +
#   theme_minimal() +
#   labs(title = paste("Estimated Gray Whale Migration Curve: Year", target_year),
#        x = "Days into Season",
#        y = "Estimated Expected Abundance (mean_N)")



##############################################
# 
# # --- 5. Inspect Results ---
# params.1.stan <- c("S1", "beta_S1", "S2", "P", "phi",
#                    "sigma_proc_P", "Corrected_Est", "Max")
# 
# # --- Get Summaries for Specific Global Parameters ---
# stan_global_summary <- fit_stan$summary(
#   variables = params.1.stan,
#   default_summary_measures(), 
#   default_convergence_measures(),
#   extra_quantiles = ~quantile2(., probs = c(0.025, 0.975)))
# 
# print(stan_global_summary)
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
# stan.S1 <- stan_global_summary[grep("^S1\\[", stan_global_summary$variable),] %>%
#   mutate(MCMC = "Stan",
#          year.idx = seq(1, n_year))
# jags.S1 <- jags.global.summary[grep("^S1\\[", jags.global.summary$variable),] %>%
#   mutate(MCMC = "Jags",
#          year.idx = seq(1, n_year))
# 
# S1.stan.jags <- rbind(stan.S1, jags.S1)
# 
# ggplot(S1.stan.jags) +
#   geom_point(aes(x = year.idx, y = mean, color = MCMC)) +
#   geom_errorbar(aes(x = year.idx, ymin = q2.5, ymax = q97.5)) +
#   theme(legend.position = "top")
# 
# s_jags <- jm.out$jm$sims.list$S1      # or S2
# m_j <- apply(s_jags, 2, mean); se_j <- apply(s_jags, 2, sd)/sqrt(coda::effectiveSize(s_jags))
# sj  <- fit_stan$summary("S1")
# m_s <- sj$mean; se_s <- sj$sd / sqrt(sj$ess_bulk)
# 
# z <- (m_s - m_j) / sqrt(se_j^2 + se_s^2)     # standardised difference
# rel <- (m_s - m_j) / m_j
# 
# round(summary(z), 2)
# sum(abs(z) > 2)                # how many exceed 2 MCSE
# binom.test(sum(m_s > m_j), length(m_s))$p.value   # sign test for systematic bias
# plot(m_j, m_s); abline(0, 1, col = "red")
# 
# sum(m_s > m_j)                                     # want ~17, not 34 or 0
# binom.test(sum(m_s > m_j), length(m_s))$p.value    # want non-significant
# summary((m_s - m_j) / sqrt(se_j^2 + se_s^2))       # want |z| mostly < 2
# 
# 
# stan.Nhats <- stan_global_summary[grep("^Corrected_Est", stan_global_summary$variable),] %>%
#   mutate(MCMC = "Stan",
#          year.idx = seq(1, n_year))
# jags.Nhats <- jags.global.summary[grep("^Corrected.Est", jags.global.summary$variable),] %>%
#   mutate(MCMC = "Jags",
#          year.idx = seq(1, n_year))
# 
# Nhats.stan.jags <- rbind(stan.Nhats, jags.Nhats)
# 
# ggplot(Nhats.stan.jags) +
#   geom_point(aes(x = year.idx, y = mean, color = MCMC)) +
#   geom_errorbar(aes(x = year.idx, ymin = q2.5, ymax = q97.5)) +
#   theme(legend.position = "top")
# 
# s_jags <- jm.out$jm$sims.list$Corrected.Est      # or S2
# m_j <- apply(s_jags, 2, mean); se_j <- apply(s_jags, 2, sd)/sqrt(coda::effectiveSize(s_jags))
# sj  <- fit_stan$summary("Corrected_Est")
# m_s <- sj$mean; se_s <- sj$sd / sqrt(sj$ess_bulk)
# 
# z <- (m_s - m_j) / sqrt(se_j^2 + se_s^2)     # standardised difference
# rel <- (m_s - m_j) / m_j
# 
# round(summary(z), 2)
# sum(abs(z) > 2)                # how many exceed 2 MCSE
# #systematic bias
# plot(m_j, m_s); abline(0, 1, col = "red")
# 
# sum(m_s > m_j)                                     # want ~17, not 34 or 0
# binom.test(sum(m_s > m_j), length(m_s))$p.value    # want non-significant
# summary((m_s - m_j) / sqrt(se_j^2 + se_s^2))       # want |z| mostly < 2
# 
# fit_stan$summary("logit_p0")   # expect mean ~1.40, above the 1.39 anchor
# 
# # same comparison on Raw_Est, which has no correction factor
# s_jags <- jm.out$jm$sims.list$Raw.Est      # or S2
# m_j <- apply(s_jags, 2, mean); se_j <- apply(s_jags, 2, sd)/sqrt(coda::effectiveSize(s_jags))
# sj  <- fit_stan$summary("Raw_Est")
# m_s <- sj$mean; se_s <- sj$sd / sqrt(sj$ess_bulk)
# 
# z <- (m_s - m_j) / sqrt(se_j^2 + se_s^2)     # standardised difference
# rel <- (m_s - m_j) / m_j
# 
# round(summary(z), 2)
# sum(abs(z) > 2)                # how many exceed 2 MCSE
# #systematic bias
# plot(m_j, m_s); abline(0, 1, col = "red")
# 
# sum(m_s > m_j)                                     # want ~17, not 34 or 0
# binom.test(sum(m_s > m_j), length(m_s))$p.value    # want non-significant
# summary((m_s - m_j) / sqrt(se_j^2 + se_s^2))       # want |z| mostly < 2
# 
# mean(fit_stan$summary("Raw_Est")$mean) / mean(apply(jm.out$jm$sims.list$Raw.Est, 2, mean))
# sum(m_s > m_j)
# 
# # and check the two correction factors directly
# fit_stan$summary("corr_factor")$mean
# 
# fit_stan$summary("sigma_Obs")$mean
# median(jm.out$jm$sims.list$sd.obs)     # was 0.243
# 
# 
# # implied mean detection in each fit
# mean(plogis(jm.out$jm$sims.list$alpha))           # JAGS, before covariates
# fit_stan$summary(c("BF_Fixed","VS_Fixed","logit_p0","sigma_Obs"))$mean
# 
# # vs JAGS BF.Fixed / VS.Fixed
# summary(stan_data$bf); summary(stan_data$vs)   # raw ranges? vs shifted?
# 
# # --- Get Summaries for Year-Specific Parameters ---
# # This will print the estimates for every year automatically
# richards_summary <- fit_stan$summary(
#   variables = c("Max", "S1", "S2", "K", "P1"),
#   "mean", "median", "quantile~0.025", "quantile~0.975", "rhat"
# )
# print(richards_summary)
# 
# library(bayesplot)
# library(ggplot2)
# 
# # Extract draws for the parameters you want to check
# # (Using a few years of Max and S1 as an example)
# draws_diagnostic <- fit_stan$draws(c("mean_prob", "Max[1]", "Max[10]", "S1[1]"))
# 
# # Generate trace plots
# mcmc_trace(draws_diagnostic) +
#   theme_minimal() +
#   labs(title = "MCMC Chain Trace Plots")
# 
# # Extract all S1 year estimates
# s1_draws <- fit_stan$draws("S1")
# 
# # Plot the intervals across all years sequentially
# mcmc_intervals(s1_draws) +
#   theme_minimal() +
#   labs(title = "Posterior Estimates for S1 Across Years",
#        x = "Parameter Value", y = "Year Index")
# 
# # Compare the fixed effect coefficients
# covariate_draws <- fit_stan$draws(c("BF_Fixed", "VS_Fixed"))
# 
# mcmc_areas(covariate_draws, prob = 0.95) +
#   theme_minimal() +
#   labs(title = "Posterior Distributions of Detection Covariates",
#        x = "Effect Size")
# 
# # 1. Extract the mean_N matrix from the posterior
# # format = "draws_matrix" converts it into a clean math grid
# mean_N_posterior <- fit_stan$draws("mean_N", format = "draws_matrix")
# 
# # 2. Pick a specific year to visualize (e.g., Year 5)
# target_year <- 5
# day_columns <- paste0("mean_N[", 1:n.days, ",", target_year, "]")
# 
# # Extract only the columns belonging to that year
# year_data <- mean_N_posterior[, day_columns]
# 
# # 3. Calculate the mean and 95% Credible Interval for each day
# trajectory <- data.frame(
#   Day = 1:n.days,
#   Mean = colMeans(year_data),
#   Lower = apply(year_data, 2, quantile, probs = 0.025),
#   Upper = apply(year_data, 2, quantile, probs = 0.975)
# )
# 
# # 4. Plot the migration curve!
# ggplot(trajectory, aes(x = Day, y = Mean)) +
#   geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "blue", alpha = 0.15) +
#   geom_line(color = "blue", size = 1) +
#   theme_minimal() +
#   labs(title = paste("Estimated Gray Whale Migration Curve: Year", target_year),
#        x = "Days into Season",
#        y = "Estimated Expected Abundance (mean_N)")
# 
# 
# # --- Posterior Predictive Simulation Loop ---
# library(bayesplot)
# library(posterior)
# 
# cat("Extracting posterior draws...\n")
# # Extract parameters as flat iteration-by-column matrices
# mean_N_mat <- as.matrix(fit_stan$draws("mean_N", format = "matrix"))
# r_mat      <- as.matrix(fit_stan$draws("r", format = "matrix"))
# p_mat      <- as.matrix(fit_stan$draws("obs_prob", format = "matrix"))
# 
# S <- nrow(p_mat)       # Total MCMC iterations saved
# V <- nrow(flat_df)     # Total unique data observations
# yrep <- matrix(0, nrow = S, ncol = V)
# 
# # Pre-calculate column lookups so R doesn't have to search text inside the loop
# mean_N_cols <- match(paste0("mean_N[", flat_df$day_idx, ",", flat_df$year_idx, "]"), colnames(mean_N_mat))
# r_cols      <- match(paste0("r[", flat_df$day_idx, ",", flat_df$year_idx, "]"), colnames(r_mat))
# 
# cat("Simulating replicated datasets...\n")
# for (s in 1:S) {
#   # Grab the parameters for this specific MCMC step across all observations
#   mu_val  <- mean_N_mat[s, mean_N_cols]
#   phi_val <- r_mat[s, r_cols]
#   p_val   <- p_mat[s, ]
#   
#   # Structural safety floor
#   mu_val[mu_val < 1e-6] <- 1e-6
#   
#   # Step 1: Simulate the true unobserved abundance (N) for this iteration
#   # (NOTE: For your Poisson models, swap this line to: N_sim <- rpois(V, lambda = mu_val))
#   N_sim <- rnbinom(V, size = phi_val, mu = mu_val)
#   
#   # Step 2: Simulate the observation process to get final whale counts
#   yrep[s, ] <- rbinom(V, size = N_sim, prob = p_val)
# }
# 
# # Define your actual raw observed data vector
# y <- flat_df$n
# 
# # We look at the first 50 simulated datasets to keep the plot clean
# ppc_dens_overlay(y, yrep[1:50, ]) +
#   theme_minimal() +
#   labs(title = "Posterior Predictive Overlays",
#        x = "Whale Count Value", y = "Density")
# 
# # 1. Does the model accurately predict the maximum number of whales seen in a single survey?
# plot_max <- ppc_stat(y, yrep, stat = "max") + 
#   labs(title = "Checking Maximum Values")
# 
# # 2. Does the model accurately handle zero-inflation (proportion of zero counts)?
# prop_zero <- function(x) mean(x == 0)
# plot_zero <- ppc_stat(y, yrep, stat = "prop_zero") + 
#   labs(title = "Checking Proportion of Zeros")
# 
# # View them side by side
# library(gridExtra)
# grid.arrange(plot_max, plot_zero, ncol = 2)
# 
# # Because looking at thousands of rows at once is overwhelming, 
# # let's look at a subset of 100 observations to see the fit clearly.
# subset_idx <- 1:100
# 
# ppc_intervals(y[subset_idx], yrep[, subset_idx]) +
#   theme_minimal() +
#   labs(title = "Predictive Intervals vs. Observed Counts (Subset)",
#        x = "Observation Index", y = "Whale Count")