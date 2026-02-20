#Jags_ModelComparison
#
# Compares Jags models using Richards function using LOOIC, Pareto-K, 
# rank-normalized Rhat (Rhat), and effective sample sizes (ESS).
# 

rm(list = ls())
source("Granite_Canyon_Counts_fcns.R")
library(tidyverse)
library(posterior)
library(ggplot2)
library(bayesplot)
library(loo)
library(rstanarm)

min.dur <- 60
YEAR <- 2026

# Model name IDs
# model.names <- c("1", "2", "3", "4",
#                  "5", "6", "7", "8")

model.names <- c( "6", "7", "8")

# model IDs in the manuscript is in the same order as above but the numbers are
# different:
model.ID <- c(1:length(model.names))

# Select only parameters for the function. 
# ^ means the subsequent letter is the beginning of a string
# \\. means a literal period
# | means "or"
# ^P\\b(?!\\.\\w*) means P followed by nothing, including a period or a white space
# \\[ means a literal bracket
# \\b is a word boundary
# (?!...) is a negative look ahead, which requires the perl = TRUE argument in R's grep.
# params <- "^VS\\.Fixed|^BF\\.Fixed|^Max\\[|^S1\\[|^S2\\[|^P\\b(?!\\.\\w*)|^OBS\\.RF\\["

# For constant Max models (Model 7 and Model 8)
# params.1 <- "^VS\\.Fixed|^BF\\.Fixed|^Max|^S1\\[|^S2\\[|^P\\b(?!\\.\\w*)|^OBS\\.RF\\["

# including hyperparameters and year-specific Max (all models except 7 and 8.
#params.2 <- "^VS\\.Fixed|^BF\\.Fixed|^Max\\[|^S1|^S2|^P|^OBS\\.RF\\["
params.a1 <- "^VS\\.Fixed|^BF\\.Fixed|^Max\\[|^S1|^S2|^P|^alpha\\["

#max.Rhat.big <- list()
prop.big.Rhat <- n.params <- n.big.Rhat <- n.bad.Pareto <- prop.bad.Pareto <- LOOIC <- vector(mode = "numeric", length = length(model.names))

min.ESS.bulk <- min.ESS.tail <- vector(mode = "numeric", length = length(model.names))

k <- 6
for (k in 1:length(model.names)){
  .out <- readRDS(paste0("RData/JAGS_Richards_HSSM_v", 
                         model.names[k], "a1_1968to", YEAR, "_min", 
                         min.dur, "_NoBUGS.rds"))
  
  new.Rhat <- rank.normalized.R.hat(.out$jm$samples, 
                                    params = params.a1, 
                                    MCMC.params = .out$MCMC.params)
  
  n.params[k] <- length(new.Rhat)
  #max.Rhat <- lapply(.out[[k]]$jm$Rhat, FUN = max, na.rm = T) %>%
  #  unlist()
  max.Rhat.big <- new.Rhat[which(new.Rhat > 1.01)]
  
  # data.array <- .out$jags.input$jags.data$n
  # data.array[,2,which(.out$jags.input$jags.data$n.station == 1)] <- NA
  # data.array[,2,which(.out$jags.input$jags.data$n.station == 1)] <- NA
  
  LOOIC.n <- compute.LOOIC(loglik.array = .out$jm$sims.list$log.lkhd,
                           #data.array = data.array,
                           MCMC.params = .out$MCMC.params)
  
  n.big.Rhat[k] <- length(max.Rhat.big)
  
  prop.big.Rhat[k] <- 100 * (length(max.Rhat.big)/n.params[k])
  
  LOOIC[k] <- LOOIC.n$loo.out$estimates["looic", "Estimate"]
  
  n.bad.Pareto[k] <- sum(LOOIC.n$loo.out$pointwise[,5] > 0.7)
  
  prop.bad.Pareto[k] <- 100 * (n.bad.Pareto[k]/nrow(LOOIC.n$loo.out$pointwise))
  
  summary.posterior <- .out$posterior.summary #summarise_draws(post)
  
  summary.posterior %>%
    select(variable, ess_bulk) %>%
    na.omit() %>%
    arrange(ess_bulk) -> ESS.bulk
  
  min.ESS.bulk[k] <- min(ESS.bulk$ess_bulk)
  
  summary.posterior %>%
    select(variable, ess_tail) %>%
    na.omit() %>%
    arrange(ess_tail) -> ESS.tail
  
  min.ESS.tail[k] <- min(ESS.tail$ess_tail)
  
}

out.table <- data.frame(model = model.ID,
                        LOOIC = LOOIC,
                        n.params = n.params,
                        p.big.Rhat = prop.big.Rhat,
                        #LOOIC = LOOIC,
                        n.big.Rhat = n.big.Rhat,
                        n.bad.Pareto = n.bad.Pareto,
                        p.bad.Pareto = prop.bad.Pareto,
                        min.ESS.bulk = min.ESS.bulk,
                        min.ESS.tail = min.ESS.tail) %>%  
  arrange(desc(min.ESS.bulk))  %>%
  mutate(dLOOIC = LOOIC - min(LOOIC[min.ESS.bulk == max(min.ESS.bulk)])) %>%
  
  select(model, dLOOIC, n.params, n.big.Rhat, 
         p.big.Rhat, n.bad.Pareto, p.bad.Pareto, 
         min.ESS.bulk, min.ESS.tail, LOOIC)

saveRDS(out.table,
        file = paste0("RData/Richards_ModelComparison_", YEAR, ".rds"))
#d.t.2 <- Sys.time() - t.2 # 1.68 min

