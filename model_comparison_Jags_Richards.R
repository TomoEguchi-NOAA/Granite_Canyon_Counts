# model_comparison
# Extracts Rhat, LOOIC, and Pareto K statistics from all runs and creates a table
# to compare performance of all models with Richards function
# 

all.vers <- c(2, 15, 17, 18, 5, 16, 19, 20, 1, 13, 21, 22, 3, 12, 23, 24, 10, 11,
              25, 26, 9, 14, 27, 28)

model <- vector(mode = "character", length = length(all.vers))
n.big.Rhat <- n.bad.Pareto <- prop.bad.Pareto <- LOOIC <- vector(mode = "numeric", length = length(all.vers))

min.dur <- 60
years <- c(2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024, 2025)
data.dir <- "RData/V2.1_Feb2025"
max.day <- 100

for (k in 1:length(all.vers)){
  model[k] <- paste0("v", all.vers[k], "a")
  model.name <- paste0("Richards_Nmixture_", model[k]) 
  jags.model <- paste0("models/model_", model.name, ".txt")
  
  out.file.name <- paste0("RData/JAGS_", model.name, 
                          "_1968to", max(years), 
                          "_min", min.dur,
                          "_NoBUGS.rds")
  jm.out <- readRDS(file = out.file.name)
  data.array <- jm.out$jags.input$jags.data$n
  data.array[,2,which(jm.out$jags.input$jags.data$n.station == 1)] <- NA
  data.array[,2,which(jm.out$jags.input$jags.data$n.station == 1)] <- NA
  
  LOOIC.n <- compute.LOOIC(loglik.array = jm.out$jm$sims.list$log.lkhd,
                           data.array = data.array,
                           MCMC.params = MCMC.params)
  
  # There are some (< 0.5%) bad ones. I should look at which ones are not fitting well.
  
  # Compute new Rhat (rank-normalized Rhat as per Vehtari et al. 2021)
  params <- "^VS\\.Fixed|^BF\\.Fixed|^Max|^S1|^S2|^P1|^P2|^K|^OBS\\.RF\\["
  new.Rhat <- rank.normalized.R.hat(jm.out$jm$samples, params)
  #max.Rhat <- lapply(.out[[k]]$jm$Rhat, FUN = max, na.rm = T) %>%
  #  unlist()
  max.new.Rhat.big <- new.Rhat[which(new.Rhat > 1.01)]
  
  n.big.Rhat[k] <- length(max.new.Rhat.big)
  LOOIC[k] <- LOOIC.n$loo.out$estimates["looic", "Estimate"]
  n.bad.Pareto[k] <- sum(LOOIC.n$loo.out$pointwise[,5] > 0.7)
  prop.bad.Pareto[k] <- 100 * (n.bad.Pareto[k]/nrow(LOOIC.n$loo.out$pointwise))
}

out.table <- data.frame(model = model,
                        LOOIC = LOOIC,
                        n.big.Rhat = n.big.Rhat,
                        n.bad.Pareto = n.bad.Pareto,
                        p.bad.Pareto = prop.bad.Pareto)


write.csv(out.table, file = "Data/Richards_model_comp.csv", quote = FALSE,
          row.names = FALSE)
