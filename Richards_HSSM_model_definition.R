#Richards_HSSM_model_definitions
# A function to create JAGS models for Richards_HSSM abundance estimation
# 
# This function creates a text file of a model with user input of which
# parameters of Richards function are time-specific. If any parameter is fixed 
# to a constant, provide the value. 
# 
# Richards function:
# 
# # M1 <- (1 + (2 * exp(K) - 1) * exp((1/S1) * (P - d))) ^ (-1/exp(K))
# M2 <- (1 + (2 * exp(K) - 1) * exp((1/S2) * (P - d))) ^ (-1/exp(K))
# N <- min + (max - min) * (M1 * M2)
# 
# Provide either "time", a constant, or the parameter name to K, S1, S2, and P. 
# If "time", parameters will be specified as time-specific. For example, K = "time"
# will result in K[y]. To fix the parameter, provide K = 1, for example. If a 
# constant but unknown over time, then specify as K = "K", for example. No K[y]
# is allowed at the moment. No priors are set for time-varying K in the code.
# 
# The current convention is:
# M5: P = "time", Max = "time", S1 = "time", S2 = "time", K = 1
# M6: P = "time", Max = "time", S1 = "S1", S2 = "S2", K = 1
# M7: P = "time", Max = "time", S1 = "time", S2 = "S2", K = 1
# M8: P = "time", Max = "time", S2 = "time", S1 = "S1", K = 1
# 
# Poisson likelihood (Poisson): a1
# Negative binomial likelihood (NegBin): a2
# 
# Example:
# model.name <- Richards_HSSM_model_definition(K = 1, 
#                                              S1 = "time", 
#                                              S2 = "time, 
#                                              P = P, 
#                                              lkhd = "NegBin", 
#                                              name = "M5a2")
#
# It returns the model file name

Richards_HSSM_model_definition <- function(K = 1, 
                                           S1 = "time", 
                                           S2 = "time", 
                                           P = "time", 
                                           Max = "time",
                                           lkhd = "Poisson", 
                                           name){
  filename <- paste0("models/model_Richards_HSSM_", name, ".jags")
  file.id <- file(filename, open = "wt")
  
  S1 <- ifelse(tolower(S1) == "time",  "-S1[y]", S1)
  S2 <- ifelse(tolower(S2) == "time", "S2[y]", S2)
  P <- ifelse(tolower(P) == "time", "P[y]", P)
  K <- ifelse(tolower(K) == "time", "K[y]", K)
  Max <- ifelse(tolower(Max) == "time", "Max[y]", Max)
  
  if (length(grep("M5", name)) > 0){
    S2.txt <- paste("for (y in 1:n.year){", "\n",
		                "      S2[y] ~ dgamma(S2.alpha, S2.beta)", "\n",  
	                  "   } #y", "\n")
    S1.txt <- paste("for (y in 1:n.year){", "\n",
                    "   S1[y] ~ dgamma(S1.alpha, S1.beta)", "\n",  
                    "   } #y", "\n")
  } 
  
  if (length(grep("M7", name)) > 0){
    S1.txt <- paste("for (y in 1:n.year){", "\n",
                    "      S1[y] ~ dgamma(S1.alpha, S1.beta)", "\n",
                    "   } #y", "\n")
    S2.txt <- paste0("S2 ~ dgamma(S2.alpha, S2.beta)", "\n")
  } 
  if (length(grep("M8", name)) > 0) {
    S2.txt <- paste("for (y in 1:n.year){", "\n",
                    "      S2[y] ~ dgamma(S2.alpha, S2.beta)", "\n",
                    "   } #y", "\n")
    S1.txt <- paste0("S1 ~ dgamma(S1.alpha, S1.beta)", "\n")
  }
  
  if (length(grep("M6", name)) > 0){
    S1.txt <- paste0("S1 ~ dgamma(S1.alpha, S1.beta)", "\n")
    S2.txt <- paste0("S2 ~ dgamma(S2.alpha, S2.beta)", "\n")
  }
  
  M1 <- paste("M1[t, y] <- (1 + (2 * exp(", K, ") - 1) * exp((1/(", S1, ")) * (", P, "- t))) ^ (-1/exp(", K, "))")
  M2 <- paste("M2[t, y] <- (1 + (2 * exp(", K, ") - 1) * exp((1/(", S2, ")) * (", P, "- t))) ^ (-1/exp(", K, "))")
  mean.N <- paste("mean.N[t, y] <- ", Max, " * (M1[t, y] * M2[t, y])")
  
  if (lkhd == "Poisson"){
    lkhd.txt <- paste("# Poisson Rate", "\n",
                      "            lambda[d, s, y] <- whales.available[d,s,y] * obs.prob[d,s,y]", "\n",
				              "            n[d, s, y] ~ dpois(lambda[d,s,y])", "\n")
    log.lkhd.txt <- paste("         log.lkhd[(d-1),s,y] <- logdensity.pois(n[d,s,y], lambda[d,s,y])", "\n")
  } else if (lkhd == "NegBin"){
    lkhd.txt <- paste("# 1. Calculate expected mean (kappa)", "\n", 
                      "            kappa[d, s, y] <- whales.available[d,s,y] * obs.prob[d, s, y]", "\n",
                      "            # 2. Convert kappa to the JAGS Negative Binomial probability parameter", "\n",
                      "            p_nb[d, s, y] <- r / (r + kappa[d, s, y])", "\n",
                      "            # 3. The new Negative Binomial likelihood function", "\n",
                      "            n[d, s, y] ~ dnegbin(p_nb[d, s, y], r)", "\n")
    log.lkhd.txt <- paste("log.lkhd[(d-1),s,y] <- logdensity.negbin(n[d,s,y], p_nb[d, s, y], r)", "\n")
  }
  
  model.text <- cat("model{", "\n", 
  "   for (y in 1:n.year){", "\n",    
  "      # Fix Day 1", "\n",
  "      # Use small epsilon to avoid log(0) errors if used in log-likelihood", "\n",
  "      mean.N[1, y] <- 0.0001", "\n", 
  "      # Fix Last Day", "\n",
  "      mean.N[n.days, y] <- 0.0001", "\n",
  "      for (t in 2:(n.days-1)){", "\n",
  "         ", M1, "\n",
  "         ", M2, "\n",
  "         ", mean.N, "\n", 
  "      } #t", "\n",
  "   } #y", "\n",
  "   \n",
  
  "   for (y in 1:n.year){", "\n",
  "      for (s in 1:n.station[y]){", "\n",
  "         for (d in 2:(periods[y,s]+1)){", "\n",
	"            # Calculate 'Available Whales' based on effort (watch length)", "\n",
	"            # watch.length must be in days", "\n",
  "            whales.available[d,s,y] <- mean.N[day[d,s,y], y] * watch.length[(d-1),s,y]", "\n",
  "            ", lkhd.txt, "\n",
  "            ", log.lkhd.txt, "\n",
	"            # probability is mean + covariate effects and watch.length as an offset", "\n",
	"            # bf and vs d index is -1 because d = 1 is for no whales. But, bf and vs", "\n",
  "            # don't have the records for d = 1.", "\n",
  "            # bf.1 and vs.1 indicate the raw Beaufort and visibility code, rather than ", "\n",
  "            # standardized (bf and vs).", "\n",
  "            # vs.1 - 1 is to make the minimum VS = 0 in the linear model", "\n",
  "            logit(obs.prob[d,s,y]) <- alpha[obs.fixed[d,s,y]] +", "\n",
  "                                      (BF.Fixed * bf.1[(d-1),s,y]) +", "\n",
  "                                      (VS.Fixed * (vs.1[(d-1),s,y] - 1))", "\n",
  "         } #d", "\n",
  "      } #s", "\n",
  "   } #y", "\n", 
  "   \n",
  
  "   r ~ dunif(0, 50)  # Dispersion parameter for the Negative Binomial", "\n",
  "   \n",
  
  "   # Random Effect Variance", "\n",
  "   # We restrict the standard deviation to prevent observers from drifting too far apart.", "\n",
  "   sd.obs ~ dunif(0, 1.5)  ", "\n",
  "   tau.obs <- pow(sd.obs, -2)", "\n",
  "   \n",
  
  "   # By pulling them from dnorm(0, tau.obs), they naturally shrink toward zero,", "\n",
  "   for (o in 1:n.obs.fixed){ ", "\n",
  "      alpha[o] ~ dnorm(0, tau.obs)", "\n", 
  "   }", "\n", 
  "   \n",
  
	"   # priors for Richards function parameters. ", "\n",
	"   # --- Priors for the Peak (P) Trend Model ---", "\n",
	"   beta0.P ~ dunif(30, 60)       # Global Intercept (Average peak day across all years)", "\n",
	"   beta1.P ~ dnorm(0, 0.25)      # Global Slope (How many days the peak shifts per year)", "\n",
  "   \n",
  
	"   sd.proc.P ~ dnorm(0, 1)T(0,)  # Standard deviation of the random year effects", "\n",
	"   tau.proc.P <- pow(sd.proc.P, -2)", "\n",
	
	"   # --- Apply the trend and random effect to each year ---", "\n",
	"   for (y in 1:n.year){", "\n",
	"      # 1. Calculate the expected mean peak for this specific year based on the trend", "\n",
	"      mu.P[y] <- beta0.P + (beta1.P * year.index[y])", "\n",
	"      # 2. Estimate the actual year's peak, allowing it to randomly deviate from the trend line", "\n",
  "		   P[y] ~ dnorm(mu.P[y], tau.proc.P)", "\n",
  "	  } #y", "\n",
  "   \n",
  
	"   # --- Priors for the Abundance (Max) Trend Model ---", "\n",
	"   beta0.Max ~ dnorm(7.6, 0.25)  # Global Intercept for log(Max)", "\n",
	"   beta1.Max ~ dnorm(0, 0.25)    # Global Slope (Rate of increase/decrease in Max per year)", "\n",
  "   \n",
  
	"   sd.proc.Max ~ dnorm(0, 1)T(0,)", "\n",
	"   tau.proc.Max <- pow(sd.proc.Max, -2)", "\n",
  "   \n",
  
	"   # --- Apply the trend and random effect ---", "\n",
	"   for (y in 1:n.year){", "\n",
	"	     mu.log.Max[y] <- beta0.Max + (beta1.Max * year.index[y])", "\n",
	"	     log.Max[y] ~ dnorm(mu.log.Max[y], tau.proc.Max)", "\n",
	"	     Max[y] <- exp(log.Max[y])", "\n",
	"   } #y", "\n",
		
	"   # Uniform distributions were adjusted to capture entire posterior ", "\n",
	"   # distributions of the hyper-parameters. ", "\n",
	"   # Hyper parameters for S1 and S2 were changed following suggestions by Gemini 2026-02-03", "\n",
	"   S1.alpha ~ dnorm(10, 0.1)T(0,)  #dunif(0.1, 50)", "\n",
	"   S1.beta ~ dgamma(1, 1)   #dunif(0.01, 5)", "\n",
	
	"   S2.alpha ~ dnorm(10, 0.1)T(0,)  #dunif(0.1, 50)", "\n",
	"   S2.beta ~ dgamma(1, 1)     #dunif(0.01, 10)", "\n",
  "  ", S1.txt, "\n",
  "  ", S2.txt, "\n",
  "   \n",
  
	"   ## Beaufort and visibility	", "\n",
	"   ## Gemini: The precision parameters were changed from 0.001 to 0.25 2026-02-03", "\n",
	"   BF.Fixed ~ dnorm(0,0.25) ", "\n",
	"   VS.Fixed ~ dnorm(0,0.25) ", "\n",
  "   \n",
  
	"   ### Summaries, Abundance Estimates, and Other Derived Quantities ", "\n",
	"   for(t in 1:n.year){", "\n",
	"	     Raw.Est[t] <- sum(mean.N[1:n.days, t])", "\n",
	"	     # multiply raw estimates by correction factor for nighttime passage rates (below)", "\n",
	"	     Corrected.Est[t] <- Raw.Est[t] * corr.factor ", "\n",
  "   }#t", "\n",
  "   \n",
  
	"   # Correction factor for nighttime passage rates (Perryman et al. 1999 Marine Mammal Science):", "\n",
	"   corr.factor ~ dnorm(mean.corr, tau.corr)", "\n",
	"   mean.corr <- 1.0875", "\n",
	"   sd.corr <- 0.03625", "\n",
	"   tau.corr <- pow(sd.corr,-2)", "\n",
  "} #model", "\n", file = file.id)
  
  close(file.id)
  #sink(file.id)
  return(filename)
  
}
