#Richards_HSSM_model_definitions
# A function to create JAGS models for Richards_HSSM abundance estimation
# 
# This function creates a text file of a model with user input of which
# parameters of Richards function are time-specific. If any parameter is fixed 
# to a constant, provide the value. 
# 
# Richards function:
# 
# C1 <- (1 + (2 * exp(K) - 1) * exp((1/S1) * (P - d))) ^ (-1/exp(K))
# C2 <- (1 + (2 * exp(K) - 1) * exp((1/S2) * (P - d))) ^ (-1/exp(K))
# N <- min + (max - min) * (C1 * C2)
# 
# Provide either "time", a constant, or the parameter name to S1 and S2. 
# If "time", parameters will be specified as time-specific. For example, S1 = "time"
# will result in S1[y]. To fix the parameter, provide S1 = "S1". P and Max
# are assumed to be time/season specific and K = 1.
# 
# The current convention is:
# M1: P = "time", Max = "time", S1 = "time", S2 = "time", K = 1
# M2: P = "time", Max = "time", S1 = "S1", S2 = "S2", K = 1
# M3: P = "time", Max = "time", S1 = "time", S2 = "S2", K = 1
# M4: P = "time", Max = "time", S1 = "S1", S2 = "time", K = 1
# 
# Poisson likelihood (Poisson): a1
# Negative binomial likelihood (NegBin): a2
# 
# Example:
# model.name <- Richards_HSSM_model_definition(K = 1, 
#                                              S1 = "time", 
#                                              S2 = "time, 
#                                              P = "time, 
#                                              lkhd = "NegBin")
#
# It returns the model file name

Richards_HSSM_model_definition <- function(K = 1, 
                                           S1 = "time", 
                                           S2 = "time", 
                                           P = "time", 
                                           Max = "time",
                                           lkhd = "Poisson"){
  Run.Time <- as.character(Sys.time())
  
  if (P == "time" & Max == "time" & S1 == "time" & S2 == "time") M <- "M1"
  if (P == "time" & Max == "time" & S1 == "S1" & S2 == "S2") M <- "M2"
  if (P == "time" & Max == "time" & S1 == "time" & S2 == "S2") M <- "M3"
  if (P == "time" & Max == "time" & S1 == "S1" & S2 == "time") M <- "M4"
  
  if (tolower(lkhd) == "poisson") name <- paste0(M, "a1")
  if (tolower(lkhd) == "negbin") name <- paste0(M, "a2")
  
  filename <- paste0("models/model_Richards_HSSM_", name, ".jags")
  file.id <- file(filename, open = "wt")
  
  S1.par <- ifelse(tolower(S1) == "time",  "-S1[y]", paste0("-", S1))
  S2.par <- ifelse(tolower(S2) == "time", "S2[y]", S2)
  P.par <- ifelse(tolower(P) == "time", "P[y]", P)
  K.par <- ifelse(tolower(K) == "time", "K[y]", K)
  Max.par <- ifelse(tolower(Max) == "time", "Max[y]", Max)
  
  if (tolower(S2) == "time"){
    S2.txt <- paste("for (y in 1:n.year){", "\n",
		                "      S2[y] ~ dgamma(S2.alpha, S2.beta)", "\n",  
	                  "   } #y", "\n")
  } else {
    S2.txt <- paste0("S2 ~ dgamma(S2.alpha, S2.beta)", "\n")
  }
  
  if (tolower(S1) == "time"){
    S1.txt <- paste("for (y in 1:n.year){", "\n",
                    "   S1[y] ~ dgamma(S1.alpha, S1.beta)", "\n",  
                    "   } #y", "\n")
  } else {
    S1.txt <- paste0("S1 ~ dgamma(S1.alpha, S1.beta)", "\n")  
  }
 
  C1.txt <- paste("C1[t, y] <- (1 + (2 * exp(", K.par, ") - 1) * exp((1/(", S1.par, ")) * (", P.par, "- t))) ^ (-1/exp(", K.par, "))")
  C2.txt <- paste("C2[t, y] <- (1 + (2 * exp(", K.par, ") - 1) * exp((1/(", S2.par, ")) * (", P.par, "- t))) ^ (-1/exp(", K.par, "))")
  mean.N <- paste("mean.N[t, y] <- ", Max.par, " * (C1[t, y] * C2[t, y])")
  
  if (tolower(lkhd) == "poisson"){
    lkhd.txt <- paste("# Poisson Rate", "\n",
                      "            lambda[d, s, y] <- whales.available[d,s,y] * obs.prob[d,s,y]", "\n",
				              "            n[d, s, y] ~ dpois(lambda[d,s,y])", "\n")
    log.lkhd.txt <- paste("         log.lkhd[(d-1),s,y] <- logdensity.pois(n[d,s,y], lambda[d,s,y])", "\n")
  } else if (tolower(lkhd) == "negbin"){
    lkhd.txt <- paste("# 1. Calculate expected mean (kappa)", "\n", 
                      "            kappa[d, s, y] <- whales.available[d,s,y] * obs.prob[d, s, y]", "\n",
                      "            # 2. Convert kappa to the JAGS Negative Binomial probability parameter", "\n",
                      "            p_nb[d, s, y] <- r / (r + kappa[d, s, y])", "\n",
                      "            # 3. The new Negative Binomial likelihood function", "\n",
                      "            n[d, s, y] ~ dnegbin(p_nb[d, s, y], r)", "\n")
    log.lkhd.txt <- paste("log.lkhd[(d-1),s,y] <- logdensity.negbin(n[d,s,y], p_nb[d, s, y], r)", "\n")
  }
  
  model.text <- cat("# Created:", Run.Time, "\n",
                    "# Craeted by: Richards_HSSM_model_definition.R", "\n",
                    "# A JAGS model for Hierarchical State Space Modeling using", "\n",
                    "# the Richards function to estimate gray whale abundance from", "\n",
                    "# count data at Granite Canyon, CA. ", "\n\n", 
                    
                    "# The current model naming convention is:", "\n",
                    "# M1: P = 'time', Max = 'time', S1 = 'time', S2 = 'time, K = 1", "\n",
                    "# M2: P = 'time', Max = 'time', S1 = 'S1', S2 = 'S2', K = 1", "\n",
                    "# M3: P = 'time', Max = 'time', S1 = 'time', S2 = 'S2', K = 1", "\n",
                    "# M4: P = 'time', Max = 'time', S1 = 'S1', S2 = time, K = 1", "\n\n",
                    
                    "# Poisson likelihood (Poisson): a1", "\n",
                    "# Negative binomial likelihood (NegBin): a2", "\n",
                    "\n",
                    
  "model{", "\n", 
  "   for (y in 1:n.year){", "\n",    
  "      # Fix Day 1", "\n",
  "      # Use small epsilon to avoid log(0) errors if used in log-likelihood", "\n",
  "      mean.N[1, y] <- 0.0001", "\n", 
  "      # Fix Last Day", "\n",
  "      mean.N[n.days, y] <- 0.0001", "\n",
  "      for (t in 2:(n.days-1)){", "\n",
  "         ", C1.txt, "\n",
  "         ", C2.txt, "\n",
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
	"            # probability is observer effects + covariate effects", "\n",
	"            # bf and vs d index is -1 because when d = 1, there no whales (i.e., day 1).", "\n",
  "            # But, bf and vs don't have the records for d = 1.", "\n",
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
  
  "   # Dispersion parameter for the Negative Binomial", "\n", 
  "   r ~ dunif(0, 50)  ", "\n",
  "   \n",
  
  "   # Random Effect Variance", "\n",
  "   # We restrict the standard deviation to prevent observers from drifting too far apart.", "\n",
  "   sd.obs ~ dunif(0, 1.5)  ", "\n",
  "   tau.obs <- pow(sd.obs, -2)", "\n",
  "   \n",
  
  "   # By pulling them from dnorm(0, tau.obs), they naturally shrink toward zero,", "\n",
  "   for (o in 1:n.obs.fixed){ ", "\n",
  "      alpha[o] ~ dnorm(0, tau.obs)", "\n", 
  "   } #o", "\n", 
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
	"      mu.P[y] <- beta0.P + (beta1.P * year.index[y])", "\n\n",
  
	"      # 2. Estimate the actual year's peak, allowing it to randomly deviate from the trend line", "\n",
  "		   P[y] ~ dnorm(mu.P[y], tau.proc.P)", "\n",
  "	  } #y", "\n",
  "   \n",
  
	"   # --- Priors for the Max Trend Model ---", "\n",
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
  "  ", S1.txt, 
  "  ", S2.txt, 
  "   \n",
  
	"   # Beaufort and visibility	", "\n",
	"   # Gemini: The precision parameters were changed from 0.001 to 0.25 2026-02-03", "\n",
	"   BF.Fixed ~ dnorm(0,0.25) ", "\n",
	"   VS.Fixed ~ dnorm(0,0.25) ", "\n",
  "   \n",
  
	"   # Summaries, Abundance Estimates, and Other Derived Quantities ", "\n",
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
  
  return(filename)
  
}
