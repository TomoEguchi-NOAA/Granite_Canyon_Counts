n.year <- 8
time.com <- time.sp <- covariate.com <- covariate.sp <- matrix(nrow = 90, ncol = n.year)

for(d in 1:90){ #the full migration is considered to be 90 days starting on Dec 1
  for(t in 1:n.year){
    time.com[d,t] <- d #weird, just creating a vector of 1:90
    time.sp[d,t] <- d
  }
}

mean.time.com <- sd.time.com <- mean.time.sp <- sd.time.sp <- vector(mode = "numeric", length = n.year)
for(t in 1:n.year){
  mean.time.com[t]<-mean(time.com[1:90,t]) # this will always be 45.5 unless the 1:90 time frame changes
  sd.time.com[t]<-sd(time.com[1:90,t]) # and this will always be 26.1247 unless the time frame changes
  mean.time.sp[t]<-mean(time.sp[1:90,t])
  sd.time.sp[t]<-sd(time.sp[1:90,t])
}#t

mean.beta <- beta.sigma <- beta.tau <-  vector(mode = "numeric", length = 3)
beta <- matrix(nrow = 3, ncol = n.year)
for (l in 1:3){
  mean.beta[l] <- rnorm(1, 0,100) # means are specified as N(0,10)
  beta.sigma[l] <- runif(1, 0,10) # sd's are specified as U(0,10)
  beta.tau[l] <- beta.sigma[l] ^(-2)
  
  # Annual 'beta' perameters distributed around hyper-parameters
  for(t in 1:n.year){ #different betas for each year
    beta[l,t] <- rnorm(1, mean.beta[l], beta.sigma[l])
  }#t
}#l


com <- sp <- Common <- matrix(nrow = 90, ncol = n.year)
X <- array(dim = c(90, 3, n.year))
for(d in 1:90){ #the full migration is considered to be 90 days starting on Dec 1
  for(t in 1:n.year){

    covariate.com[d,t]<-(time.com[d,t] - mean.time.com[t])/sd.time.com[t] # this makes a 'covariate', which is just a straight line
    
    covariate.sp[d,t]<-(time.sp[d,t]-mean.time.sp[t])/sd.time.sp[t] # this makes a 'covariate', which is just a straight line
    
    for(l in 1:3){
      # this makes three vectors related to day 'd' - one a flat line, 
      # one a straight increasing line, one a curve - which are then 
      # multiplied by a, b, and c in the Common model
      X[d,l,t] <- covariate.com[d,t] ^ (l-1) 
    }#l
  }#t
  
  for(t in 1:n.year){ # Put all of the above together to calculate the Common model estimate for each day
    
    # X is the same across all years, and has the different shapes for 
    # each part of the polynomial (a * 1, b*d, c*d^2). Then it's multiplied 
    # by the different a,b,c from each year.
    
    com[d,t] <- beta[,t] %*% X[d,,t] #inprod(beta[,t], X[d,,t]) 
    # inprod does the calculation for each day / each year across all of the betas 
    
    # NOTE: Model is inverse-logged (exp) as below for the full season summation
    Common[d,t] <- exp(com[d,t]) # Both models are on the log scale, and then added in log space to the effort correction
    
  }#t
}#d
