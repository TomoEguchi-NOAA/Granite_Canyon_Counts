
#library(R2jags)
library(abind)
library(R2WinBUGS)

Data.dir <- "C:/Users/tomoe/OneDrive/Documents/R/Granite_Canyon_Counts/Formatted Data Files for 2019 Analysis/Data"
WinBUGS.dir <- paste0(Sys.getenv("HOME"), "/WinBUGS14")

#Number of watch periods in each year's survey
periods <-c(136, 135, 164, 178, 179, 151)

#Watch start times, as fraction of a day
begin <- as.matrix(read.table(paste0(Data.dir, "/begin.txt"),
                              header=T, 
                              nrows = max(periods)))
#watch end times
end <- as.matrix(read.table(paste0(Data.dir, "/end.txt"),
                            header=T, 
                            nrows = max(periods)))

#whale counts
n <- as.matrix(read.table(paste0(Data.dir, "/n.txt"),
                          header=T, 
                          nrows = max(periods)))
dim(n) <- c(179,2,6) #convert this back to a 3D array
#n <- abind(n,array(0,dim=c(2,2,6)), along=1) #add two trailing 0s to the end of the sightings array (this is for the day 1 and day 90 zero-whale anchor points)


n1 <- as.matrix(read.table(paste0(Data.dir, "/n1.txt"),
                           header=T, 
                           nrows = max(periods)))      #These aren't needed; used in the WinBUGS GUI version previously. Here we can just re-use n
dim(n1) <- c(179,2,6) #convert this back to a 3D array


n2 <- as.matrix(read.table(paste0(Data.dir, "/n2.txt"),
                           header=T, 
                           nrows = max(periods)))      #These aren't needed
dim(n2) <- c(179,2,6) #convert this back to a 3D array


u <- as.matrix(read.table(paste0(Data.dir, "/u.txt"),
                          header=T, 
                          nrows = max(periods)))
dim(u) <- c(179,2,6) #convert this back to a 3D array
#u <- abind(u,array(0,dim=c(2,2,6)), along=1) #add two trailing 0s to the end of the effort on/off array
#for(i in 1:length(periods)){ #Place 1's for 'effort on' for the two periods following the end of true watches (this is for the day 1 and day 90 zero-whale anchor points)
#  u[(periods[i]+1):(periods[i]+2),,i] <- 1
#}


#visibility
vs <- as.matrix(read.table(paste0(Data.dir, "/vs.txt"),
                           header=T, 
                           nrows = max(periods)))
#vs <- rbind(vs,matrix(NA,nrow=2,ncol=6)) #Add two trailing NAs (this is for the day 1 and day 90 zero-whale anchor points)

#beaufort
bf <- as.matrix(read.table(paste0(Data.dir, "/bf.txt"),
                           header=T, 
                           nrows = max(periods)))
#bf <- rbind(bf,matrix(NA,nrow=2,ncol=6)) #Add two trailing NAs (this is for the day 1 and day 90 zero-whale anchor points)

#observer numbers
obs <- as.matrix(read.table(paste0(Data.dir, "/obs.txt"),
                            header=T, 
                            nrows = max(periods)))
dim(obs) <- c(179,2,6) #convert this back to a 3D array
#obs <- abind(obs,array(0,dim=c(2,2,6)), along=1) #add two trailing 0s to the end of the effort on/off array
#for(i in 1:length(periods)){ #Place 36s for 'no observer' for the two periods following the end of true watches (this is for the day 1 and day 90 zero-whale anchor points)
#  obs[(periods[i]+1):(periods[i]+2),,i] <- 36 #this will force it to the mean observation probability with no observer effect
#}

N_inits <- as.matrix(read.table(paste0(Data.dir, "/N_inits.txt"),
                                header=T,
                                nrows = (max(periods)+2)))

for(i in 1:length(periods)){
  N_inits[(periods[i]+1):(periods[i]+2),i] <- NA #we're going to make N a partially observed data object with anchor points at day 1 and 90
}

N <- matrix(NA,nrow=181,ncol=6) #The 'data' has to be the inverse of the inits, with NAs for all of the estimated Ns, and 0s for the days 1 and 90
for(i in 1:length(periods)){
  N[(periods[i]+1):(periods[i]+2),i] <- 0 #True number of whales passing fixed at 0 for day 1 and 90
}


#N_inits<-N_inits[1:179,] # These removals are necessary if you don't add the anchor points to the DATA inputs above
#N_inits[137:138,1] <- NA
#N_inits[136:137,2] <- NA
#N_inits[165:166,3] <- NA
#N_inits[179,4] <- NA
#N_inits[152:153,6] <- NA
# ^N inits were set as 2x observed n. John added two trailing 0s to N in the model, and in the N inits (which had to be removed above). You should revisit this to decide if N will be fixed to 0 in the model, or if it's better to do that with faux observations of 0 animals


#specify survey days associated with each watch period
day <- floor(begin)
t <- round((begin+end)/2)
#Add a couple of extra rows of NAs to the end of the day index reference to match up with the fixed 0s in N (above), assigning them to days 1 and 90
day <- rbind(as.matrix(day),matrix(NA,nrow=2,ncol=6))
for(i in 1:length(periods)){ #Set the anchor points: days 1 and 90
  day[(periods[i]+1):(periods[i]+2),i] <- c(1,90)
}

t <- rbind(as.matrix(t),matrix(NA,nrow=2,ncol=6))
for(i in 1:length(periods)){ #Set the anchor points: days 1 and 90
  t[(periods[i]+1):(periods[i]+2),i] <- c(1,90)
}



#### Two separate models each fit to data, and then a third replicate of data to select 'best fit' (same structure as Durban et al)

sink("GW_Nmix_Orig.bugs")
cat("

    model{
      
  
  # count i
  # station (Trailer) s
  # year t
  
  
  ### Linking Lambda from N Mixture to Common Model
  
    for(t in 1:n.year){
      for(d in 1:90){ #this is sub-indexed because the total number of shifts is different among years, and the begin / end times are formatted to have NAs where there were no watches
        #day[j,t] <- floor(begin[j,t]) #this is the day that each watch / count is associated with (I ended up inputing this as data instead of begin/end)
        #log(lambda[d,t]) <- log(Watch.Length) + com[d,t] #this is sub-indexed by day. There are multiple watches per day,
                                                              #all of which are poisson distributed around a mean (lambda) whales per watch
                                                              #The com[day,year] common model is the total number of whales passing on a given day.
                                                              #By adding the log(Watch.Length), we're multiplying that by the fraction of a day represented by a watch
                                                              #So the end result is expected whales per watch (lambda) = expected whales per day (f - common model) * watch period
                                                              #In this case all watch periods are the same length, but this could be adjusted for different watch lengths (i.e. watch.length[j,t] <- end[j,t] - begin[j,t])
                                                              #Lastly, each watch is poisson distributed around the daily expected whales per watch
  
  
        # Fit the two separate models to replicate count data:
        
#        log(lambda.com[d,t]) <- log(Watch.Length) + com[d,t] #common model daily lambda fit
#        log(lambda.sp[d,t]) <- log(Watch.Length) + sp[d,t] #specific model daily lambda fit
#        log(lambda[d,t]) <- log(Watch.Length) + selected[d,t] # the 'best-fit' model, selected through the deily z switch
        
        com.tmp[d,t] <- com[d,t] # coded with an intermediate 'temp' vector, identical to John's code (but with dsum() instead of cut())
        com.cut[d,t] <- cut(com.tmp[d,t]) #to prevent feedback (supposedly analogous to cut() function in BUGS)
        
        sp.tmp[d,t] <- sp[d,t]
        sp.cut[d,t] <- cut(sp.tmp[d,t])
      
  
        #selected[d,t] <- (z[d,t]*com[d,t]) + (1-z[d,t])*sp[d,t] #where com[] is the common seasonal curve model, sp[] is the specific spline model, and z[d,t] is the dbern(0.5) daily switch to select between model types. Also consider a simple annual switch, z[t] instead of selecting daily

        z[d,t] ~ dbern(0.5) #uninformative prior for model switch

      }#d
    }#t
    
    
    # Below mirrors John's code, indexing lambda by periods instead of by days, with the reference vector to link it to the appropriate model day
    for(t in 1:n.year){
      for(j in 1:(periods[t]+2)){
      
        selected[j,t] <- (z[day[j,t],t]*com.cut[day[j,t],t]) + (1-z[day[j,t],t])*sp.cut[day[j,t],t] #where com[] is the common seasonal curve model, sp[] is the specific spline model, and z[d,t] is the dbern(0.5) daily switch to select between model types. Also consider a simple annual switch, z[t] instead of selecting daily

        
        log(lambda.com[j,t]) <- log(Watch.Length[j,t]) + com[day[j,t],t] #common model daily lambda fit
        log(lambda.sp[j,t]) <- log(Watch.Length[j,t]) + sp[day[j,t],t] #specific model daily lambda fit
        
        #log(lambda[j,t]) <- log(Watch.Length[j,t]) + selected[day[j,t],t] # the 'best-fit' model, selected through the deily z switch
        log(lambda[j,t]) <- log(Watch.Length[j,t]) + selected[j,t] # the 'best-fit' model, selected through the deily z switch

      }#j
    }#t
  
  ### N Mixture process
  
  for(t in 1:n.year){
    for(j in 1:periods[t]){ #NOT periods + 2 because those day 1 & 90 anchor points are fixed, so N shouldn't be estimated for those, only for the days with true watch periods
      for(s in 1:n.station){
      
        n[j,s,t] ~ dbin(obs.prob[j,s,t],N[j,t]) # N mixture model - two observations of 'true' N (passing whales), with different observation probabilities
        n.com[j,s,t] ~ dbin(obs.prob.com[j,s,t],N.com[j,t]) # N mixture model for common model fit
        n.sp[j,s,t] ~ dbin(obs.prob.sp[j,s,t],N.sp[j,t]) # N mixture model for specific model fit

      }#s    
    }#j
  }#t

  for(t in 1:n.year){
    for(j in 1:(periods[t]+2)){ #periods + 2 to include the added anchor points of 0 whales at days 1 and 90


      # Replicate John's code:
      N[j,t] ~ dpois(lambda[j,t]) # Then, the 'true' N's from a given day are distributed around a poisson mean whales per watch (lambda)
      N.com[j,t] ~ dpois(lambda.com[j,t]) # Then, the 'true' N's from a given day are distributed around a poisson mean whales per watch (lambda)
      N.sp[j,t] ~ dpois(lambda.sp[j,t]) # Then, the 'true' N's from a given day are distributed around a poisson mean whales per watch (lambda)

    }#j
  }#t
  
  
  ### Observation probability
  # It's necessary to have three separate observation probabilities for the three different models (common, spline, and selected)
  # The reason: the common model will force a 'normal' seasonal curve, even if the number of whales in the middle of the season drops substantially
  # If the spline model is then forced to share the same observation probability, it will overestimate the number of whales (because the seasonal curve would lead to a very low sighting probability during those days)
  for(t in 1:n.year){
    for(j in 1:periods[t]){
      for(s in 1:n.station){

        #Final obs prob    = observer 1/0 * (base obs prob + BF on/off *Fixed effect of BF * BF + VS on/off *Fixed effect of VS * VS    + Obs effect on/off * Random effect of observer)

        #The below replicates Durban et al 2016 code:
        # Selected model:
        logit(prob[j,s,t]) <- logit(mean.prob) + (BF.Switch*BF.Fixed*bf[j,t]) + (VS.Switch*VS.Fixed*vs[j,t]) + (OBS.Switch*OBS.RF[obs[j,s,t]])
        # the u data is whether there were observers on watch. 0 counts are often associated with years/shifts with no second observer. So if u=0, it will fix observation probability at 0
        obs.prob[j,s,t] <- u[j,s,t]*prob[j,s,t]
        
        # Common (com) model:
        logit(prob.com[j,s,t]) <- logit(mean.prob.com) + (BF.Switch.com*BF.Fixed.com*bf[j,t]) + (VS.Switch.com*VS.Fixed.com*vs[j,t]) + (OBS.Switch.com*OBS.RF.com[obs[j,s,t]])
        obs.prob.com[j,s,t] <- u[j,s,t]*prob.com[j,s,t]
        
        # Spline (sp) model:
        logit(prob.sp[j,s,t]) <- logit(mean.prob.sp) + (BF.Switch.sp*BF.Fixed.sp*bf[j,t]) + (VS.Switch.sp*VS.Fixed.sp*vs[j,t]) + (OBS.Switch.sp*OBS.RF.sp[obs[j,s,t]])
        obs.prob.sp[j,s,t] <- u[j,s,t]*prob.sp[j,s,t]
        
        
      }#s
    }#j
  }#t
  
  #Uninformative prior for mean.prob
  mean.prob ~ dunif(0,1)
    mean.prob.com ~ dunif(0,1)
    mean.prob.sp ~ dunif(0,1)
  
  
  ### Specification of terms within observation probability linear model
  
  ## Observer random effect
  # SLECTED MODEL

   for(o in 1:n.obs){
    OBS.RF[o] ~ dnorm(0,tau.Obs)
  }#o
  
  #Uninformative prior for tau.Obs
  sigma.Obs ~ dunif(0,2)
  tau.Obs <- pow(sigma.Obs,-2)
  
  # COMMON MODEL

   for(o in 1:n.obs){
    OBS.RF.com[o] ~ dnorm(0,tau.Obs.com)
  }#o
  
  #Uninformative prior for tau.Obs
  sigma.Obs.com ~ dunif(0,2)
  tau.Obs.com <- pow(sigma.Obs.com,-2)
  
  # SPECIFIC MODEL

   for(o in 1:n.obs){
    OBS.RF.sp[o] ~ dnorm(0,tau.Obs.sp)
  }#o
  
  #Uninformative prior for tau.Obs
  sigma.Obs.sp ~ dunif(0,2)
  tau.Obs.sp <- pow(sigma.Obs.sp,-2)
  
  OBS.Switch ~ dbern(0.5)
    OBS.Switch.com ~ dbern(0.5)
    OBS.Switch.sp ~ dbern(0.5)

  
  ## Beaufort
  BF.Switch ~ dbern(0.5) #uninformative prior for the BF.Switch, which determines whether to include the effect of beaufort conditions (multiply by 0 or 1)
    BF.Switch.com ~ dbern(0.5)
    BF.Switch.sp ~ dbern(0.5)
  
  #Below is the single fixed effect multiplied by BF rating in the updated obs prob equation:
  BF.Fixed ~ dnorm(0,0.01)
    BF.Fixed.com ~ dnorm(0,0.01)
    BF.Fixed.sp ~ dnorm(0,0.01)
  
  ## Visibility
  VS.Switch ~ dbern(0.5) #uninformative prior for the VS.Switch
    VS.Switch.com ~ dbern(0.5)
    VS.Switch.sp ~ dbern(0.5)
    
  #Below is the single fixed effect multiplied by VS rating in the updated obs prob equation:
  VS.Fixed ~ dnorm(0,0.01)
    VS.Fixed.com ~ dnorm(0,0.01)
    VS.Fixed.sp ~ dnorm(0,0.01)
  
  
  
  ### Seasonal Curve Models
  
  ## Model 1, normal curve Common Model with shared hyper-parameters
  
  # Hyper parameter prior specification
  # These correspond to 'beta' parameters a, b, and c in the 'Common' model, Durban et al 2015
  for (l in 1:3){
    mean.beta[l] ~ dnorm(0,0.01) # means are specified as N(0,10)
    beta.sigma[l] ~ dunif(0,10) # sd's are specified as U(0,10)
    beta.tau[l] <- pow(beta.sigma[l],-2)
    
    # Annual 'beta' perameters distributed around hyper-parameters
    for(t in 1:n.year){ #different betas for each year
      beta[l,t] ~ dnorm(mean.beta[l],beta.tau[l])
    }#t
  }#l
  
  # Calculate the seasonal migration curve effect for the 'Common' model
  # This is some wild math/code from John Durban's original code, but it works:
  
  # mean and sd of the time vector below
      for(t in 1:n.year){
        mean.time.com[t]<-mean(time.com[1:90,t]) # this will always be 45.5 unless the 1:90 time frame changes
        sd.time.com[t]<-sd(time.com[1:90,t]) # and this will always be 26.1247 unless the time frame changes
        mean.time.sp[t]<-mean(time.sp[1:90,t])
        sd.time.sp[t]<-sd(time.sp[1:90,t])
      }#t
  
  
  for(d in 1:90){ #the full migration is considered to be 90 days starting on Dec 1
    for(t in 1:n.year){
      time.com[d,t] <- d #weird, just creating a vector of 1:90
      time.sp[d,t] <- d
      covariate.com[d,t]<-(time.com[d,t]-mean.time.com[t])/sd.time.com[t] # this makes a 'covariate', which is just a straight line
      covariate.sp[d,t]<-(time.sp[d,t]-mean.time.sp[t])/sd.time.sp[t] # this makes a 'covariate', which is just a straight line


      for(l in 1:3){
        X[d,l,t] <- pow(covariate.com[d,t],l-1) # this makes three vectors related to day 'd' - one a flat line, one a straight increasing line, one a curve - which are then multiplied by a, b, and c in the Common model
      }#l
    }#t
    for(t in 1:n.year){ # Put all of the above together to calculate the Common model estimate for each day
      
      com[d,t]<-inprod(beta[,t],X[d,,t]) # X is the same across all years, and has the different shapes for each part of the polynomial (a * 1, b*d, c*d^2). Then it's multiplied by the different a,b,c from each year. 
      # inprod does the calculation for each day / each year across all of the betas 
      
      # NOTE: Model is inverse-logged (exp) as below for the full season summation
      log(Common[d,t]) <- com[d,t] # Both models are on the log scale, and then added in log space to the effort correction
      
    }#t
  }#d
  
  
  ## Model 2, spline fit Specific Model
  for(t in 1:n.year){  
    for(d in 1:90){
      for(k in 1:n.knots){
        Z1[d,k,t]<-pow(uZ1[d,k,t], 1)
        uZ1[d,k,t]<-(covariate.sp[d,t]-knot[k])*step(covariate.sp[d,t]-knot[k])
      }#k
    }#d
  }#t
  
    for(d in 1:90){
      for (l in 1:2){
        for(t in 1:n.year){
          X.sp[d,l,t] <- pow(covariate.sp[d,t],l-1) # this makes three vectors related to day 'd' - one a flat line, one a straight increasing line (similar to the common model above but a second degree polynomial instead of a third degree poly)
        }#t
      }#l
    }#d
    
  for(t in 1:n.year){  
    for(k in 1:n.knots){
      b.sp[k,t]~dnorm(0,tau.b.sp[t]) #annual regression coefficients for each spline knot
    }#k
    
    tau.b.sp[t]<-pow(sd.b.sp[t],-2)
    sd.b.sp[t]~dunif(0,10) #uniform prior on regression coefficient SD, as per Durban et al
  
    for(l in 1:2){
      beta.sp[l,t]~dnorm(0,0.01) #N(0,10) prior for S0 and S1 coefficients, as per Durban et al
    }#l

  
    for(d in 1:90){
      sp[d,t] <- inprod(beta.sp[,t],X.sp[d,,t]) + inprod(b.sp[,t],Z1[d,,t]) # multiplying splines across days to make the penalized spline model fit

      # NOTE: Model is inverse-logged (exp) as below for the full season summation
      log(Specific[d,t]) <- sp[d,t] # Both models are on the log scale, and then added in log space to the effort correction
    }#d
  }#t
  
  ### Summaries, Abundance Estimates, and Other Derived Quantities
  
  ## Seasonal Abundance Estimate:
  for(t in 1:n.year){
    for(d in 1:90){
      # Daily estimate, based on either the common model (f) or the specific model (sp). For days with observations, this should have a somewhat confident z estimate. For days with no watches, I assume it will be balanced 50/50 between the two models
      log(Daily.Est[d,t]) <- z[d,t]*com.cut[d,t] + (1-z[d,t])*sp.cut[d,t] #where sp[] is the specific model, and z[j,t] is the dbern(0.5) daily switch to select between model types. Also consider a simple annual switch, z[t] instead of selecting daily

      #Daily.Est[d,t] <- z[d,t]*Common[d,t] + (1-z[d,t])*Specific[d,t]

    }#d
    raw.unrounded[t] <- sum(Daily.Est[1:90,t])
    Raw.Est[t] <- round(raw.unrounded[t])
    #Raw.Est[t] <- sum(Daily.Est[1:90,t])
    Corrected.Est[t] <- Raw.Est[t]*corr.factor # multiply raw estimates by correction factor for nighttime passage rates (below)
  }#t
  
  # Correction factor for nighttime passage rates:
  corr.factor~dnorm(mean.corr,tau.corr)
  mean.corr<-1.0875
  sd.corr<-0.03625
  tau.corr<-pow(sd.corr,-2)
  
  
  
    }#model
    ",fill = TRUE)
sink()




#Shorten data to first 4 years only to replicate the analysis in Durban et al 2016:
n.short <- n[,,1:4]
obs.short <- obs[,,1:4]
periods.short <- periods[1:4]
u.short <- u[,,1:4]
vs.short <- vs[,1:4]
bf.short <- bf[,1:4]
day.short <-day[,1:4]
t.short <- t[,1:4]
N.short <- N[,1:4]
N_inits.short <- N_inits[,1:4]
N_inits.short[which(N_inits.short == 0,arr.ind = T)] <- 1

Watch.Length <- rbind(end,matrix(NA,nrow=2,ncol=6)) - 
  rbind(begin,matrix(NA,nrow=2,ncol=6))
for(i in 1:length(periods)){ #Place 36s for 'no observer' for the two periods following the end of true watches (this is for the day 1 and day 90 zero-whale anchor points)
  Watch.Length[(periods[i]+1):(periods[i]+2),i] <- 1 #this will force it to the mean observation probability with no observer effect
}
Watch.Length.short <- Watch.Length[,1:4]


jags.data.short <- list(n=n.short,
                        n.com=n.short,
                        n.sp=n.short,
                        n.station = dim(n.short)[2],
                        n.year = dim(n.short)[3],
                        n.obs = max(obs.short),
                        periods = periods.short,
                        obs=obs.short,
                        #Watch.Length = 0.0625,
                        u=u.short,
                        vs=vs.short,
                        bf=bf.short,
                        day=t.short,
                        #day=day.short,
                        N=N.short,
                        N.com=N.short,
                        N.sp=N.short,
                        knot=c(-1.46,-1.26,-1.02,-0.78,-0.58,
                               -0.34,-0.10,0.10,0.34,0.57,
                               0.78,1.02,1.26,1.46),
                        n.knots=14,
                        Watch.Length=Watch.Length.short)

jags.inits.short <- function() list(mean.prob = 0.5,
                                    BF.Fixed = 0,
                                    VS.Fixed = 0,
                                    mean.prob.sp = 0.5,
                                    BF.Fixed.sp = 0,
                                    VS.Fixed.sp = 0,
                                    mean.prob.com = 0.5,
                                    BF.Fixed.com = 0,
                                    VS.Fixed.com = 0,
                                    mean.beta = c(0,0,0), #mean.beta = c(5,0.14,-3.5),
                                    beta.sigma = c(0.5,0.5,0.5),#beta.sigma = c(7,7,7),
                                    BF.Switch = 1,
                                    VS.Switch = 1,
                                    OBS.Switch = 1,
                                    sigma.Obs = 1,
                                    BF.Switch.sp = 1,
                                    VS.Switch.sp = 1,
                                    OBS.Switch.sp = 1,
                                    sigma.Obs.sp = 1,
                                    BF.Switch.com = 1,
                                    VS.Switch.com = 1,
                                    OBS.Switch.com = 1,
                                    sigma.Obs.com = 1,
                                    N = N_inits.short,
                                    N.com = N_inits.short,
                                    N.sp = N_inits.short,
                                    #z = matrix(1,nrow=90,ncol=6),
                                    beta.sp = array(data=0,dim=c(2,4)),
                                    sd.b.sp = c(1,1,1,1),
                                    z=matrix(1,nrow=90,ncol=4))

# jags.data <- list(n=n,
#                   n.com=n,
#                   n.sp=n,
#                   n.station = dim(n)[2],
#                   n.year = dim(n)[3],
#                   n.obs = max(obs),
#                   periods = periods,
#                   obs=obs,
#                   #Watch.Length = 0.0625,
#                   u=u,
#                   vs=vs,
#                   bf=bf,
#                   #day=day,
#                   day=t,
#                   N=N,
#                   N.com=N,
#                   N.sp=N,
#                   knot=c(-1.46,-1.26,-1.02,-0.78,-0.58,-0.34,-0.10,0.10,0.34,0.57,0.78,1.02,1.26,1.46),
#                   n.knots=14,
#                   #begin=begin,
#                   #end=end,
#                   Watch.Length=Watch.Length)
# 
# jags.inits <- function() list(mean.prob = 0.5,
#                               BF.Fixed = 0,
#                               VS.Fixed = 0,
#                               mean.prob.sp = 0.5,
#                               BF.Fixed.sp = 0,
#                               VS.Fixed.sp = 0,
#                               mean.prob.com = 0.5,
#                               BF.Fixed.com = 0,
#                               VS.Fixed.com = 0,
#                               mean.beta = c(0,0,0), #mean.beta = c(5,0.14,-3.5),
#                               beta.sigma = c(1,1,1),#beta.sigma = c(7,7,7),
#                               BF.Switch = 1,
#                               VS.Switch = 1,
#                               OBS.Switch = 1,
#                               sigma.Obs = 1,
#                               BF.Switch.sp = 1,
#                               VS.Switch.sp = 1,
#                               OBS.Switch.sp = 1,
#                               sigma.Obs.sp = 1,
#                               BF.Switch.com = 1,
#                               VS.Switch.com = 1,
#                               OBS.Switch.com = 1,
#                               sigma.Obs.com = 1,
#                               N = N_inits,
#                               N.com = N_inits,
#                               N.sp = N_inits,
#                               #z = matrix(1,nrow=90,ncol=6),
#                               beta.sp = array(data=0,dim=c(2,6)),
#                               sd.b.sp = c(1,1,1,1,1,1),
#                               z=matrix(1,nrow=90,ncol=6))



#### To run 2006-2019 data, load .RDS object and name it jags.data, then update intitial values:
# Files are 2006-2019_GC_Formatted_Data and 2006-2019_GC_N_inits


# WinBUGS gives errors when N inits are set to 0. 
#Try setting them to 1 instead (seems to work):
#N_inits[which(N_inits==0,arr.ind = T)] <- 1

# jags.inits <- function() list(mean.prob = 0.5,
#                               BF.Fixed = 0,
#                               VS.Fixed = 0,
#                               mean.prob.sp = 0.5,
#                               BF.Fixed.sp = 0,
#                               VS.Fixed.sp = 0,
#                               mean.prob.com = 0.5,
#                               BF.Fixed.com = 0,
#                               VS.Fixed.com = 0,
#                               mean.beta = c(0,0,0), #mean.beta = c(5,0.14,-3.5),
#                               beta.sigma = c(1,1,1),#beta.sigma = c(7,7,7),
#                               BF.Switch = 1,
#                               VS.Switch = 1,
#                               OBS.Switch = 1,
#                               sigma.Obs = 1,
#                               BF.Switch.sp = 1,
#                               VS.Switch.sp = 1,
#                               OBS.Switch.sp = 1,
#                               sigma.Obs.sp = 1,
#                               BF.Switch.com = 1,
#                               VS.Switch.com = 1,
#                               OBS.Switch.com = 1,
#                               sigma.Obs.com = 1,
#                               N = N_inits,
#                               N.com = N_inits,
#                               N.sp = N_inits,
#                               #z = matrix(1,nrow=90,ncol=6),
#                               beta.sp = array(data=0,dim=c(2,7)),
#                               sd.b.sp = c(1,1,1,1,1,1,1),
#                               z=matrix(1,nrow=90,ncol=7))



parameters <- c("lambda","OBS.RF","OBS.Switch",
                "BF.Switch","BF.Fixed","VS.Switch",
                "VS.Fixed","mean.prob","mean.prob.com",
                "mean.prob.sp","BF.Fixed.com","BF.Fixed.sp",
                "VS.Fixed.com","VS.Fixed.sp",
                "Corrected.Est","Raw.Est","z","com","sp",
                "Daily.Est","mean.beta","beta.sigma",
                "beta","beta.sp","b.sp","sd.b.sp")

ni <- 100000
nt <- 80
nb <- 60000
nc <- 3


#GW_Nmix <- jags(jags.data, inits=jags.inits, parameters, "GW_Nmix_Orig.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

#GW_Nmix_short <- jags(jags.data.short, inits=jags.inits.short, parameters, "GW_Nmix_Orig.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())


library(R2WinBUGS)

#Run time: 
Start_Time<-Sys.time()

GW_Nmix <- bugs(data = jags.data.short,
                inits = jags.inits.short,
                parameters = parameters,
                model.file="GW_Nmix_Orig.bugs",
                n.chains = nc,
                n.iter = ni, n.burnin = nb, n.thin = nt,
                debug=T,
                bugs.directory = WinBUGS.dir)
#"C:/Users/joshua.stewart/Desktop/Gray Whale Abundance Estimates/WinBUGS14/")

Run_Time <- Sys.time() - Start_Time
#save.image("GW BUGS 7yr 100k.RData")
save.image("RData/GW BUGS 4yr 100k v0.RData")

# library(R2jags)
# 
# GW_Nmix <- jags(jags.data.short, jags.inits.short, parameters, "GW_Nmix_Orig_noCUT.bugs", 
#                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
# 
# 
# 
# 
# attach.jags(GW_Nmix_short)
# 
# apply(Corrected.Est,2,median)
# 
# 
# par(mfrow=c(2,2))
# year <- 4
# 
# #John's plots:
# 
# UCIs <- apply(Daily.Est[,,year],2,quantile,0.975)
# LCIs <- apply(Daily.Est[,,year],2,quantile,0.025)
# plot(apply(exp(sp[,,year]),2,median),type='l',ylim=c(0,max(UCIs)+50), xlab="Days since 1 December", ylab = "Whales per day")
# lines(apply(exp(com[,,year]),2,median),type='l',lty=2)
# segments(x0=1:90,y0=LCIs,y1=UCIs)
# 
# 
# 
# #Com vs Sp
# plot(apply(exp(sp[,,year]),2,quantile,0.975),type='l',lty=2)
# lines(apply(exp(sp[,,year]),2,median),type='l')
# lines(apply(exp(sp[,,year]),2,quantile,0.025),type='l',lty=2)
# 
# plot(apply(exp(com[,,year]),2,quantile,0.975),type='l',lty=2)
# lines(apply(exp(com[,,year]),2,median),type='l')
# lines(apply(exp(com[,,year]),2,quantile,0.025),type='l',lty=2)
# 
# plot(apply(Daily.Est[,,year],2,quantile,0.975),type='l',lty=2)
# lines(apply(Daily.Est[,,year],2,median),type='l')
# lines(apply(Daily.Est[,,year],2,quantile,0.025),type='l',lty=2)
# 
# median(apply(exp(sp[,,4]),1,sum))
# 
# # PARAMETER NAME TRANSFERS:
# 
# # JOSH     JOHN
# #    Spline 
# # b.sp     b1 
# # beta.sp  beta1
# # X.sp     X1
# 
