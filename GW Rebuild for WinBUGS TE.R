
rm(list=ls())
#library(R2jags)
library(abind)
library(R2WinBUGS)

#Number of watch periods in each year's survey
periods <-c(136, 135, 164, 178, 179, 151)

#Watch start times, as fraction of a day
begin <- as.matrix(read.table("Formatted Data Files for 2019 Analysis/Data/begin.txt", 
                              header=T))

#watch end times
end <- as.matrix(read.table("Formatted Data Files for 2019 Analysis/Data/end.txt", header=T))

#whale counts
n <- as.matrix(read.table("Formatted Data Files for 2019 Analysis/Data/n.txt",header=T))
dim(n) <- c(179,2,6) #convert this back to a 3D array
#n <- abind(n,array(0,dim=c(2,2,6)), along=1) #add two trailing 0s to the end of the sightings array (this is for the day 1 and day 90 zero-whale anchor points)


n1 <- as.matrix(read.table("Formatted Data Files for 2019 Analysis/Data/n1.txt",header=T))      #These aren't needed; used in the WinBUGS GUI version previously. Here we can just re-use n
dim(n1) <- c(179,2,6) #convert this back to a 3D array


n2 <- as.matrix(read.table("Formatted Data Files for 2019 Analysis/Data/n2.txt",header=T))      #These aren't needed
dim(n2) <- c(179,2,6) #convert this back to a 3D array


u <- as.matrix(read.table("Formatted Data Files for 2019 Analysis/Data/u.txt",header=T))
dim(u) <- c(179,2,6) #convert this back to a 3D array
#u <- abind(u,array(0,dim=c(2,2,6)), along=1) #add two trailing 0s to the end of the effort on/off array
#for(i in 1:length(periods)){ #Place 1's for 'effort on' for the two periods following the end of true watches (this is for the day 1 and day 90 zero-whale anchor points)
#  u[(periods[i]+1):(periods[i]+2),,i] <- 1
#}


#visibility
vs <- as.matrix(read.table("Formatted Data Files for 2019 Analysis/Data/vs.txt",header=T))
#vs <- rbind(vs,matrix(NA,nrow=2,ncol=6)) #Add two trailing NAs (this is for the day 1 and day 90 zero-whale anchor points)

#beaufort
bf <- as.matrix(read.table("Formatted Data Files for 2019 Analysis/Data/bf.txt",header=T))
#bf <- rbind(bf,matrix(NA,nrow=2,ncol=6)) #Add two trailing NAs (this is for the day 1 and day 90 zero-whale anchor points)

#observer numbers
obs <- as.matrix(read.table("Formatted Data Files for 2019 Analysis/Data/obs.txt",header=T))
dim(obs) <- c(179,2,6) #convert this back to a 3D array
#obs <- abind(obs,array(0,dim=c(2,2,6)), along=1) #add two trailing 0s to the end of the effort on/off array
#for(i in 1:length(periods)){ #Place 36s for 'no observer' for the two periods following the end of true watches (this is for the day 1 and day 90 zero-whale anchor points)
#  obs[(periods[i]+1):(periods[i]+2),,i] <- 36 #this will force it to the mean observation probability with no observer effect
#}

N_inits <- as.matrix(read.table("Formatted Data Files for 2019 Analysis/Initial Values/N_inits.txt",header=T))

for(i in 1:length(periods)){
  N_inits[(periods[i]+1):(periods[i]+2),i] <- NA #we're going to make N a partially observed data object with anchor points at day 1 and 90
}

#The 'data' has to be the inverse of the inits, 
# with NAs for all of the estimated Ns, and 0s for the days 1 and 90
N <- matrix(NA,nrow=181,ncol=6) 

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

#### Two separate models each fit to data, 
# and then a third replicate of data to select
# 'best fit' (same structure as Durban et al)

# Shorten data to first 4 years only to replicate 
# the analysis in Durban et al 2016:
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

Watch.Length <- rbind(end,matrix(NA,nrow=2,ncol=6)) - rbind(begin,matrix(NA,nrow=2,ncol=6))
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
                        knot=c(-1.46,-1.26,-1.02,-0.78,
                               -0.58,-0.34,-0.10,0.10,
                               0.34,0.57,0.78,1.02,1.26,1.46),
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

jags.data <- list(n=n,
                  n.com=n,
                  n.sp=n,
                  n.station = dim(n)[2],
                  n.year = dim(n)[3],
                  n.obs = max(obs),
                  periods = periods,
                  obs=obs,
                  #Watch.Length = 0.0625,
                  u=u,
                  vs=vs,
                  bf=bf,
                  #day=day,
                  day=t,
                  N=N,
                  N.com=N,
                  N.sp=N,
                  knot=c(-1.46,-1.26,-1.02,-0.78,
                         -0.58,-0.34,-0.10,0.10,
                         0.34,0.57,0.78,1.02,1.26,1.46),
                  n.knots=14,
                  #begin=begin,
                  #end=end,
                  Watch.Length=Watch.Length)

jags.inits <- function() list(mean.prob = 0.5,
                              BF.Fixed = 0,
                              VS.Fixed = 0,
                              mean.prob.sp = 0.5,
                              BF.Fixed.sp = 0,
                              VS.Fixed.sp = 0,
                              mean.prob.com = 0.5,
                              BF.Fixed.com = 0,
                              VS.Fixed.com = 0,
                              mean.beta = c(0,0,0), #mean.beta = c(5,0.14,-3.5),
                              beta.sigma = c(1,1,1),#beta.sigma = c(7,7,7),
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
                              N = N_inits,
                              N.com = N_inits,
                              N.sp = N_inits,
                              #z = matrix(1,nrow=90,ncol=6),
                              beta.sp = array(data=0,dim=c(2,6)),
                              sd.b.sp = c(1,1,1,1,1,1),
                              z=matrix(1,nrow=90,ncol=6))

#### To run 2006-2019 data, load .RDS object and name it jags.data, 
# then update intitial values:
# Files are 2006-2019_GC_Formatted_Data and 2006-2019_GC_N_inits

# WinBUGS gives errors when N inits are set to 0. 
# Try setting them to 1 instead (seems to work):
N_inits[which(N_inits==0,arr.ind = T)] <- 1

parameters <- c("lambda","OBS.RF","OBS.Switch",
                "BF.Switch","BF.Fixed","VS.Switch",
                "VS.Fixed","mean.prob","mean.prob.com",
                "mean.prob.sp","BF.Fixed.com",
                "BF.Fixed.sp","VS.Fixed.com",
                "VS.Fixed.sp",
                "Corrected.Est","Raw.Est","z",
                "com","sp","Daily.Est","mean.beta",
                "beta.sigma","beta","beta.sp","b.sp","sd.b.sp")

ni <- 100000
nt <- 80
nb <- 60000
nc <- 3

# ni <- 5000
# nb <- 1000
# nt <- 20



#Run time: 
Start_Time<-Sys.time()

GW_Nmix <- bugs(data = jags.data,
                inits = jags.inits,
                parameters = parameters,
                model.file="GW_Nmix_Orig.bugs",
                n.chains = nc,
                n.iter = ni, n.burnin = nb, n.thin = nt,
                debug=F,
                bugs.directory = "C:/Users/tomo.eguchi/Documents/WinBUGS14")

Run_Time <- Sys.time() - Start_Time
save.image("T:/Eguchi/R_out/GW BUGS 7yr 100k.RData")

# library(R2jags)
# 
# GW_Nmix <- jags(jags.data.short, 
#                 jags.inits.short, parameters, 
#                 "GW_Nmix_Orig_noCUT.bugs", 
#                 n.chains = nc, n.thin = nt, 
#                 n.iter = ni, n.burnin = nb, 
#                 working.directory = getwd())
# 
# 
# 

#attach.jags(GW_Nmix_short)

#apply(Corrected.Est,2,median)


par(mfrow=c(2,2))
year <- 4

#John's plots:

UCIs <- apply(Daily.Est[,,year],2,quantile,0.975)
LCIs <- apply(Daily.Est[,,year],2,quantile,0.025)
plot(apply(exp(sp[,,year]),2,median), 
     type='l',
     ylim=c(0,max(UCIs)+50), 
     xlab="Days since 1 December", 
     ylab = "Whales per day")
lines(apply(exp(com[,,year]),2,median),type='l',lty=2)
segments(x0=1:90,y0=LCIs,y1=UCIs)

#Com vs Sp
plot(apply(exp(sp[,,year]),2,quantile,0.975),type='l',lty=2)
lines(apply(exp(sp[,,year]),2,median),type='l')
lines(apply(exp(sp[,,year]),2,quantile,0.025),type='l',lty=2)

plot(apply(exp(com[,,year]),2,quantile,0.975),type='l',lty=2)
lines(apply(exp(com[,,year]),2,median),type='l')
lines(apply(exp(com[,,year]),2,quantile,0.025),type='l',lty=2)

plot(apply(Daily.Est[,,year],2,quantile,0.975),type='l',lty=2)
lines(apply(Daily.Est[,,year],2,median),type='l')
lines(apply(Daily.Est[,,year],2,quantile,0.025),type='l',lty=2)

median(apply(exp(sp[,,4]),1,sum))

# PARAMETER NAME TRANSFERS:

# JOSH     JOHN
#    Spline 
# b.sp     b1 
# beta.sp  beta1
# X.sp     X1

