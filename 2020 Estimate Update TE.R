
library(R2jags)
library(abind)
library(dplyr)


# load 2019 formatted data
GC_2019 <- readRDS("RData/GC_2019_FinalData.rds")  # TE

# TE: This is poor programming. Hard coding these numbers are not good.
periods <-c(136, 135, 164, 178, 179, 151, 196)


ObsList <- read.csv("Observer list.csv")

match(unique(GC_2019$obs),ObsList$FirstOfPobs)
NewObs<-unique(GC_2019$obs)[which(is.na(match(unique(GC_2019$obs),ObsList$FirstOfPobs)))]

NewObsList <- rbind(as.matrix(ObsList),matrix(c(NewObs,42:50),ncol=2))

write.csv(NewObsList, "2019 Observer list.csv",row.names = F)

obs2019 <- as.numeric(plyr::mapvalues(GC_2019$obs, from=NewObsList[,1], to=NewObsList[,2]))


begin <- as.matrix(read.table("Data/begin.txt",header=T))

end <- as.matrix(read.table("Data/end.txt",header=T))


n <- as.matrix(read.table("Data/n.txt",header=T))
dim(n) <- c(179,2,6) #convert this back to a 3D array
#n <- abind(n,array(0,dim=c(2,2,6)), along=1) #add two trailing 0s to the end of the sightings array (this is for the day 1 and day 90 zero-whale anchor points)


n1 <- as.matrix(read.table("Data/n1.txt",header=T))      #These aren't needed
dim(n1) <- c(179,2,6) #convert this back to a 3D array


n2 <- as.matrix(read.table("Data/n2.txt",header=T))      #These aren't needed
dim(n2) <- c(179,2,6) #convert this back to a 3D array

u <- as.matrix(read.table("Data/u.txt",header=T))
dim(u) <- c(179,2,6) #convert this back to a 3D array
#u <- abind(u,array(0,dim=c(2,2,6)), along=1) #add two trailing 0s to the end of the effort on/off array
#for(i in 1:length(periods)){ #Place 1's for 'effort on' for the two periods following the end of true watches (this is for the day 1 and day 90 zero-whale anchor points)
#  u[(periods[i]+1):(periods[i]+2),,i] <- 1
#}



vs <- as.matrix(read.table("Data/vs.txt",header=T))
#vs <- rbind(vs,matrix(NA,nrow=2,ncol=6)) #Add two trailing NAs (this is for the day 1 and day 90 zero-whale anchor points)

bf <- as.matrix(read.table("Data/bf.txt",header=T))
#bf <- rbind(bf,matrix(NA,nrow=2,ncol=6)) #Add two trailing NAs (this is for the day 1 and day 90 zero-whale anchor points)

obs <- as.matrix(read.table("Data/obs.txt",header=T))
dim(obs) <- c(179,2,6) #convert this back to a 3D array
#obs <- abind(obs,array(0,dim=c(2,2,6)), along=1) #add two trailing 0s to the end of the effort on/off array
#for(i in 1:length(periods)){ #Place 36s for 'no observer' for the two periods following the end of true watches (this is for the day 1 and day 90 zero-whale anchor points)
#  obs[(periods[i]+1):(periods[i]+2),,i] <- 36 #this will force it to the mean observation probability with no observer effect
#}



# Add 2019 data onto previous data

# n
n2019 <- matrix(c(GC_2019$n,rep(0,length(GC_2019$n))), nrow=length(GC_2019$n), ncol=2)

n.all <- abind(abind(n,array(0,dim=c((dim(n2019)[1]-dim(n)[1]),2,6)),along=1),n2019,along=3)

# u
u2019 <- matrix(c(rep(1,length(GC_2019$n)),rep(0,length(GC_2019$n))), nrow=length(GC_2019$n), ncol=2)

u.all <- abind(abind(u,array(0,dim=c((dim(n2019)[1]-dim(n)[1]),2,6)),along=1),u2019,along=3)

# vs
vs.all <- cbind(rbind(vs,matrix(NA,nrow=dim(GC_2019)[1]-dim(vs)[1],ncol=6)),GC_2019$vs)

# bf
bf.all <- cbind(rbind(bf,matrix(NA,nrow=dim(GC_2019)[1]-dim(bf)[1],ncol=6)),GC_2019$bf)

# begin
begin.all <- cbind(rbind(begin,matrix(NA,nrow=dim(GC_2019)[1]-dim(begin)[1],ncol=6)),GC_2019$begin)

# end
end.all <- cbind(rbind(end,matrix(NA,nrow=dim(GC_2019)[1]-dim(end)[1],ncol=6)),GC_2019$end)

# obs
obs.2st.2019 <- matrix(c(obs2019,rep(36,length(obs2019))), nrow=length(obs2019), ncol=2)

obs.all <- abind(abind(obs,array(36,dim=c((dim(obs.2st.2019)[1]-dim(obs)[1]),2,6)),along=1),obs.2st.2019,along=3)


Watch.Length <- rbind(end.all,matrix(NA,nrow=2,ncol=7)) - rbind(begin.all,matrix(NA,nrow=2,ncol=7))
for(i in 1:length(periods)){ #Place 36s for 'no observer' for the two periods following the end of true watches (this is for the day 1 and day 90 zero-whale anchor points)
  Watch.Length[(periods[i]+1):(periods[i]+2),i] <- 1 #this will force it to the mean observation probability with no observer effect
}

# N inits
N_inits_john <- as.matrix(read.table("Data/N_inits.txt",header=T))
N_inits_temp <- rbind(N_inits_john, matrix(NA,nrow=17,ncol=6))
N_inits_2019 <- c(n.all[,1,7]*2,NA,NA)
N_inits <- cbind(N_inits_temp,N_inits_2019)  
  
#N_inits_temp <- n.all[,1,]*2

#N_inits <- rbind(N_inits_temp,matrix(NA,nrow=2,ncol=7))

for(i in 1:length(periods)){
  N_inits[(periods[i]+1):(periods[i]+2),i] <- NA #we're going to make N a partially observed data object with anchor points at day 1 and 90
#  N_inits[(periods[i]+3):dim(N_inits)[1],i] <- NA # this will give an error, subscript out of bounds on the final year, but it still did its job
}

N <- matrix(NA,nrow=dim(N_inits)[1],ncol=dim(N_inits)[2]) #The 'data' has to be the inverse of the inits, with NAs for all of the estimated Ns, and 0s for the days 1 and 90
for(i in 1:length(periods)){
  N[(periods[i]+1):(periods[i]+2),i] <- 0 #True number of whales passing fixed at 0 for day 1 and 90
}


day <- floor(begin.all)
t <- round((begin.all+end.all)/2)
#Add a couple of extra rows of NAs to the end of the day index reference to match up with the fixed 0s in N (above), assigning them to days 1 and 90
day <- rbind(as.matrix(day),matrix(NA,nrow=2,ncol=7))
for(i in 1:length(periods)){ #Set the anchor points: days 1 and 90
  day[(periods[i]+1):(periods[i]+2),i] <- c(1,90)
}

t <- rbind(as.matrix(t),matrix(NA,nrow=2,ncol=7))
for(i in 1:length(periods)){ #Set the anchor points: days 1 and 90
  t[(periods[i]+1):(periods[i]+2),i] <- c(1,90)
}


# Set up JAGS data objects and run the model:


jags.data <- list(n=n.all,
                  n.com=n.all,
                  n.sp=n.all,
                  n.station = dim(n.all)[2],
                  n.year = dim(n.all)[3],
                  n.obs = max(obs.all),
                  periods = periods,
                  obs=obs.all,
                  #Watch.Length = 0.0625,
                  u=u.all,
                  vs=vs.all,
                  bf=bf.all,
                  #day=day,
                  day=t,
                  N=N,
                  N.com=N,
                  N.sp=N,
                  knot=c(-1.46,-1.26,-1.02,-0.78,-0.58,-0.34,-0.10,0.10,0.34,0.57,0.78,1.02,1.26,1.46),
                  n.knots=14,
                  #begin=begin,
                  #end=end,
                  Watch.Length=Watch.Length)

jags.inits <- function() list(mean.prob = 0.5,
                              BF.Fixed = 0,
                              VS.Fixed = 0,
                              mean.beta = c(0,0,0), #mean.beta = c(5,0.14,-3.5),
                              beta.sigma = c(0.5,0.5,0.5),#beta.sigma = c(7,7,7),
                              BF.Switch = 0,
                              VS.Switch = 0,
                              OBS.Switch = 0,
                              sigma.Obs = 1,
                              N = N_inits,
                              N.com = N_inits,
                              N.sp = N_inits,
                              #z = matrix(1,nrow=90,ncol=6),
                              beta.sp = array(data=0,dim=c(2,7)),
                              sd.b.sp = c(1,1,1,1,1,1,1),
                              z=matrix(1,nrow=90,ncol=7))


library(R2jags)

parameters <- c("lambda","OBS.RF","OBS.Switch","BF.Switch","BF.Fixed","VS.Switch","VS.Fixed","mean.prob","mean.prob.com","mean.prob.sp","BF.Fixed.com","BF.Fixed.sp","VS.Fixed.com","VS.Fixed.sp",
                "Corrected.Est","Raw.Est","z","com","sp","Daily.Est","mean.beta","beta.sigma","beta","beta.sp","b.sp","sd.b.sp")

ni <- 20000
nt <- 20
nb <- 10000
nc <- 3


GW_Nmix_2019 <- jags(jags.data, inits=jags.inits, parameters, "GW_Nmix_Orig.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())


attach.jags(GW_Nmix_2019)

apply(Corrected.Est,2,median)


year <- 7

#John's plots:

UCIs <- apply(Daily.Est[,,year],2,quantile,0.975)
LCIs <- apply(Daily.Est[,,year],2,quantile,0.025)
plot(apply(exp(sp[,,year]),2,median),type='l',ylim=c(0,max(UCIs)+50))
lines(apply(exp(com[,,year]),2,median),type='l',lty=2)
segments(x0=1:90,y0=LCIs,y1=UCIs)

