# Runs Laake's example code
# 
# # ERAnalysis library had to be rebuilt for the current version of R. To do so,
# (Thanks to Jeff Laake!)
# 1. Install devtools package
# 2. Copy ERAnalysis folder into R working directory (the inner ERAnalysis folder 
# in his GitHub repository: https://github.com/jlaake/ERAnalysis/tree/master/
# 3. Create a new project - existing folder - select ERAnalysis
# 4. Click on the "Build" tab
# 5. Click "Install"
# 
# It may complain about not having certain libraries. Install them and re-install.
# 
# 
rm(list=ls())
library(ERAnalysis)  

#Run example code from the help file
#
# The recent survey data 1987 and after are stored in ERSurveyData and those data
# are processed by the ERAbund program to produce files of sightings and effort.
# The sightings files are split into Primary observer and Secondary observer sightings.
# Primary observer sightings are whales that are not travelling North and are defined by
# those when EXPERIMENT==1 (single observer) or a designated LOCATION when EXPERIMENT==2.
#  For surveys 2000/2001 and 2001/2002, the primary observer was at LOCATION=="N"
# and for all other years, LOCATION=="S".
#
# Based on the projected timing of the passage of the whale (t241) perpendicular to the 
# watch station, the sighting was either contained in the watch (on effort) or not (off effort).
# The dataframe Primary contains all of the on effort sightings and PrimaryOff contains all
# of the off-effort sightings.  The code below shows that the counts of primary
# sighting records matches the total counts of the on and off effort sightings split into
# the 2 dataframes.
#
data(PrimaryOff)
data(Primary)
data(ERSurveyData)
NorthYears=c(2000,2001)
rowSums(with(ERSurveyData[ERSurveyData$EFLAG==3 & ERSurveyData$TRAVELDIR!="N"&
                            (ERSurveyData$EXPERIMENT==1 | 
                               (ERSurveyData$EXPERIMENT==2& 
                                  ((ERSurveyData$Start.year%in%NorthYears & 
                                      ERSurveyData$LOCATION=="N") | 
                                     (!ERSurveyData$Start.year%in%NorthYears & 
                                        ERSurveyData$LOCATION=="S")))),],
             table(Start.year,LOCATION)))

table(Primary$Start.year) + table(PrimaryOff$Start.year)


# Likewise, the secondary sightings are those with EXPERIMENT==2 but the LOCATION that
# is not designated as primary.  The following show that the counts match for ERSurveyData
# and the dataframe SecondarySightings that was created by ERAbund
data(SecondarySightings)
{rowSums(with(ERSurveyData[ERSurveyData$EFLAG==3&ERSurveyData$TRAVELDIR!="N"&
                             (ERSurveyData$EXPERIMENT==2& ( 
                               (ERSurveyData$Start.year%in%NorthYears & ERSurveyData$LOCATION=="S") | 
                                 (!ERSurveyData$Start.year%in%NorthYears & (ERSurveyData$LOCATION=="N" | 
                                                                              ERSurveyData$LOCATION=="NP")))),], 
              table(Start.year,LOCATION)))}

table(SecondarySightings$Start.year)

{rowSums(with(ERSurveyData[ERSurveyData$EFLAG==3&ERSurveyData$TRAVELDIR!="N"&
                             ERSurveyData$VISCODE<=4&ERSurveyData$WINDFORCE<=4&
                             (ERSurveyData$EXPERIMENT==1 | 
                                (ERSurveyData$EXPERIMENT==2& ( 
                                  (ERSurveyData$Start.year%in%NorthYears & ERSurveyData$LOCATION=="N") | 
                                    (!ERSurveyData$Start.year%in%NorthYears & ERSurveyData$LOCATION=="S")))),],
              table(Start.year,LOCATION)))}



# The data in PrimarySightings are all southbound sightings for all years in which visibility and beaufort
# are less than or equal to 4. Below the counts are shown for the 2 dataframes for
# recent surveys since 1987/88.
table(Primary$Start.year[Primary$vis<=4 & Primary$beaufort<=4])
data(PrimarySightings)
table(PrimarySightings$Start.year[PrimarySightings$Start.year>=1987])

# The following code produces the series of abundance estimates used in the paper and
# some comparative naive estimates.  It will also optionally (commented out) compute
# the variance-covariance matrix for the estimates.
#
# Get sightings and effort data
data(PrimarySightings)
data(PrimaryEffort)

# The data from observer MAS in 1995 is excluded because this observer
# did not participate in the double count in that year. 
cat("\nMAS 1995/96 sightings = ",nrow(PrimarySightings[(PrimarySightings$Observer=="MAS"&PrimarySightings$Start.year==1995),]))
cat("\nMAS hours of effort in 1995/96= ",24*sum(PrimaryEffort$effort[(PrimaryEffort$Observer=="MAS"&PrimaryEffort$Start.year==1995)]))

PrimarySightings=PrimarySightings[!(PrimarySightings$Observer=="MAS"&PrimarySightings$Start.year==1995),]
PrimaryEffort=PrimaryEffort[!(PrimaryEffort$Observer=="MAS"&PrimaryEffort$Start.year==1995),]
table(PrimarySightings$Start.year)
# Define arguments used in the analysis
all.years=unique(PrimaryEffort$Start.year)
recent.years=all.years[all.years>=1987]
early.years=all.years[all.years<1987]
final.time=sapply(tapply(floor(PrimaryEffort$time),PrimaryEffort$Start.year,max),function(x) ifelse(x>90,100,90))
lower.time=rep(0,length(final.time))
fn=1.0817
se.fn=0.0338

# Effort and sightings prior to 1987 were filtered for an entire watch if vis or beaufort 
# exceeded 4 at any time during the watch.  This is done for surveys starting in 1987 with the
# Use variable which is set to FALSE for all effort records in a watch if at any time the vis or
# beaufort exceeded 4 during the watch.
# Here are the hours of effort that are excluded (FALSE) and included (TRUE) by each year
# Note that for most years <1987 there are no records with Use==FALSE because the filtered records
# were excluded at the time the dataframe was constructed. The only exception is for 1978 in which  
# one watch (5 hours) was missing a beaufort value so it was excluded.
tapply(PrimaryEffort$effort,list(PrimaryEffort$Use,PrimaryEffort$Start.year),sum)*24

# These are the number of sightings that were included/excluded based on Use
Sightings=PrimarySightings
Sightings=merge(Sightings,subset(PrimaryEffort,select=c("key","Use")))
tapply(Sightings$Start.year,list(Sightings$Use,Sightings$Start.year),length)

# Filter effort and sightings and store in dataframes Effort and Sightings
Effort=PrimaryEffort[PrimaryEffort$Use,]  
Sightings=PrimarySightings
Sightings$seq=1:nrow(Sightings)
Sightings=merge(Sightings,subset(Effort,select=c("key")))
Sightings=Sightings[order(Sightings$seq),]

# Using the filtered data, compute simple minded abundance estimates by treating the 
# sampled periods (watches) as a random sample of a period of a specified number
# of days (e.g., 4 days).  It uses raw counts with no correction for pod size or missed pods.
period=4
# compute the fraction of each period that was sampled (eg. 12 hours / (4 days *24 hrs/day))
sampled.fraction=with(Effort,
                      {
                        Day=as.numeric(as.Date(Date)-as.Date(paste(Start.year,"-12-01",sep="")))
                        tapply(effort,list(Start.year,cut(Day,seq(0,100,period))),sum)/period
                      })
# compute the number of whales counted in each period
whales.counted=with(Sightings,
                    {
                      Day=as.numeric(as.Date(Date)-as.Date(paste(Start.year,"-12-01",sep="")))
                      tapply(podsize,list(Start.year,cut(Day,seq(0,100,period))),sum)
                    })

# Compute simple minded population estimate and plot it
period.estimate=apply(whales.counted/sampled.fraction,1,sum,na.rm=TRUE)
#pdf("SimpleMindedEstimates.pdf")
plot(all.years,period.estimate,xlab="Survey year",ylab="Estimated gray whale population size")
#dev.off()

# Compute naive estimates of abundance for the 23 surveys.  These use the uncorrected
# counts of whales from the primary observer during watches in which neither Beaufort nor
# vis exceeded 4.  For each year a gam with a smooth over time is fitted and this is
# used to predict total abundance throughout the migration from the counts of whales
# during the sampled periods.  There is no correction for missed pods or for measurement
# error in podsize. Each fitted migration gam is plotted with the observed values and
# saved in the file NaiveMigration.pdf.
#pdf("NaiveMigration.pdf")
naive.abundance.models=vector("list",23)
i=0
for (year in all.years)
{
  i=i+1
  primary=Sightings[Sightings$Start.year==year,]
  primary$Start.year=factor(primary$Start.year)
  ern=subset(Effort,
             subset=as.character(Start.year)==year,
             select=c("Start.year","key","begin","end","effort","time","vis","beaufort"))
  
  ern$Start.year=factor(ern$Start.year)
  naive.abundance.models[[i]]=estimate.abundance(spar=NULL,
                                                 dpar=NULL,
                                                 gsS=gsS,
                                                 effort=ern, 
                                                 sightings=primary, 
                                                 final.time=final.time[i],
                                                 lower.time=lower.time[i],
                                                 gformula=~s(time),
                                                 dformula=NULL)
}
#dev.off()

Nhat.naive=sapply(naive.abundance.models,function(x) x$Total)
#pdf("NaiveMigrationEstimates.pdf")
plot(all.years,
     Nhat.naive,
     xlab="Survey year",
     ylab="Estimated gray whale population size",pch=16)
#dev.off()

# Define set of models to be evaluated for detection
models=c("podsize+Dist+Observer",
         "podsize+Dist+Observer+beaufort",
         "podsize+Dist+Observer+vis",
         "podsize+Dist+Observer+Vis")

# Create time series of estimates based on previous approach using Reilly pod size
# correction method but using 1978 data for surveys <=1987 and 1992-1994 aerial data
# for surveys >=1992
data(add.cf.reilly)
data(add.cf.laake)
Sightings$corrected.podsize[Sightings$Start.year<=1987]=reilly.cf(Sightings$podsize[Sightings$Start.year<=1987],
                                                                  add.cf.reilly)
Sightings$corrected.podsize[Sightings$Start.year>1987]=reilly.cf(Sightings$podsize[Sightings$Start.year>1987],
                                                                 add.cf.laake)
#pdf("ReillyApproach.pdf")
reilly.estimates=compute.series(models,
                                naive.abundance.models,
                                sightings=Sightings,
                                effort=Effort,TruePS=FALSE)

Nhat.reilly=reilly.estimates$Nhat
Nhat.reilly[1:15]=Nhat.naive[1:15]*(Nhat.reilly[16]/Nhat.naive[16])
avg.Reilly.podsize=tapply(Sightings$corrected.podsize,Sightings$Start.year,mean)
#dev.off()

#pdf("ReillyEstimates.pdf")
plot(all.years,Nhat.reilly,
     xlab="Survey year",ylab="Estimated gray whale population size",pch=16)
#dev.off()

#Next compute the series of abundance estimates for the most recent 8 years by
# fitting and selecting the best detection model but not applying the pod size correction.
# From those 8 estimates and the naive estimates, compute an average ratio and 
# apply it to generate the estimates for the first 15 surveys prior to 1987.
Sightings$corrected.podsize=Sightings$podsize
abundance.estimates.nops.correction=compute.series(models, 
                                                   naive.abundance.models,
                                                   sightings=Sightings,
                                                   effort=Effort,
                                                   TruePS=FALSE)
Sightings$corrected.podsize=NULL
cat("\nCorrection factor for missed pods without pod size correction\n")
print(cbind(Year=all.years,ratio=abundance.estimates.nops.correction$Nhat/Nhat.naive/fn))

# Next compute the series of abundance estimates for the most recent 8 years by
# fitting and selecting the best detection model.  From those 8 estimates and the
# naive estimates, compute an average ratio and apply it to generate the estimates
# for the first 15 surveys prior to 1987. Note with hessian=TRUE, the analysis can
# take about 30-60 minutes to complete. (TE: Takes about 6 minutes now. But added
# the if-else. 2023-08-31)
if (!file.exists("RData/Laake_abundance_estimates.rds")){
#  pdf("Migration.pdf")
  abundance.estimates=compute.series(models,naive.abundance.models,
                                     sightings=Sightings,
                                     effort=Effort,
                                     hessian=TRUE)
#  dev.off()
 saveRDS(abundance.estimates,
          file = "RData/Laake_abundance_estimates.rds")  
} else {
  abundance.estimates <- readRDS("RData/Laake_abundance_estimates.rds")
}

#pdf("CurrentEstimates.pdf")
plot(all.years,
     abundance.estimates$Nhat,xlab="Survey year", ylab="Estimated gray whale population size")
#dev.off()
cat("\nCorrection factor for pod size error\n")
print(cbind(Year=all.years, ratio=abundance.estimates$Nhat/abundance.estimates.nops.correction$Nhat))

# Print out estimates for Table 7 in Laake et al. -- pod size and detection  
for(i in 1:8)
  print(cbind(paste(signif(abundance.estimates$det[[i]][[1]]$par,digits=3)," (",
                    signif(sqrt(diag(solve(abundance.estimates$det[[i]][[1]]$model$hessian))),digits=3),")",sep="")))

# Next compute E(S) and it's std error for Table 7
#
ES.var=function(par,hessian,gsS,nmax=20)
{
  spar=par[1:2]
  vc=solve(hessian)[1:2,1:2]
  sparest=spar
  sparest[1]=spar[1]*0.999
  s1.low=expected.podsize(sparest,gsS=gsS,nmax=nmax)$ES
  sparest[1]=spar[1]*1.001
  s1.hi=expected.podsize(sparest,gsS=gsS,nmax=nmax)$ES
  sparest=spar
  sparest[2]=spar[2]*0.999
  s2.low=expected.podsize(sparest,gsS=gsS,nmax=nmax)$ES
  sparest[2]=spar[2]*1.001
  s2.hi=expected.podsize(sparest,gsS=gsS,nmax=nmax)$ES
  deriv.s=c((s1.hi-s1.low)/(.002*spar[1]),(s2.hi-s2.low)/(.002*spar[2]))
  return(list(ES=expected.podsize(spar,gsS=gsS,nmax=nmax)$ES,se=sqrt(t(deriv.s)%*%vc%*%deriv.s)))
}

data(gsS)

# in podsize.calibration.matrix.r
expected.podsize=function(spar,gsS,nmax=20)
{
  # Compute true pod size distribution and ES
  fS=gammad(spar,nmax)
  ES=sum(fS*(1:nmax))
  # Compute conditional distribution and conditional expectation for
  # true pod size given an observed pod size
  fSs=fS*gsS/matrix(colSums(fS*gsS),nrow=nmax,ncol=nmax,byrow=TRUE)
  ESs=apply(fSs,2,function(x) sum((1:nmax)*x))
  return(list(fS=fS,ES=ES,fSs=fSs,ESs=ESs))
}

# in psfit.gamma.R
gammad <- function(par,nmax,True=NULL,shape=TRUE)
{
  # This function creates a discretized gamma function from 1 to nmax
  #
  # Arguments:
  #
  #  par  - parameters for gamma
  #  nmax - maximum pod size
  #  True - values of true pod size for a range of pod sizes (eg 4+) in which
  #          a common distribution is being fitted with a relationship on
  #          either the shape (if shape=TRUE), or the scale=1/rate
  #  shape - TRUE/FALSE; see above
  #
  # The parameters:
  #  True = NULL then log(shape)=par[1], log(rate)=par[2] and scale=1/rate
  #  True is not null and !shape then log(shape)=par[1], log(rate)=-par[2]-par[3]*True
  #  True is not null and shape=TRUE then log(shape)=par[2]+par[3]*True, log(rate)=par[1]
  if(is.null(True))
  {
    pp=pgamma(1:nmax,shape=exp(par[1]),rate=exp(par[2]))
    ps=diff(c(0,pp))/pp[nmax]
    ps[ps==0]=1e-12
    ps=ps/sum(ps)
  }
  else
  {
    if(!shape)
    {
      pp= do.call("rbind",lapply(True, function(x) pgamma(1:nmax,shape=exp(par[1]),rate=exp(-par[2]-par[3]*x))))
      ps=t(apply(cbind(rep(0,nrow(pp)),pp),1,diff))/pp[,nmax]
    }
    else
    {
      pp= do.call("rbind",lapply(True, function(x) pgamma(1:nmax,shape=exp(par[2]+par[3]*x),rate=exp(par[1]))))
      ps=t(apply(cbind(rep(0,nrow(pp)),pp),1,diff))/pp[,nmax]
    }
    ps[ps==0]=1e-12
    ps=ps/rowSums(ps)
  }
  return(ps)
}


for(i in 1:8)
  print( ES.var(abundance.estimates$detection.models[[i]][[1]]$model$par,
                abundance.estimates$detection.models[[i]][[1]]$model$hessian,gsS))

# Error in expected.podsize(sparest, gsS = gsS, nmax = nmax) :
#   could not find function "expected.podsize"
#   
#   The function is in podsize.calibration.matrix.r. The package directs the working
#   directory to Laake's desktop... that may be a problem. 

# NOT RUN: Compute the var-cov matrix of the estimates.  These are not run by default
# because it takes a very long time to run. To run copy code from documentation into R.
## Not run: 
#
# Load the fitted mixed-effects gamma model with a random pod effect -- this was the most parsimonious model
data(gamma.pod)
# #
# ps.results=gamma.pod
# ps.results$vc=solve(ps.results$hessian)
# abundance.vc=compute.series.var(abundance.estimates,naive.abundance.models,
#                                 ps.results=ps.results,delta=c(0.001,0.01))
# write.table(log(abundance.estimates$Nhat),"nhat.txt")
# 
# lnvc=log(1+abundance.vc$vc/outer(abundance.estimates$Nhat,
#                                  abundance.estimates$Nhat,"*"))
# write.table(lnvc,"vc.txt")
# # The following produces Fig 4 in Laake et al.
# pdf("EstimatesWithCL.pdf")
# gw.abund <-
#   structure(list(year = c(1967, 1968, 1969, 1970, 1971, 1972, 1973,
#                           1974, 1975, 1976, 1977, 1978, 1979, 1984, 1985, 1987, 1992, 1993,
#                           1995, 1997, 2000, 2001, 2006), 
#                  N.old = c(13776, 12869, 13431, 11416,
#                            10406, 16098, 15960, 13812, 15481, 16317, 17996, 13971, 17447,
#                            22862, 21444, 22250, 18844, 24638, 24065, 29758, 19448, 18178,
#                            20110), 
#                  se.N.old = c(1082, 708, 758, 590, 614, 834, 872, 781,
#                               930, 818, 1249, 753, 984, 1379, 1120, 1115, 1190, 1475, 1393,
#                               3122, 1882, 1780, 1766)), 
#             .Names = c("year","N.old", "se.N.old"), row.names = c(NA, -23L), class = "data.frame")
#             
#                                                                                                                              
# conf.int=function(abundance, CV, alpha=0.05, digits=2, prt=FALSE)
# {
#   # Computes confidence intervals based on lognormal distr.
#   # JMB / NMML / 11 Sep 2008
#   
#   if (alpha <0 || alpha > .999) stop("alpha must be in (0,1)")
#   z = round(abs(qnorm(alpha/2)),2)
#   if (prt) cat("N:",abundance,"  cv:",CV,"  alpha:",alpha,"  z:",z,"\n")
#   C <- exp(z * sqrt(log(1 + CV^2)))
#   SL <- round(abundance/C,digits)
#   SU <- round(abundance * C,digits)
#   data.frame(SL,SU)
# }
# 
# gwplot=function (N,se,year,bar.wid=.5,p=NULL)
# {
#   #   Function written by JMB with mods by JLL
#   require("ggplot2")
#   # Plot gray whale abund. data using ggplot2 package
#   if (length(N) != length(se) | length(se) != length(year)) stop("error in vector lengths")
#   cv <- se/N 
#   yy <- conf.int(N,cv)
#   y.max <- yy$SU ; y.min <- yy$SL
#   yy <- cbind(y.min,y.max)
#   gwd <- data.frame(year,N) 
#   gwd <- cbind(gwd,yy)
#   names(gwd) <- c("Year","Abundance","y.min","y.max")
#   if(is.null(p))
#   {
#     limits <- aes(ymax = y.max, ymin = y.min)   
#     p <- ggplot(gwd,aes(x=Year,y=Abundance))
#     p <- p + geom_point(colour="darkblue",data=gwd) +
#       geom_errorbar(limits,colour="blue",width=bar.wid,data=gwd)
#   }
#   else
#   {
#     limits <- aes(ymax = y.max, ymin = y.min)   
#     p <- p + geom_point(colour="black",data=gwd) +
#       geom_errorbar(limits,colour="black",linetype=2,width=bar.wid,data=gwd)
#   }
#   p
# }
# 
# xx=gwplot(abundance.estimates$Nhat,
#           abundance.vc$se,
#           gw.abund$year)
# gwplot(gw.abund$N.old,
#        gw.abund$se.N.old,
#        gw.abund$year-.5,
#        p=xx)
# 
# dev.off()
# 
# # Alternate runs using different matching weights
# abundance.estimates.lo=compute.series(models,naive.abundance.models,
#                                       sightings=Sightings,
#                                       effort=Effort,hessian=TRUE,twt=0.11,dwt=3.02)
# 
# write.table(log(abundance.estimates.lo$Nhat),"nhat_lo.txt")
# 
# abundance.estimates.hi=compute.series(models,
#                                       naive.abundance.models,
#                                       sightings=Sightings,
#                                       effort=Effort,
#                                       hessian=TRUE,twt=0.27,dwt=5.06)
# 
# write.table(log(abundance.estimates.hi$Nhat),"nhat_hi.txt")
# 
# abundance.vc.lo=compute.series.var(abundance.estimates.lo,
#                                    naive.abundance.models,
#                                    ps.results=ps.results,
#                                    twt=0.11,dwt=3.02,delta=c(0.001,0.01))
# 
# lnvc=log(1+abundance.vc.lo$vc/outer(abundance.estimates.lo$Nhat,
#                                     abundance.estimates.lo$Nhat,"*"))
# 
# write.table(lnvc,"vc_lo.txt")
# 
# abundance.vc.hi=compute.series.var(abundance.estimates.hi,
#                                    naive.abundance.models,
#                                    ps.results=ps.results,
#                                    twt=0.27,dwt=5.06,delta=c(0.001,0.02))
# 
# lnvc=log(1+abundance.vc.hi$vc/outer(abundance.estimates.hi$Nhat,
#                                     abundance.estimates.hi$Nhat,"*"))
# write.table(lnvc,"vc_hi.txt")

## End(Not run)
## 

# Plot for each year average observed pod size and previously used corrected pod size 
data(gsS)
#pdf("PodsizeComparisons.pdf")
avg.podsize=tapply(Sightings$podsize, Sightings$Start.year,mean)
plot(all.years, 
     avg.podsize, 
     ylim=c(1,3), 
     xlab="Survey Year (1 Dec)",
     ylab="Average gray whale pod size",pch=1)

points(all.years, avg.Reilly.podsize, pch=2)

legend(1970,1.25,pch=c(1,2),legend=c("Observed average","Reilly correction"))

# Next plot the averages for the last 8 years after linking pods, Reilly corrected value and then the expected value
# from the fitted gamma distribution
plot(recent.years, 
     abundance.estimates$summary.df$MeanPS, 
     ylim=c(1,3), 
     xlab="Survey Year (1 Dec)",
     ylab="Average gray whale pod size",pch=1)

# Error in expected.podsize(x[[1]]$par[1:2], gsS, 20) : 
# could not find function "expected.podsize"
# Error in gammad(spar, nmax) : could not find function "gammad"
# These functions are in podsize.calibration.matrix.r and psfit.gamma.R, respectively.
# So, I don't know why these errors return... I loaded the functions manually in
# the workspace and it ran fine. 
expected.ps=sapply(abundance.estimates$detection.models,
                   function(x) expected.podsize(x[[1]]$par[1:2],gsS,20)$ES)


points(recent.years,expected.ps, pch=16)
points(recent.years,avg.Reilly.podsize[16:23], pch=2)
legend(1990,1.35,pch=c(1,2,16),legend=c("Observed average","Reilly correction","Expected pod size"))
#dev.off()


# Plot distributions of pod size in calibration data and for true pod size from that year
# This is Figure 5 in Laake et al.
data(PodsizeCalibrationTable)
obs.1992=with(PodsizeCalibrationTable[substr(PodsizeCalibrationTable$key,1,11)=="Aerial_1992",],
              table(cut(True,c(1,2,3,4,20),right=FALSE)))
obs.1993=with(PodsizeCalibrationTable[substr(PodsizeCalibrationTable$key,1,11)=="Aerial_1993",],
              table(cut(True,c(1,2,3,4,20),right=FALSE)))
obs.1997=with(PodsizeCalibrationTable[substr(PodsizeCalibrationTable$key,1,8)=="Tracking",],
              table(cut(True,c(1,2,3,4,20),right=FALSE)))
obs.1992=obs.1992/sum(obs.1992)
obs.1993=obs.1993/sum(obs.1993)
obs.1997=obs.1997/sum(obs.1997)
ps.1992=gammad(c(-.073,-.3474),20)
ps.1993=gammad(c(-.07,-.474),20)
ps.1997=gammad(c(-.598,-.674),20)
ps.1992=c(ps.1992[1:3],sum(ps.1992[4:20]))
ps.1993=c(ps.1993[1:3],sum(ps.1993[4:20]))
ps.1997=c(ps.1997[1:3],sum(ps.1997[4:20]))
par(mfrow=c(3,1))
barplot(rbind(obs.1992,ps.1992),beside=TRUE,space=c(0,.1),
        main="1992-93",ylab="Proportion")
barplot(rbind(obs.1993,ps.1993),beside=TRUE,space=c(0,.1),
        main="1993-94",ylab="Proportion")
barplot(rbind(obs.1997,ps.1997),beside=TRUE,space=c(0,.1),
        main="1997-98",ylab="Proportion",
        xlab="True pod size",
        legend.text=c("Calibration","Estimated"),args.legend=list(x=8,y=.8))
# Code to construct Figure 2
par(mfrow=c(2,2))
data(PodsizeCalibrationTable)
psdf=PodsizeCalibrationTable
ps1=as.matrix(psdf[psdf$True==1,3:22])
ps2=as.matrix(psdf[psdf$True==2,3:22])
ps3=as.matrix(psdf[psdf$True==3,3:22])
psplus=as.matrix(psdf[psdf$True>=4,3:22])
gpar=gamma.pod$par
ps.results=gamma.pod
expect.ps=create.gsS(ps.results,True=psdf$True[psdf$True>=4][1:21])
expect.ps=colMeans(expect.ps)
barplot(rbind(gsS[1,],colMeans(ps1/rowSums(ps1))),beside=TRUE,main="True size = 1",legend.text=c("exp","obs"))
barplot(rbind(gsS[2,],colMeans(ps2/rowSums(ps2))),beside=TRUE,main="True size = 2",legend.text=c("exp","obs"))
barplot(rbind(gsS[3,],colMeans(ps3/rowSums(ps3))),beside=TRUE,main="True size = 3",legend.text=c("exp","obs"))
barplot(rbind(expect.ps,colMeans(psplus/rowSums(psplus))),
        beside=TRUE,
        main=paste("True size = 4+",sep=""),
        legend.text=c("exp","obs"))

