# Compares Laake et al.'s analysis and proposed Richards' function analysis
# 
# 
# This R script runs Laake et al's model and computes abudnance estimates. Then,
# extracts gray whale sightings data from the ERAnalysis package and runs the 
# Richards' function JAGS models
# 
# Laake's model is found in the ERAnalysis package.
# 

library(ERAnalysis)  
library(tidyverse)
library(ggplot2)

data(PrimaryOff)
data(Primary)
data(ERSurveyData)
NorthYears=c(2000,2001)

# The data in PrimarySightings are all southbound sightings for all years in which visibility and beaufort
# are less than or equal to 4. Below the counts are shown for the 2 dataframes for
# recent surveys since 1987/88.
data(PrimarySightings)
data(PrimaryEffort)

# The following code produces the series of abundance estimates used in the paper and
# some comparative naive estimates.  


# The data from observer MAS in 1995 is excluded because this observer
# did not participate in the double count in that year. 
cat("\nMAS 1995/96 sightings = ",
    nrow(PrimarySightings[(PrimarySightings$Observer=="MAS" & PrimarySightings$Start.year==1995),]))
cat("\nMAS hours of effort in 1995/96= ",
    24*sum(PrimaryEffort$effort[(PrimaryEffort$Observer=="MAS" & PrimaryEffort$Start.year==1995)]))

PrimarySightings=PrimarySightings[!(PrimarySightings$Observer=="MAS" & 
                                      PrimarySightings$Start.year==1995),]
PrimaryEffort=PrimaryEffort[!(PrimaryEffort$Observer=="MAS" & 
                                PrimaryEffort$Start.year==1995),]

# Define arguments used in the analysis
all.years=unique(PrimaryEffort$Start.year)
recent.years=all.years[all.years>=1987]
early.years=all.years[all.years<1987]
final.time=sapply(tapply(floor(PrimaryEffort$time),
                         PrimaryEffort$Start.year,max),
                  function(x) ifelse(x>90,100,90))
lower.time=rep(0,length(final.time))
fn=1.0817
se.fn=0.0338

# Effort and sightings prior to 1987 were filtered for an entire watch if vis or beaufort 
# exceeded 4 at any time during the watch.  This is done for surveys starting in 1987 with the
# Use variable which is set to FALSE for all effort records in a watch if at any time the vis or
# beaufort exceeded 4 during the watch.
# 
# Here are the hours of effort that are excluded (FALSE) and included (TRUE) by each year
# Note that for most years <1987 there are no records with Use==FALSE because the filtered records
# were excluded at the time the dataframe was constructed. The only exception is for 1978 in which  
# one watch (5 hours) was missing a beaufort value so it was excluded.
tapply(PrimaryEffort$effort, 
       list(PrimaryEffort$Use, PrimaryEffort$Start.year),sum)*24

# These are the number of sightings that were included/excluded based on Use
Sightings=PrimarySightings
Sightings=merge(Sightings, subset(PrimaryEffort,select=c("key","Use")))
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
period.estimate=apply(whales.counted/sampled.fraction, 1, sum, na.rm=TRUE)

#plot(all.years,period.estimate,xlab="Survey year",ylab="Estimated gray whale population size")


# Compute naive estimates of abundance for the 23 surveys.  These use the uncorrected
# counts of whales from the primary observer during watches in which neither Beaufort nor
# vis exceeded 4.  For each year a gam with a smooth over time is fitted and this is
# used to predict total abundance throughout the migration from the counts of whales
# during the sampled periods.  There is no correction for missed pods or for measurement
# error in podsize. Each fitted migration gam is plotted with the observed values and
# saved in the file NaiveMigration.pdf.

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

Nhat.naive=sapply(naive.abundance.models,function(x) x$Total)
all.consec.years <- data.frame(Year = seq(from = min(all.years), to = max(all.years)))
Nhat.naive %>% 
  as.data.frame() %>%
  rownames_to_column(var = "Year_char") %>% #-> tmp
  rename(Nhat.naive = ".") %>% #-> tmp
  mutate(Year = as.numeric(Year_char)) %>% #-> tmp
  select(-Year_char) %>%
  right_join(all.consec.years, by = "Year") %>%
  arrange(Year) -> Nhat.df

p.Nhat.naive <- ggplot(Nhat.df) +
  geom_point(aes(x = Year, y = Nhat.naive)) +
  xlab("Survey year") +
  ylab("Estimated gray whale population size (naive)")

# plot(all.years,
#      Nhat.naive,
#      xlab="Survey year",
#      ylab="Estimated gray whale population size",pch=16)

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

if (!file.exists("RData/reilly_estimates.rds")){
  reilly.estimates=compute.series(models,
                                  naive.abundance.models,
                                  sightings=Sightings,
                                  effort=Effort,TruePS=FALSE)
  saveRDS(reilly.estimates,
          file = "RData/reilly_estimates.rds")  
} else {
  reilly.estimates <- readRDS("RData/reilly_estimates.rds")
}

Nhat.reilly=reilly.estimates$Nhat
Nhat.reilly[1:15]=Nhat.naive[1:15]*(Nhat.reilly[16]/Nhat.naive[16])
avg.Reilly.podsize=tapply(Sightings$corrected.podsize,Sightings$Start.year,mean)

Nhat.reilly %>%
  as.data.frame() %>%
  rownames_to_column(var = "Year_char") %>%
  rename(Nhat.reilly = ".") %>% #-> tmp
  mutate(Year = as.numeric(Year_char)) %>%
  select(-Year_char) %>%
  right_join(Nhat.df, by = "Year") %>%
  arrange(Year) -> Nhat.df

p.Nhat.reilly <- ggplot(Nhat.df) +
  geom_point(aes(x = Year, y = Nhat.reilly)) +
  xlab("Survey year") +
  ylab("Reilly-corrected estimated abundane of gray whales")

#plot(all.years,Nhat.reilly,xlab="Survey year",ylab="Estimated gray whale population size",pch=16)

#Next compute the series of abundance estimates for the most recent 8 years by
# fitting and selecting the best detection model but not applying the pod size correction.
# From those 8 estimates and the naive estimates, compute an average ratio and 
# apply it to generate the estimates for the first 15 surveys prior to 1987.
Sightings$corrected.podsize=Sightings$podsize

if (!file.exists("RData/Nhat_nops_correction.rds")){
  abundance.estimates.nops.correction=compute.series(models, 
                                                     naive.abundance.models,
                                                     sightings=Sightings,
                                                     effort=Effort,
                                                     TruePS=FALSE)
  saveRDS(abundance.estimates.nops.correction,
          "RData/Nhat_nops_correction.rds")  
} else {
  abundance.estimates.nops.correction <- readRDS("RData/Nhat_nops_correction.rds")
}

Sightings$corrected.podsize=NULL
# cat("\nCorrection factor for missed pods without pod size correction\n")
# print(cbind(Year=all.years, ratio=abundance.estimates.nops.correction$Nhat/Nhat.naive/fn))

# Next compute the series of abundance estimates for the most recent 8 years by
# fitting and selecting the best detection model.  From those 8 estimates and the
# naive estimates, compute an average ratio and apply it to generate the estimates
# for the first 15 surveys prior to 1987. Note with hessian=TRUE, the analysis can
# take about 30-60 minutes to complete. (TE: Takes about 6 minutes now. But added
# the if-else. 2023-08-31)
if (!file.exists("RData/Laake_abundance_estimates.rds")){
  abundance.estimates=compute.series(models,
                                     naive.abundance.models,
                                     sightings=Sightings,
                                     effort=Effort,
                                     hessian=TRUE)

  saveRDS(abundance.estimates,
          file = "RData/Laake_abundance_estimates.rds")  
} else {
  abundance.estimates <- readRDS("RData/Laake_abundance_estimates.rds")
}

abundance.estimates$Nhat %>%
  as.data.frame() %>%
  rownames_to_column(var = "Year_char") %>%
  rename(Nhat.corrected = ".") %>% #-> tmp
  mutate(Year = as.numeric(Year_char)) %>%
  select(-Year_char) %>%
  right_join(Nhat.df, by = "Year") %>%
  arrange(Year) -> Nhat.df

p.Nhat.corrected <- ggplot(Nhat.df) +
  geom_point(aes(x = Year, y = Nhat.corrected)) +
  xlab("Survey year") +
  ylab("Corrected estimated abundane of gray whales")

# plot(all.years,
#      abundance.estimates$Nhat,xlab="Survey year",
#      ylab="Estimated gray whale population size")

cat("\nCorrection factor for pod size error\n")
print(cbind(Year=all.years, 
            ratio=abundance.estimates$Nhat/abundance.estimates.nops.correction$Nhat))

# Print out estimates for Table 7 in Laake et al. -- pod size and detection  
for(i in 1:8)
  print(cbind(paste(signif(abundance.estimates$det[[i]][[1]]$par,digits=3)," (",
                    signif(sqrt(diag(solve(abundance.estimates$det[[i]][[1]]$model$hessian))),
                           digits=3),")",sep="")))

# Next compute E(S) and it's std error for Table 7
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

# a function from podsize.calibration.matrix.r
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

# A function from psfit.gamma.R
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

data(gsS)

for(i in 1:8)
  print( ES.var(abundance.estimates$detection.models[[i]][[1]]$model$par,
                abundance.estimates$detection.models[[i]][[1]]$model$hessian,gsS))

# Plot for each year average observed pod size and previously used corrected pod size 
avg.podsize=tapply(Sightings$podsize, Sightings$Start.year,mean)
avg.podsize %>%
  as.data.frame() %>%
  rownames_to_column(var = "Year_char") %>%
  rename(avg.podsize = ".") %>% #-> tmp
  mutate(Year = as.numeric(Year_char)) %>%
  select(-Year_char) %>%
  right_join(Nhat.df, by = "Year") %>%
  arrange(Year) -> Nhat.df

avg.Reilly.podsize %>%
  as.data.frame() %>%
  rownames_to_column(var = "Year_char") %>%
  rename(avg.Reilly.podsize = ".") %>% #-> tmp
  mutate(Year = as.numeric(Year_char)) %>%
  select(-Year_char) %>%
  right_join(Nhat.df, by = "Year") %>%
  arrange(Year) -> Nhat.df

p.avg.podsize <- ggplot(Nhat.df) +
  geom_point(aes(x = Year, y = avg.podsize),
             color = "red") +
  geom_point(aes(x = Year, y = avg.Reilly.podsize),
             color = "blue") +
  xlab("Survey year") +
  ylab("Average gray whale pod size (red = average, blue = Reilly correction")

#legend(1970,1.25,pch=c(1,2),legend=c("Observed average","Reilly correction"))

# Next plot the averages for the last 8 years after linking pods, Reilly corrected value and then the expected value
# from the fitted gamma distribution

pod.size.df <- data.frame(Year = rep(recent.years, 3),
                          Pod.size = c(abundance.estimates$summary.df$MeanPS,
                                       sapply(abundance.estimates$detection.models,
                                              function(x) expected.podsize(x[[1]]$par[1:2],gsS,20)$ES),
                                       avg.Reilly.podsize[16:23]),
                          Method = c(rep("Observed average", length(recent.years)),
                                     rep("Reilly correction", length(recent.years)),
                                     rep("Expected pod size", length(recent.years))))

p.pod.size <- ggplot(pod.size.df) +
  geom_point(aes(x = Year, y = Pod.size,
                 color = Method)) +
  xlab("Survey year") + 
  ylab("Average gray whale pod size")
  
# plot(recent.years, 
#      abundance.estimates$summary.df$MeanPS, 
#      ylim=c(1,3), 
#      xlab="Survey Year (1 Dec)",
#      ylab="Average gray whale pod size",pch=1)

# expected.ps=sapply(abundance.estimates$detection.models,
#                    function(x) expected.podsize(x[[1]]$par[1:2],gsS,20)$ES)
# 
# points(recent.years,expected.ps, pch=16)
# points(recent.years,avg.Reilly.podsize[16:23], pch=2)
# legend(1990,1.35,pch=c(1,2,16),legend=c("Observed average", 
# "Reilly correction","Expected pod size"))

# Ran the same models on newer data since 2006.
# I do not have raw count data from 2007 - 2014. Surveys were conducted in
# 2007/2008, 2009/2010, and 2010/2011. For those years, only daily sum
# are avaialable. Raw counts are available from 2014/2015, 2015/2016, 2019/2020,
# 2021/2022, and 2022/2023 

# this file contains all necessary inputs:
data.0 <- readRDS("RData/2006-2019_GC_Formatted_Data.RDS")

# output from Ver2.0 extraction
years <- c("2015", "2016", "2020", "2022", "2023")


