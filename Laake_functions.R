compute.sampling.multiplier=function(model,data,lower=0,upper=90)
{
   TotalArea=integrate(integrate.gam,lower=lower,upper=upper,model=model,vis=1,
            beauf=0)$val
   int.areas=apply(data[,c("begin","end","vis","beaufort")],1,function(x)
            integrate(integrate.gam,lower=x[1],upper=x[2],model=model,vis=as.vector(x[3]),
            beauf=as.vector(x[4]))$val)
   return(TotalArea/sum(int.areas))
}
integrate.gam=function(x,model,vis,beauf)
{
# Function for computing area under portions of the gam migration model
#
# Arguments:
#    x     -  vector of times to be evaluated within range of migration curve
#  model   -  gam model fitted to migration curve data
#  vis     - vis values to use for prediction if needed
#  beauf   - beaufort values to use for prediction if needed
#
 newdata=data.frame(time=x,vis=vis,beauf=beauf)
 return(predict(model,newdata=newdata,type="response"))
}

compute.series=function(models, naive.abundance,sightings=NULL,effort=NULL,
                        gsS=NULL,best=TRUE,Match=NULL,cutoff=4,
                        final.time=c(rep(90,20),100,100,90),
                        lower.time=rep(0,23), lcut=0.2, mcut=1.0, 
                        twt=0.18, dwt=3.95, pwt=0.05,
                        crittype="ellipse",
                        DistBreaks=c(0,1,2,3,4,20),Use=TRUE,
                        hessian=FALSE,debug=FALSE,TruePS=TRUE,
                        recent.years=c(1987,1992,1993,1995,1997,2000,2001,2006),fn=1.0817)
{
# Define function to select best detection model and return list of models
# within a cutoff delta AICc value.
select.detection.models=function(x,models,cutoff=4)
{
   mod=vector("list",length(models))      
   for(i in 1:length(models))
      mod[[i]]=io.glm(x,as.formula(paste("seen~",models[i],sep="")))
   AICValues=sapply(mod,function(x) AIC(x))
   DeltaAIC=AICValues-min(AICValues)
   mod.numbers=(1:length(models))[order(DeltaAIC)]
   return(mod[mod.numbers[sort(DeltaAIC)<cutoff]])
}
# Get sightings and effort
# Assign day value and match up to Use field if Use=TRUE and select only those
# sightings with Use=TRUE.
PrimarySightings=NULL
PrimaryEffort=NULL
if(is.null(sightings))
{
   data(PrimarySightings,envir=environment()) #,package="ERAnalysis",envir=environment())
   sightings=PrimarySightings
   PrimarySightings=NULL
}
if(is.null(effort))
{
   data(PrimaryEffort,envir=environment()) #,package="ERAnalysis",envir=environment())
   effort=PrimaryEffort
   PrimaryEffort=NULL
}
sightings$day=as.Date(sightings$Date)-as.Date(paste(sightings$Start.year,"-12-01",sep=""))
sightings$seq=1:nrow(sightings)
if(Use)
{
   sightings=merge(sightings,subset(effort,select=c("key","Use")),by="key")
   sightings=sightings[sightings$Use,]
   sightings=sightings[order(sightings$seq),]
}
# if lcut >0, link sightings
if(lcut>0)
sightings=do.call("rbind",lapply(split(sightings, paste(sightings$Start.year,
        sightings$day, sep = "_")),
         function(x) return(linkdf(x,lcut=lcut,twt=twt,dwt=dwt,crittype=crittype))))
# Next unless the dataframe was provided create the Match dataframe with the 
#   specified parameters
if(is.null(Match))
   Match=create.match(lcut=lcut,mcut=mcut,twt=twt,dwt=dwt,pwt=pwt,crittype=crittype)
# Link Match to effort to obtain Use values if Use==TRUE
# in some cases the watch differs between primary and secondary matched detections
# when the t241 was at a boundary so Use value determined by Primary; also, in some
# cases there was no effort for primary because it was removed because
# there was no effort for the watch that met the beaufort/vis conditions.  In that
# case there is no matching record in primary effort so Use is set to FALSE.
if(Use)
{
   effort$key1=paste(effort$Date,effort$watch,sep="_")
   UseWatch=as.data.frame(table(effort$key1,effort$Use))
   UseWatch=UseWatch[UseWatch$Freq>0,1:2]
   names(UseWatch)=c("key","Use")
   Match=merge(Match,UseWatch,all.x=TRUE)
   Match$Use=factor(Match$Use,levels=c(TRUE,FALSE))
   Match$Use[is.na(Match$Use)]=FALSE
   Match=Match[order(Match$seq),]
   Match$Use=rep(Match$Use[Match$station=="P"],each=2)
   Match=Match[Match$Use==TRUE,]
}
# Run through the set of models for each of the recent years and select the
# best set of models within a delta aic "cutoff"
initial.models=vector("list",length(recent.years))
i=0
for (year in recent.years)
{
   i=i+1
   zz=Match[Match$Start.year==year,]
   zz$Start.year=factor(zz$Start.year)
   zz$Dist=cut(zz$distance,DistBreaks)
   zz$Observer=factor(zz$Observer)
   zz$Vis=cut(zz$vis,c(0,3,6))
   initial.models[[i]]=select.detection.models(zz,models,cutoff)
   print(summary(initial.models[[i]][[1]]))
   cat("\n",length(initial.models[[i]]))
}
# Run through the selected set of models and exclude any in which vis or
# beaufort effects are positive.
formulae=vector("list",length(recent.years))
for (i in 1:length(recent.years))
{
   for (j in 1:length(initial.models[[i]]))
   {
     est=coef(initial.models[[i]][[j]])
     vis.pos=grep("vis",names(est))
     Vis.pos=grep("Vis",names(est))
     beauf.pos=grep("beaufort",names(est))
     if(length(vis.pos)==0 &length(Vis.pos)==0 &length(beauf.pos)==0 )
       formulae[[i]]=c(formulae[[i]],formula(initial.models[[i]][[j]]))
     else
     {
        if(length(vis.pos)!=0 && est[vis.pos]<0)
          formulae[[i]]=c(formulae[[i]],formula(initial.models[[i]][[j]]))
        if(length(Vis.pos)!=0 && est[Vis.pos]<0)
          formulae[[i]]=c(formulae[[i]],formula(initial.models[[i]][[j]]))
        if(length(beauf.pos)!=0 && est[beauf.pos]<0)
          formulae[[i]]=c(formulae[[i]],formula(initial.models[[i]][[j]]))
     }
   }
}
# For each of the 8 recent years, compute the abundance estimate for each of the
# selected detection models unless best==TRUE
detection.models=vector("list",8)
abundance.models=vector("list",8)
#   Exclude effort and sightings for MAS in 1995 because she had no match data for 1995
sightings=sightings[!(sightings$Observer=="MAS"&sightings$Start.year==1995),]
effort=effort[!(effort$Observer=="MAS"&effort$Start.year==1995),]
i=0
if(is.null(gsS))data(gsS,envir=environment()) #,package="ERAnalysis",envir=environment())
for (year in recent.years)
{
   i=i+1
#  Get match data and primary sightings for the year and set up variables
   zz=Match[Match$Start.year==year,]
   zz$Start.year=factor(zz$Start.year)
   primary=sightings[sightings$Start.year==year,]
   primary$Start.year=factor(primary$Start.year)
   primary$Dist=cut(primary$distance,DistBreaks)
   primary$Vis=cut(primary$vis,c(0,3,6))
   zz$Dist=cut(zz$distance,DistBreaks)
   zz$Vis=cut(zz$vis,c(0,3,6))
   zz$Observer=factor(zz$Observer)
   primary$Observer=factor(primary$Observer,levels=levels(zz$Observer))
#  Get effort data depending ion value of Use
   Start.year=NULL
   if(Use)
   {
      ern=subset(effort,subset=Use & as.character(Start.year)==year,
          select=c("Start.year","key","begin","end","effort","time","vis","beaufort"))
   } else
   {
      ern=subset(effort, as.character(Start.year)==year,
          select=c("Start.year","key","begin","end","effort","time","vis","beaufort"))
   }
   ern$Start.year=factor(ern$Start.year)
   if(i==1)
   {
      if(is.null(primary$corrected.podsize))   
         summary.df=data.frame(Year=year,nMatch=nrow(zz)/2,nPrimary=nrow(primary),
                            Whales=sum(primary$podsize),MeanPS=mean(primary$podsize),
                            HrsEffort=sum(ern$effort)*24)
      else
         summary.df=data.frame(Year=year,nMatch=nrow(zz)/2,nPrimary=nrow(primary),
                            Whales=sum(primary$podsize),MeanPS=mean(primary$podsize),CMeanPS=mean(primary$corrected.podsize),
                            HrsEffort=sum(ern$effort)*24)
   }
   else
   { 
      if(is.null(primary$corrected.podsize))   
         summary.df=rbind(summary.df, 
                          data.frame(Year=year,nMatch=nrow(zz)/2, 
                                     nPrimary=nrow(primary),
                                     Whales=sum(primary$podsize),
                                     MeanPS=mean(primary$podsize),
                                     HrsEffort=sum(ern$effort)*24))
      else
        summary.df=rbind(summary.df,
                         data.frame(Year=year,
                                    nMatch=nrow(zz)/2,
                                    nPrimary=nrow(primary),
                                    Whales=sum(primary$podsize),
                                    MeanPS=mean(primary$podsize),
                                    CMeanPS=mean(primary$corrected.podsize),
                                    HrsEffort=sum(ern$effort)*24))
   }
   cat("\n\nSurvey year               : ", year)
   cat("\n# Match records           : ", nrow(zz)/2)   
   cat("\n# Primary sightings       : ", nrow(primary))   
   cat("\n# Total whales            : ", sum(primary$podsize))   
   cat("\nMean observed pod size     : ", mean(primary$podsize))   
   if(!is.null(primary$corrected.podsize))   
   cat("\nMean corrected pod size    : ", mean(primary$corrected.podsize))   
   cat("\nHours of effort           : ", sum(ern$effort)*24)   
#  Loop over each detection model or just the best model
   nmodels=length(formulae[[i]])
   if(best) nmodels=1
   abundance.models[[i]]=vector("list",nmodels)
   detection.models[[i]]=vector("list",nmodels)
   for(j in 1:nmodels)
   {
     dformula=as.formula(paste("~", as.character(formulae[[i]][[j]])[3],sep=""))
#    If TruePS then replace podsize in formula with True which represents the unknown 
#      true podsize which is handled specifically by the code.
     if(TruePS & length(grep("podsize",as.character(dformula)))!=0)
        dtformula=as.formula(paste(sub("podsize","True",as.character(dformula)),collapse=""))
     else
        dtformula=dformula
#    Fit the detection/pod size model
     detection.models[[i]][[j]]=fit.missed.pods(formula=dtformula,pbyyear=TRUE,debug=debug,hessian=hessian,
          data=zz,primary=primary,gsS=gsS)
#    Compute the abundance estimate for this year; depends on whether TruePS=TRUE and
#    whether pod size was corrected to allow use of Reilly approach
     if(is.null(primary$corrected.podsize) | TruePS)
     {
       abundance.models[[i]][[j]]=estimate.abundance(spar=detection.models[[i]][[j]]$par[1:2],
         dpar=detection.models[[i]][[j]]$par[3:length(detection.models[[i]][[j]]$par)],gsS=gsS,effort=ern,
         sightings=primary, final.time=final.time[15+i],lower.time=lower.time[15+i],
         gformula=~s(time),dformula=dtformula)
     }
     else
     {
       abundance.models[[i]][[j]]=estimate.abundance(spar=NULL,
         dpar=detection.models[[i]][[j]]$par[3:length(detection.models[[i]][[j]]$par)],
         gsS=gsS,effort=ern,
         sightings=primary, final.time=final.time[15+i],lower.time=lower.time[15+i],
         gformula=~s(time),dformula=dtformula)
     }
     cat("\nEstimated whale abundance w/o fn: ", abundance.models[[i]][[j]]$Total)
   }
}
# Compute ratio of naive and corrected abundances for 1987-2006 (eqn 23)
Nbest=sapply(abundance.models,function(x) x[[1]]$Total)
Nnaive=sapply(naive.abundance[16:23],function(x)x$Total)
ratio=sum(Nbest)/sum(Nnaive)
# Compute series of estimates for 1967-2006 without nighttime correction factor (eqn 24)  
Nhat=c(sapply(naive.abundance[1:15],
              function(x)x$Total)*ratio,
       sapply(abundance.models,
              function(x) x[[1]]$Total))

# Apply nighttime correction factor (eqn 29)
Nhat=fn*Nhat
summary.df$Nhat=Nhat[16:23]
return(list(summary.df=summary.df,
            Match=Match,
            formulae=formulae,
            detection.models=detection.models,
            abundance.models=abundance.models,
            ratio=ratio,
            Nhat=Nhat,
            TruePS=TruePS))
}

compute.series.var=function(results,naive.abundance,ps.results,Use=TRUE,
                            lcut=0.2,twt=0.18,dwt=3.95,crittype="ellipse",
                            final.time=c(rep(90,20),100,100,90), lower.time=rep(0,23),
                            DistBreaks=c(0,1,2,3,4,20),
                            recent.years=c(1987,1992,1993,1995,1997,2000,2001,2006),
                            fn=1.0817,se.fn=0.0338,debug=FALSE,delta=c(0.001,0.01))
{
# Get PrimarySightings and effort
data(PrimarySightings,envir=environment()) #,package="ERAnalysis",envir=environment())
data(PrimaryEffort,envir=environment()) #,package="ERAnalysis",envir=environment())
data(gsS,envir=environment()) #,package="ERAnalysis",envir=environment())
# Assign day value and match up to Use field if Use=TRUE
PrimarySightings$day=as.Date(PrimarySightings$Date)-as.Date(paste(PrimarySightings$Start.year,"-12-01",sep=""))
PrimarySightings$seq=1:nrow(PrimarySightings)
if(Use)
{
   PrimarySightings=merge(PrimarySightings,subset(PrimaryEffort,
                                                  select=c("key","Use")),by="key")
   PrimarySightings=PrimarySightings[PrimarySightings$Use,]
   PrimarySightings=PrimarySightings[order(PrimarySightings$seq),]
}
# if lcut >0, link PrimarySightings
if(lcut>0)
PrimarySightings=do.call("rbind",lapply(split(PrimarySightings, paste(PrimarySightings$Start.year,
        PrimarySightings$day, sep = "_")),
         function(x) return(linkdf(x,lcut=lcut,twt=twt,dwt=dwt,crittype=crittype))))
#   Exclude effort and sightings for MAS in 1885 because she had no match data for 1995
PrimarySightings=PrimarySightings[!(PrimarySightings$Observer=="MAS"&PrimarySightings$Start.year==1995),]
PrimaryEffort=PrimaryEffort[!(PrimaryEffort$Observer=="MAS"&PrimaryEffort$Start.year==1995),]
#
# For each of the 8 recent years, compute the delta method variance due to variation in abundance estimate for each of the
# selected detection models unless best==TRUE
cat("Computing var1 component\n")
npar=0
for (i in 1:length(recent.years))
  npar=npar+ length(results$detection.models[[i]][[1]]$par)
vc.theta=matrix(0,nrow=npar,ncol=npar)
partial=matrix(0,nrow=23,ncol=npar)
i=0
jpos=0
Ns=sapply(results$abundance.models,function(x) x[[1]]$Total)
NaiveNs.late=sapply(naive.abundance[16:23],function(x)x$Total)
NaiveNs.early=sapply(naive.abundance[1:15],function(x)x$Total)
for (year in recent.years)
{
   i=i+1
#  Get effort data depending on value of Use
   Start.year=NULL
   if(Use)
   {
      ern=subset(PrimaryEffort,subset=Use & as.character(Start.year)==year,
          select=c("Start.year","key","begin","end","effort","time","vis","beaufort"))
   } else
   {
      ern=subset(PrimaryEffort, as.character(Start.year)==year,
          select=c("Start.year","key","begin","end","effort","time","vis","beaufort"))
   }
   ern$Start.year=factor(ern$Start.year)
   primary=PrimarySightings[PrimarySightings$Start.year==year,]
   primary$Start.year=factor(primary$Start.year)
   primary$Dist=cut(primary$distance,DistBreaks)
   primary$Vis=cut(primary$vis,c(0,3,6))
   primary$Observer=factor(primary$Observer,levels(factor(results$Match$Observer[results$Match$Start.year==year])))
   primary$seen=1
#  Loop over each parameter for group size and detection model and 
#  compute first derivatives of abundance estimates
   par.est=results$detection.models[[i]][[1]]$par
   vc.theta[(jpos+1):(jpos+length(par.est)),(jpos+1):(jpos+length(par.est))]=solve(results$detection.models[[i]][[1]]$model$hessian)
   dformula=as.formula(paste("~", as.character(results$formulae[[i]][[1]])[3],sep=""))
   if(results$TruePS & length(grep("podsize",as.character(dformula)))!=0)
      dtformula=as.formula(paste(sub("podsize","True",as.character(dformula)),collapse=""))
   else
      dtformula=dformula
   for (j in 1:length(par.est))
   {
     par.values=par.est
     par.values[j]=par.est[j]*(1-delta[1])
     nlower=estimate.abundance(spar=par.values[1:2],
       dpar=par.values[3:length(par.values)],gsS=gsS,effort=ern,
       sightings=primary, final.time=final.time[15+i],lower.time=lower.time[15+i],
       gformula=~s(time),dformula=dtformula,plotit=FALSE)$Total
     par.values[j]=par.est[j]*(1+delta[1])
     nupper=estimate.abundance(spar=par.values[1:2],
       dpar=par.values[3:length(par.values)],gsS=gsS,effort=ern,
       sightings=primary, final.time=final.time[15+i],lower.time=lower.time[15+i],
       gformula=~s(time),dformula=dtformula,plotit=FALSE)$Total
     partial[i+15,j+jpos]=(nupper-nlower)/(2*delta[1]*par.est[j])
     NewNs=Ns
     NewNs[i]=nlower
     ratio=sum(NewNs)/sum(NaiveNs.late)
     earlyNs.lower=ratio*NaiveNs.early
     NewNs[i]=nupper
     ratio=sum(NewNs)/sum(NaiveNs.late)
     earlyNs.upper=ratio*NaiveNs.early
     partial[1:15,j+jpos]=(earlyNs.upper-earlyNs.lower)/(2*delta[1]*par.est[j])
   }
   jpos=jpos+length(par.est)
}
# This is the first part of the v-c matrix for the N series
var1=partial%*%vc.theta%*%t(partial)
cat("Completed var1 component\n")
# Next compute the v-c matrix due to variation in pod size calibration data
par.est=ps.results$par
vc.theta=ps.results$vc
partial=matrix(0,nrow=23,ncol=length(par.est))
cat("Computing var2 component\n")
for (j in 1:length(par.est))
{
cat("Parameter ",j,"\n")
#    par-delta*par
  par.values=par.est
  par.values[j]=par.est[j]*(1-delta[2])
  nmax=20
  ps.results$par=par.values
  gsS=create.gsS(ps.results,nmax=20)  
  Nlower=rep(0,23)
  i=0
  for (year in recent.years)
  {
    i=i+1
#   Get effort data depending on value of Use
    if(Use)
    {
      ern=subset(PrimaryEffort,subset=Use & as.character(Start.year)==year,
          select=c("Start.year","key","begin","end","effort","time","vis","beaufort"))
    } else
    {
      ern=subset(PrimaryEffort, as.character(Start.year)==year,
          select=c("Start.year","key","begin","end","effort","time","vis","beaufort"))
    }
    ern$Start.year=factor(ern$Start.year)
    zz=results$Match[results$Match$Start.year==year,]
    zz$Start.year=factor(zz$Start.year)
    zz$Dist=cut(zz$distance,DistBreaks)
    zz$Vis=cut(zz$vis,c(0,3,6))
    zz$Observer=factor(zz$Observer)
    primary=PrimarySightings[PrimarySightings$Start.year==year,]
    primary$Start.year=factor(primary$Start.year)
    primary$Dist=cut(primary$distance,DistBreaks)
    primary$Vis=cut(primary$vis,c(0,3,6))
    primary$Observer=factor(primary$Observer,levels(factor(zz$Observer)))
    primary$seen=1
    dformula=as.formula(paste("~", as.character(results$formulae[[i]][[1]])[3],sep=""))
    if(results$TruePS & length(grep("podsize",as.character(dformula)))!=0)
      dtformula=as.formula(paste(sub("podsize","True",as.character(dformula)),collapse=""))
    else
      dtformula=dformula
    dpar=results$detection.models[[i]][[1]]$par
    ddpar=fit.missed.pods(formula=dtformula,pbyyear=TRUE,debug=FALSE,hessian=FALSE,par=dpar,
          data=zz,primary=primary,gsS=gsS)
    if(debug)
    {
       cat("\nconvergence =",ddpar$model$convergence)
       cat("\nvalue =",ddpar$model$value)
       cat("\npar =",ddpar$model$par)
       cat("\ncount =",ddpar$model$counts)
    }
    ddpar=ddpar$par
    spar=ddpar[1:2]
    ddpar=ddpar[3:length(ddpar)]
    Nlower[i+15]=estimate.abundance(spar=spar,dpar=ddpar,gsS=gsS,effort=ern,
       sightings=primary, final.time=final.time[15+i],lower.time=lower.time[15+i],
       gformula=~s(time),dformula=dtformula,plotit=FALSE)$Total
  }
    if(debug)cat("\nNlower=",Nlower)
    par.values=par.est
    par.values[j]=par.est[j]*(1+delta[2])
    ps.results$par=par.values
    gsS=create.gsS(ps.results,nmax=20)  
    Nupper=rep(0,23)
    i=0
    for (year in recent.years)
    {
      i=i+1
#     Get effort data depending on value of Use
      if(Use)
      {
        ern=subset(PrimaryEffort,subset=Use & as.character(Start.year)==year,
          select=c("Start.year","key","begin","end","effort","time","vis","beaufort"))
      } else
      {
        ern=subset(PrimaryEffort, as.character(Start.year)==year,
          select=c("Start.year","key","begin","end","effort","time","vis","beaufort"))
      }
      ern$Start.year=factor(ern$Start.year)
      zz=results$Match[results$Match$Start.year==year,]
      zz$Start.year=factor(zz$Start.year)
      zz$Dist=cut(zz$distance,DistBreaks)
      zz$Vis=cut(zz$vis,c(0,3,6))
      zz$Observer=factor(zz$Observer)
      primary=PrimarySightings[PrimarySightings$Start.year==year,]
      primary$Start.year=factor(primary$Start.year)
      primary$Dist=cut(primary$distance,DistBreaks)
      primary$Vis=cut(primary$vis,c(0,3,6))
      primary$Observer=factor(primary$Observer,levels(factor(zz$Observer)))
      primary$seen=1
      dformula=as.formula(paste("~", as.character(results$formulae[[i]][[1]])[3],sep=""))
      if(results$TruePS & length(grep("podsize",as.character(dformula)))!=0)
        dtformula=as.formula(paste(sub("podsize","True",as.character(dformula)),collapse=""))
      else
        dtformula=dformula
      dpar=results$detection.models[[i]][[1]]$par
      ddpar=fit.missed.pods(formula=dtformula,pbyyear=TRUE,debug=FALSE,hessian=FALSE,par=dpar,
          data=zz,primary=primary,gsS=gsS)
      if(debug)
      {
         cat("\nconvergence =",ddpar$model$convergence)
         cat("\nvalue =",ddpar$model$value)
         cat("\npar =",ddpar$model$par)
         cat("\ncount =",ddpar$model$counts)
      }
      ddpar=ddpar$par
      spar=ddpar[1:2]
      ddpar=ddpar[3:length(ddpar)]
      Nupper[i+15]=estimate.abundance(spar=spar,dpar=ddpar,gsS=gsS,effort=ern,
         sightings=primary, final.time=final.time[15+i],lower.time=lower.time[15+i],
         gformula=~s(time),dformula=dtformula,plotit=FALSE)$Total
    }
    if(debug) cat("\nNupper=",Nupper)
    partial[16:23,j]=(Nupper[16:23]-Nlower[16:23])/(2*delta[2]*par.est[j])
    ratio=sum(Nlower[16:23])/sum(NaiveNs.late)
    earlyNs.lower=ratio*NaiveNs.early
    ratio=sum(Nupper[16:23])/sum(NaiveNs.late)
    earlyNs.upper=ratio*NaiveNs.early
    partial[1:15,j]=(earlyNs.upper-earlyNs.lower)/(2*delta[2]*par.est[j])
  }
if(debug)write.table(partial,"partial2.txt")
var2=partial%*%vc.theta%*%t(partial)
cat("Completed var2 component\n")
# V-c matrix component from fitting gam and residual variance around gam
cat("Computing var3 component\n")
var3=matrix(0,nrow=23,ncol=23)
var.Nhat=sapply(results$abundance.models,function(x) x[[1]]$var.Total)
diag(var3[16:23,16:23])=var.Nhat
ratio=sum(Ns)/sum(NaiveNs.late)
sigma.ratio=(var(Ns)+ratio^2*var(NaiveNs.late)-2*ratio*cov(Ns,NaiveNs.late))/(length(recent.years)*mean(NaiveNs.late)^2)
var.NaiveNs=sapply(naive.abundance, function(x) x$var.Total)
var.NaiveNs.early=var.NaiveNs[1:15]
var.NaiveNs.late=var.NaiveNs[16:23]
var3[1:15,1:15]=sigma.ratio*outer(NaiveNs.early,NaiveNs.early,"*")
diag(var3[1:15,1:15])=ratio^2*NaiveNs.early^2*
              (sigma.ratio*(1+length(recent.years))/ratio^2+var.NaiveNs.early/NaiveNs.early^2)
covar.ratio=outer(NaiveNs.early,(var.Nhat/NaiveNs.late-Ns*ratio*var.NaiveNs.late/NaiveNs.late^2),"*")
var3[1:15,16:23]=covar.ratio
var3[16:23,1:15]=t(covar.ratio)
vc=var1+var2+var3
se=sqrt(diag(vc))
# Add adjustment for nighttime correction factor; first line does covariances f^2*cov(Nhat_i,Nhat_j)
vc=vc+fn^2*vc
# Following adjusts variances
diag(vc)=(results$Nhat)^2*((se.fn/fn)^2+(se/results$Nhat)^2)
se=sqrt(diag(vc))
return(list(vc=vc,var1=var1,var2=var2,var3=var3,se=se,cormat=vc/outer(se,se,"*")))
}

create.match=function(lcut=0.2,mcut=1,twt=.18,dwt=3.95,pwt=0.05,crittype="ellipse",fdir="c:/gw")
{
  Ld2Ll  <- function(Ld) 
  {
     Ll = split(Ld, paste(Ld$Start.year,Ld$day, sep = "_"))
     return(Ll)
  }
# Extract match data and call algorithm to make links/matches
  Match.list<- extract.match.data()
  if(lcut>0)
  {
     Linked.data=do.call("rbind",lapply(Match.list,
         function(x) return(linklist(x,lcut=lcut,twt=twt,dwt=dwt,crittype=crittype))))
     Linked.list <- Ld2Ll(Linked.data)
     Matched.data=do.call("rbind",lapply(Linked.list,function(x) 
         return(matchdf(x,mcut=mcut,twt=twt,dwt=dwt,pwt=pwt,crittype=crittype))))
  }
  else
  {
     Matched.data=do.call("rbind",lapply(Match.list,function(x) 
         return(matchdf(x,mcut=mcut,twt=twt,dwt=dwt,pwt=pwt,crittype=crittype))))

  }
# Rename and select needed fields and create key field
  Matched.data$station="P"
  Matched.data$station[Matched.data$PorS==2]="S"
  Matched.data$station=factor(Matched.data$station) 
  Match=NULL
  AllPrimaryEffort=NULL
  AllSecondaryEffort=NULL
  data(Match,envir=environment()) #,package="ERAnalysis",envir=environment())
  match.names=names(Match)
  xx=subset(Matched.data,select=c(names(Matched.data)[names(Matched.data)%in%match.names],"observer"))
  xx$Date=as.character(strptime(paste(xx$Start.year,"-12-01",sep=""),format="%Y-%m-%d")+3600*24*(xx$day-1))
  xx$key=paste(xx$Date,xx$watch,sep="_")
# get observer, beaufort, vis
# first do Primary
  data(AllPrimaryEffort,envir=environment()) #,package="ERAnalysis",envir=environment())
  data(AllSecondaryEffort,envir=environment()) #,package="ERAnalysis",envir=environment())
  PrimaryEffort=AllPrimaryEffort
  xx$seq=1:nrow(xx)
  zz=merge(xx[xx$station=="P",],PrimaryEffort,by="Date")
  zz$begin.time=24*(zz$begin-floor(zz$begin))
  zz$end.time=24*(zz$end-floor(zz$end))
  zz=zz[zz$t241<=(zz$end.time+1/3600)&zz$t241>=zz$begin.time,]
  zz=zz[order(zz$seq),]
# next do Secondary
  SecondaryEffort=AllSecondaryEffort
  zz.sec=merge(xx[xx$station=="S",],SecondaryEffort,by="Date")
  zz.sec$begin.time=24*(zz.sec$begin-floor(zz.sec$begin))
  zz.sec$end.time=24*(zz.sec$end-floor(zz.sec$end))
  zz.sec=zz.sec[zz.sec$t241<=(zz.sec$end.time+1/3600)&zz.sec$t241>=zz.sec$begin.time,]
  zz.sec=zz.sec[order(zz.sec$seq),]
# create list of records in which t241 of both sightings are in effort period
  on.effort=(xx$seq[xx$station=="P"] %in% zz$seq) & (xx$seq[xx$station=="S"] %in% zz.sec$seq)
  xx=xx[rep(on.effort,each=2),]
  zz=zz[zz$seq %in% xx$seq,]
  zz.sec=zz.sec[zz.sec$seq %in% xx$seq,]
# extract observer record from effort record containing t241;  this changes observer
# on record of about 2-3% of the sightings
  xx$observer[xx$station=="P"]=as.character(zz$Observer)
  xx$observer[xx$station=="S"]=as.character(zz.sec$Observer)
  xx$observer[xx$observer=="35"]="DJR"
# add beaufort and vis
  xx$beaufort[xx$station=="P"]=zz$beaufort.y
  xx$beaufort[xx$station=="S"]=zz.sec$beaufort.y
  xx$vis[xx$station=="P"]=zz$vis.y
  xx$vis[xx$station=="S"]=zz.sec$vis.y
# create observer experience table and add hours of experience to match data
  ObserverExp=CreateObserverExperience()
  xx=merge(xx,ObserverExp,by="Date",all.x=TRUE)
  xx$observer=as.character(xx$observer)
  xx$hours=xx[cbind(1:nrow(xx),match(xx$observer,names(xx)[16:74])+15)]
  xx=xx[,c(1:15,75)]
  xx$hours=as.numeric(xx$hours)
# next construct pods per hour (as determined from primary observer) during the watch
# containing the match data
  PrimaryEffort$key=paste(PrimaryEffort$Date,PrimaryEffort$watch,sep="_")
  pphr=with(PrimaryEffort, tapply(npods,key,sum))
  pphr=pphr/with(PrimaryEffort, tapply(effort,key,sum))
  xx=merge(xx,data.frame(key=names(pphr),pphr=pphr),by="key")
# Finally add on the sex of the observer
  data(Observer,envir=environment()) #,package="ERAnalysis",envir=environment())
  Observer$key=as.character(Observer$Initials)
  Observer$key[is.na(Observer$Initials)]=Observer$Observer[is.na(Observer$Initials)]
  xx=merge(xx, subset(Observer, select = c("key", "Sex")), by.x ="observer",by.y = "key", all.x = TRUE)
  xx$Observer=factor(xx$observer)
  xx$observer=NULL
  xx=xx[order(xx$seq),]
  row.names(xx)=NULL
  return(xx)
}
CreateObserverExperience=function()
{
# Create dataframe of hours of experience for each observer
# First extract values from early and recent data
EarlyEffort=NULL
ERSurveyData=NULL
data(EarlyEffort,envir=environment()) #,package="ERAnalysis",envir=environment())
EarlyExp=data.frame(Observer=EarlyEffort$Observer,Date=substr(EarlyEffort$key,1,10),Hours=as.vector(EarlyEffort$End.date.time-EarlyEffort$Begin.date.time))
data(ERSurveyData,envir=environment()) #,package="ERAnalysis",envir=environment())
EFLAG=NULL
EXPERIMENT=NULL
Start.watch=subset(ERSurveyData,subset=EFLAG==1&EXPERIMENT%in%c(1,2),select=c("OBSERVER","DATE","ETIME"))
End.watch=subset(ERSurveyData,subset=EFLAG==5&EXPERIMENT%in%c(1,2),select=c("OBSERVER","DATE","ETIME"))
Start.watch=Start.watch[order(Start.watch$DATE,Start.watch$ETIME),]
End.watch=End.watch[order(End.watch$DATE,End.watch$ETIME),]
RecentExp=Start.watch
RecentExp$Hours=End.watch$ETIME-Start.watch$ETIME
RecentExp=subset(RecentExp,select=c("OBSERVER","DATE","Hours"))
names(RecentExp)[1:2]=c("Observer","Date")
EarlyExp$Observer=as.character(EarlyExp$Observer)
RecentExp$Observer=as.character(RecentExp$Observer)
RecentExp$Date=as.character(RecentExp$Date)
EarlyExp$Date=as.character(EarlyExp$Date)
# Merge data from each set of surveys
ObserverExp=rbind(EarlyExp,RecentExp)
ObserverExp$Observer=factor(ObserverExp$Observer)
ObserverExp$Date=factor(ObserverExp$Date)
# Apply to get sum of effort by observer and date
ObserverExp=t(tapply(ObserverExp$Hours,list(ObserverExp$Observer,ObserverExp$Date),sum))
ObserverExp[is.na(ObserverExp)]=0
# Create cumsum over dates
ObserverExp=apply(ObserverExp,2,cumsum)
# Create a 0 row for first date and then offset dates such that for
# day x the experience is through day x-1
odates=row.names(ObserverExp)
ObserverExp=rbind(rep(0,ncol(ObserverExp)),ObserverExp)
ObserverExp=ObserverExp[-nrow(ObserverExp),]
ObserverExp=as.data.frame(ObserverExp)
ObserverExp$Date=odates
row.names(ObserverExp)=NULL
# Next sum columns for DJR with 35 and CDV with 33 because they had different codes
# amongst the datasets
ObserverExp$DJR=ObserverExp$DJR+ObserverExp[,"35"]
ObserverExp$CDV=ObserverExp$CDV+ObserverExp[,"33"]
return(ObserverExp)
}


estimate.abundance <- function(spar, dpar, gsS, effort, sightings, 
                               dformula=~True, gformula=~s(time),nmax=20,
                               pod=FALSE,plotit=TRUE,anchor=TRUE,
                               show.anchor=FALSE,sp=NULL,final.time=90,
                               lower.time=0,do.mult=TRUE,pool=TRUE,...)
{
# Function to compute observation-specific detection probability
# from formula parameters and sightings data
probs=function(j)
{
  sightings$True=j
  dmat=model.matrix(dformula,sightings)
  plogis(dmat%*%dpar)
}
# Order sightings by Start.year and create list of years
sightings=sightings[order(sightings$Start.year),]
years=levels(factor(sightings$Start.year))
Nhat.whales=NULL
Nhat=NULL
for(year in years)
{
   cat("\n\nSurvey year               : ", year)
   cat("\n# Primary sightings       : ", nrow(sightings[sightings$Start.year==year,]))   
   cat("\n# Total whales            : ", sum(sightings$podsize[sightings$Start.year==year]))   
   cat("\nMean observed podsize     : ", mean(sightings$podsize[sightings$Start.year==year]))   
   cat("\nHours of effort           : ", sum(effort$effort)*24)   
}
# Depending on values of dpar and spar construct estimates of
# number of whales and number of pods passing for each observation
if(!is.null(dpar)&!is.null(dformula)&!is.null(spar))
{
   # Compute matrix of detection probabilities for each observation (rows) at each
   # of 1 to nmax possible true pod sizes
   ps=sapply(1:nmax, probs)
   # Corrections for pod size and detection probability
   # Compute probability function for True pod size
   i=0                  
   for(year in years)
   {
      i=i+1
      fS=gammad(spar[((i-1)*2+1):(i*2)],nmax)
      # Compute conditional detection probability function for true size given observed size
      fSs=t(t(fS*gsS)/colSums(fS*gsS))
      # Compute estimate of expected number of whales represented by observed whales
      Nhat.whales=c(Nhat.whales, 
                    rowSums(t(fSs[,sightings$podsize[sightings$Start.year==year]]*(1:nmax))/ps[sightings$Start.year==year]))
      # Same as above for pods rather than whales
      Nhat=c(Nhat,rowSums(t(fSs[,sightings$podsize[sightings$Start.year==year]])/ps[sightings$Start.year==year]))
   }
} else
{
  # Corrections for detection probability but
  # no further corrections for pod size but pod size
  # may have been corrected prior as in reilly.cf
  if(!is.null(dpar)&!is.null(dformula))
  {
    dmat=model.matrix(dformula,sightings)
    ps=plogis(dmat%*%dpar)
    if(is.null(sightings$corrected.podsize))
       Nhat.whales=sightings$podsize/ps
      else
       Nhat.whales=sightings$corrected.podsize/ps
    Nhat=1/ps
  }else
  {
    # Corrections for pod size but no corrections
    # detection probability
    if(!is.null(spar))
    {
       i=0
      for(year in years)
      {
         i=i+1                                                       
         # Compute probability function for True pod size
         fS=gammad(spar[((i-1)*2+1):(i*2)],nmax)
         # Compute conditional detection probability function for true size given observed size
         fSs=t(t(fS*gsS)/colSums(fS*gsS))
         # Compute estimate of expected number of whales represented by observed whales
         Nhat.whales=c(Nhat.whales,rowSums(t(fSs[,sightings$podsize[sightings$Start.year==year]]*(1:nmax))))
      }
      Nhat=rep(1,nrow(sightings))
    } else
    # No further corrections for either but pod size
    # may have been corrected prior as in reilly.cf
    {
      if(is.null(sightings$corrected.podsize))
         Nhat.whales=sightings$podsize
      else
         Nhat.whales=sightings$corrected.podsize
      Nhat=rep(1,nrow(sightings))
    }
  }
}
#
# Merge estimates with effort data to construct migration timing model
# to extrapolate from sampled periods to entire migration time period
#
sightings$Nhat.whales=Nhat.whales
sightings$Nhat=Nhat
if(!pod)
  est.df=data.frame(key=unique(sightings$key),
          nhat=sapply(split(sightings$Nhat.whales,factor(sightings$key)),sum))
else
  est.df=data.frame(key=unique(sightings$key),
          nhat=sapply(split(sightings$Nhat,factor(sightings$key)),sum))
ern=subset(effort,select=c("Start.year","key","begin","end","effort","time","vis","beaufort"))
ern=merge(ern,est.df,by.x="key",by.y="key",all.x=TRUE)
ern$nhat[is.na(ern$nhat)]=0
results=fit.migration.gam(ern, years=as.numeric(years), formula=gformula, pod=pod, 
                          plotit=plotit,
                          anchor=anchor,show.anchor=show.anchor,sp=sp,
                          final.time=final.time,lower.time=lower.time,
                          do.mult=do.mult,pool=pool,...)

results$options=list(anchor=anchor,pod=pod,final.time=final.time,lower.time=lower.time,
                     gformula=gformula,dformula=dformula,spar=spar,dpar=dpar,nmax=nmax)
cat("\nEstimated whale abundance : ",results$Total)
return(results)
}


extract.match.data <- function ()
{
# Extracts the relevant data from recent surveys to do matching of double observer data.
# There are no arguments but it uses the dataframe ERSruveyData and Height which
# contains the platform heights at the observation locations each year.
# Returns a list of dataframes that are divided based on survey date.
#
# added ret, angle, dist to sighting, dist to Azmuth
#
    DistAzmuth <- function (sangle, sdist){return((sin(pi * (sangle - 241)/180) * sdist))}
    data(ERSurveyData,envir=environment()) #,package="ERAnalysis",envir=environment())
    ERSurveyData$LOCATION[ERSurveyData$LOCATION=="NP"]="N"
    Height=NULL
    data(Height,envir=environment()) #,package="ERAnalysis",envir=environment())
# Extract primary observer sightings in which travel direction is not North
# Primary observers are at South location except in years 2000,2001
# This could be used to create Primary dataframe that is in package.  Currently
# Primary and SecondarySightings are constructed from ERAbund files.
    xx = ERSurveyData[ERSurveyData$EFLAG == 3 & ERSurveyData$TRAVELDIR !=
        "N", ]
    xx = rbind(xx[xx$Start.year %in% c(2000, 2001) & (xx$EXPERIMENT ==
        1 | (xx$EXPERIMENT == 2 & xx$LOCATION == "N")), ], xx[!(xx$Start.year %in%
        c(2000, 2001)) & (xx$EXPERIMENT == 1 | (xx$EXPERIMENT ==
        2 & xx$LOCATION == "S")), ])
    Primaryx = xx
    Primaryx$station = "P"
# Extract secondary observer sightings in which travel direction is not North
# Secondary observers are at North location except in years 2000,2001
    xx = ERSurveyData[ERSurveyData$EFLAG == 3 & ERSurveyData$TRAVELDIR !=
        "N", ]
    xx = rbind(xx[xx$Start.year %in% c(2000, 2001) & (xx$EXPERIMENT ==
        2 & xx$LOCATION == "S"), ], xx[!(xx$Start.year %in% c(2000,
        2001)) & (xx$EXPERIMENT == 2 & xx$LOCATION != "S"), ])
    Secondaryx = xx
    Secondaryx$station = "S"
# Create Match.data dataframe with the experiment 2 records
    Match.data = rbind(Primaryx[Primaryx$EXPERIMENT == 2, ],
        Secondaryx[Secondaryx$EXPERIMENT == 2, ])
# get reticles,angles and times using South sightings if available otherwise use North
    reticles = Match.data$SRET
    reticles[floor(reticles * 1000) == 0] = Match.data$NRET[floor(reticles *
                                                                    1000) == 0]
    reticles[floor(reticles * 1000) == 0] = NA
    angles = Match.data$SANGLE
    angles[floor(angles) == 0] = Match.data$NANGLE[floor(angles) == 0]
    times = Match.data$STIME
    times[floor(times) == 0] = Match.data$NTIME[floor(times) ==  0]
    # constuct key field to match into Height dataframe to get heights for
    # distance calculation from reticles
    Match.data$key = factor(paste(Match.data$Start.year,Match.data$LOCATION,sep = ""))
    
    #print(names(Height))
    
    heights = merge(Match.data, Height, all.x = TRUE, by = "key")$Height
    distances = RetDistK(Height = heights, Reticles = reticles)
    # convert radial distance to perpendicular distance
    Match.data$sdup = RetDistK(Height = heights, Reticles = (reticles-.05))
    Match.data$sddn = RetDistK(Height = heights, Reticles = (reticles+.05))
    Match.data$sdist = distances
    Match.data$sang2A = angles-241
    Match.data$stime = times
    Match.data$sret = reticles
    Match.data$distance = D241(angles, distances)
    Match.data$sd2A = DistAzmuth(angles, distances)
    # Compute crossing time at the line perpendicular to the coast at the shore station (abeam)
    Match.data$t241 = T241H(times, angles, distances, Speed = 6)
    # Compute day since 1 Dec
    Match.data$day = as.vector(Match.data$DATE - as.POSIXct(paste(Match.data$Start.year,
                                                                  "-12-01", sep = ""))+1)
    # Compute watch factor variable
    Match.data$watch = 1
    Match.data$watch[Match.data$t241 >= 10.5] = 2
    Match.data$watch[Match.data$t241 >= 13.5] = 3
    # Create dataframe excluding those with missing T241 or distance and modify field names
    t241=NULL
    distance=NULL
    Match.data = subset(Match.data, subset=(t241 & distance),
                        select = c("station", "day",  "watch", "t241", "distance"
                                   , "PODSIZE","stime", "sret", "sdist", "sang2A", "sd2A","sdup","sddn"
                                   , "OBSERVER", "VISCODE", "WINDFORCE", "WINDDIR"
                                   , "Start.year"))
    names(Match.data)[names(Match.data) == "PODSIZE"] = "podsize"
    names(Match.data)[names(Match.data) == "OBSERVER"] = "observer"
    names(Match.data)[names(Match.data) == "VISCODE"] = "vis"
    names(Match.data)[names(Match.data) == "WINDFORCE"] = "beaufort"
    names(Match.data)[names(Match.data) == "WINDDIR"] = "wind.direction"
    # The following will spit the Match by DATE to do the matching.
    Match.list = split(Match.data, paste(Match.data$Start.year,
                                         Match.data$day, sep = "_"))
    rm(ERSurveyData)
    return(Match.list)
}


fit.fS <-
function(par,pstable,gsS,nmax=20)
{
# Function to fit observed pod size distribution and estimate true
# pod size distribution based on gamma distribution

  ps=gammad(par,nmax)
  -sum(pstable*log(colSums(gsS*ps)))
}

fit.gamma.re=function(x,re.factor,sformula,rformula,sigma.formula=~1,initial,maxit=1000,method="BFGS",hessian=FALSE,nmax=20,z=5)
{
# Fits a gamma mixed-effects model to the pod size calibration data.
#
# Arguments:
#
#  x             - dataframe
#  re.factor     - list of one or more vectors of factor variables to split x into
#                  a list of dataframes. Length of each factor variable vector must
#                  match number of rows in x.
#  sformula      - formula for shape of gamma
#  rformula      - formula for rate of gamma
#  sigma.formula - formula for sigma of normal error
#  initial       - required vector of initial parameter values
#  maxit         - maximum number of iterations (passed to optim)
#  method        - method used by optim for optimization
#  hessian       - passed to optim, if TRUE returns hessian for v-c matrix
#  nmax         -  maximum pod size
#
#  Value: list returned with optim results
#
   mm=lapply(split(x,re.factor,drop=TRUE),function(x) model.matrix(rformula,x))
   rr=lapply(split(x,re.factor,drop=TRUE),function(x) model.matrix(sformula,x))
   ss=lapply(split(x,re.factor,drop=TRUE),function(x) model.matrix(sigma.formula,x))
   observed.size=split(x$Estimate,re.factor,drop=TRUE)
   return(optim(initial,regam.lnl,method=method,control=list(maxit=maxit),
                hessian=hessian,mm=mm,rr=rr,ss=ss,observed.size=observed.size,nmax=nmax,z=z))
}

regam.lnl=function(par,mm,rr,ss,observed.size,nmax=20,z=5)
{
#
# Compute negative log-likelihood for a random effects gamma model for
# the pod size calibration data. A normal 0,sigma error distribution is used
# for the random effect.  It is used for the log(rate) where rate is a
# gamma parameter.
#
#  Arguments:
#   par           -  parameter values
#   mm            -  list of model matrices for rate; list elements are split by random effect factor levels
#   rr            -  list of model matrices for shape; list elements are split by random effect factor levels
#   ss            -  list of model matrices for sigma; list elements are split by random effect factor levels
#   observed.size -  list of vectors of observed sizes; list elements are split by random effect factor levels
#   z             -  defines integration bounds for normal (default is 5 for integration range of -5*sigma,5*sigma)
#   nmax          -  maximum pod size
#
#  Note that the mm,rr,ss, observed.size lists need to be defined to split the data up
#  by the intersection of the factors used to define mm,rr and ss such that for each grouping
#  there is only a single sigma.
#
#  Value:  negative log-likelihood
#
  gam.reint=function(x,obs,xbeta,shape,sigma,nmax=20)
  {
#   Internal function to be used by integrate function for integration over normal random effect
#
#   Arguments:
#
#   x    - vector of values passed by integrate
#   obs  - observed size values
#   xbeta- x%*%beta where beta is parameter vector for rate of gamma
#   shape- shape parameters for gamma
#   sigma- value of sigma for normal
#   nmax - maximum pod size
#
    rate=exp(outer(xbeta,x,"+"))
    return(dnorm(x,sd=sigma)*apply(rate,2,function(rate) prod(gam.d(shape,rate,x=obs,nmax=nmax))))
  }
# The parameter vector contains sigma parameters, rate parameters and then shape parameters
# The number of columns of the design matrix ss determines number of sigma parameters
  npars=ncol(ss[[1]])
  nparr=ncol(mm[[1]])
  neglnl=0
  cat("\npar = ",par)
# Loop over each random effect
  for(i in 1:length(mm))
  {
#    Compute xbeta (x%*%beta) where x is design matrix and beta are parameters for rate
     xbeta=as.vector(mm[[i]]%*%par[(npars+1):(npars+nparr)])
#    Compute sigma
     sigma=unique(exp(as.vector(ss[[i]]%*%par[1:npars])))
#    Compute shape
     shape=unique(exp(as.vector(rr[[i]]%*%par[(npars+nparr+1):length(par)])))
     if(length(sigma)>1)stop("\nInvalid sigma values\n")
     obs=observed.size[[i]]
#    Compute negative log likelihood by integrating this portion (subset of data defined by
#    random effects) of the likelihood over the normal distribution for the random effect
#    and accumulate the value.
     neglnl=neglnl-log(integrate(gam.reint,-z*sigma,z*sigma,obs=obs,xbeta=xbeta,sigma=sigma,shape=shape,nmax=nmax,stop.on.error=FALSE)$value)
  }
  cat("\npar = ",par)
  cat("\n neglnl= ",neglnl)
  return(neglnl)
}

fit.gamma=function(x,sformula,rformula,initial,maxit=1000,method="BFGS",hessian=FALSE,nmax=20)
{
# Fits a gamma fixed-effects model to the pod size calibration data.
#
# Arguments:
#
#  x             - dataframe
#  sformula      - formula for shape of gamma
#  rformula      - formula for rate of gamma
#  initial       - required vector of initial parameter values
#  maxit         - maximum number of iterations (passed to optim)
#  method        - method used by optim for optimization
#  hessian       - passed to optim, if TRUE returns hessian for v-c matrix
#  nmax         -  maximum pod size
#
#  Value: list returned with optim results
#
   smat=model.matrix(sformula,x)
   rmat=model.matrix(rformula,x)
   return(optim(initial,gam.lnl,method=method,control=list(maxit=maxit),
                hessian=hessian,rmat=rmat,smat=smat,observed.size=x$Estimate,nmax=nmax))
}

gam.lnl=function(par,smat,rmat,observed.size,nmax=20)
{
# Compute negative log-likelihood for a fixed effects gamma model for
# the pod size calibration data.
#
#  Arguments:
#   par           -  parameter values
#   smat          -  model matrix for shape
#   rmat          -  model matrix for rate
#   observed.size -  vector of observed sizes
#   nmax          -  maximum pod size
#
#  Value:  negative log-likelihood
#
  shape=exp(smat%*%par[1:ncol(smat)])
  rate=exp(rmat%*%par[(ncol(smat)+1):length(par)])
  neglnl=-sum(log(gam.d(shape,rate,x=observed.size,nmax=nmax)))
  cat("\npar = ",par)
  cat("\n neglnl= ",neglnl)
return(neglnl)
}

gam.d=function(shape,rate,x,nmax)
{
# Computes the gamma probabilities for a vector of shape and rates matched to observed values (x).
# Because these are pod size estimates (1,2,...) the integral for x is from x-1 to x of the gamma.
# Also, to avoid dealing with an infinite number of possible values the maximum pod size
# is set (nmax) and the distribution is renormalized over the range 0 to nmax.
#
# Arguments:
#   shape        - vector of shapes
#   rate         - vector of rates
#   x            - vector of observed pod sizes
#   nmax         -  maximum pod size
#
   if(length(x)!=length(rate))stop("\n****invalid lengths****/n")
   num=apply(matrix(pgamma(c(x-1,x),shape=shape,rate=rep(rate,2)),ncol=2),1,diff)
   denom=pgamma(rep(nmax,length(x)),shape=shape,rate=rate)
#   denom[num<1e-16]=1
   return(num/denom)
}

fit.migration.gam=function(er.migdata, years, formula=~s(time), pod=FALSE, plotit=TRUE,
                           anchor=TRUE,
                           show.anchor=FALSE,sp=NULL,final.time=rep(90,length(years)),
                           lower.time=rep(0,length(years)),do.mult=TRUE,pool=TRUE,...)
{
   final=final.time
   lower=lower.time
   if(length(final)!=length(years))stop("Number of elements in final.time does not match number of years")
   if(length(lower)!=length(years))stop("Number of elements in lower.time does not match number of years")
   require(mgcv)
   ermod=vector("list",length(years))
   formula=as.formula(paste("nhat", paste(as.character(formula),collapse=""),sep=""))
   if(pod)
     ylabel="Pods per day"
   else
     ylabel="Whales per day"
   if(!pool)
   {
   i=0
   for (year in years)
   {
     i=i+1
     er=er.migdata[er.migdata$Start.year==year,]
     er$offset=log(er$effort)
     if(anchor)
     {
        if(lower[i] <= (er$begin[1]-1))
        {
           er=rbind(er[1,],er)
           er$time[1]=c(lower[i]+0.5)
           er$begin[1]=c(lower[i])
           er$end[1]=c(lower[i]+1)
           er$offset[1]=c(0)
           er$effort[1]=1
           er$vis[1]=1
           er$beaufort[1]=0           
           er$nhat[1]=0
        }
        else
           lower[i]=floor(er$begin[1])
        if(er$end[nrow(er)]<(final[i]-1))
        {
           er=rbind(er,er[1,])
           er$time[nrow(er)]=final[i]-0.5
           er$begin[nrow(er)]=final[i]-1
           er$end[nrow(er)]=final[i]
           er$offset[nrow(er)]=0
           er$effort[nrow(er)]=1
           er$nhat[nrow(er)]=0
           er$vis[nrow(er)]=1
           er$beaufort[nrow(er)]=0           

        }
        else
           final[i]=ceiling(er$end[nrow(er)])
     }
     ermod[[i]]=gam(formula,data=er,offset=offset,family=quasipoisson)
     ppd = tapply(er$nhat,
                  floor(er$time), 
                  sum)/tapply(er$effort, floor(er$time),sum)
     
     if(plotit)
     {
        Eppd=predict(ermod[[i]],type="response")
        ymax=max(c(Eppd,ppd))*1.05
        plot(er$time, Eppd, ylim=c(0,ymax),type="l",
             main=paste(year,"/",year+1,sep=""),
             xlab="Days since 1 Dec", 
             ylab=ylabel, xlim=c(0,100))
        points(as.numeric(names(ppd)), ppd)
        
        # add ggplot and save them:
        library(ggplot2)
        pred.plot <- ggplot(data.frame(time = er$time,
                                       gam.fit = Eppd)) +
          geom_line(aes(x = time, y = gam.fit)) +
          geom_point(data = data.frame(date = as.numeric(names(ppd)),
                                       ppd = ppd),
                     aes(x = date,
                         y = ppd)) + 
          xlab("Days since 1 Dec") +
          ylab("Whales per day") +
          labs(title = paste0(year,"/",year+1))
               
        ggsave(plot = pred.plot,
               filename = paste0("figures/Laake_", 
                                 paste0(year,"-",year+1), ".png"),
               device = "png",
               dpi = 600)
        
        if(anchor & show.anchor)
        {
           if(lower.time[i] <= (er$begin[1]-1)) er=er[-1,]
           if(er$end[nrow(er)-1]<(final.time[i]-1))er=er[-nrow(er),]
           mod=gam(formula,data=er,offset=offset,family=quasipoisson)
           newdata=data.frame(time=(floor(lower[i]):(final[i]-1))+.5)
           newdata$Start.year=factor(rep(year,nrow(newdata)),levels=years)
           lines((floor(lower[i]):(final[i]-1))+.5,
             predict(mod,newdata=newdata,type="response"),
             xlim=c(lower[i],final[i]),lty=2)
        }
     }
   }
   Total=vector("numeric",length(years))
   for(i in 1:length(years))
   {
        newdata=data.frame(time=(lower[i]:(final[i]-1))+.5)
        newdata$Start.year=factor(rep(years[i],nrow(newdata)),levels=years)
        if(length(grep("vis",as.character(formula)))!=0)
          newdata$vis=1
        if(length(grep("beaufort",as.character(formula)))!=0)
          newdata$beaufort=0
        Total[i]=sum(predict(ermod[[i]],newdata=newdata,type="response"))
   }
   if(do.mult)
   {
      TotalM=vector("numeric",length(years))
      mult=vector("numeric",length(years))
      i=0
      for(year in years)
      {
        i=i+1
        mult[i]=compute.sampling.multiplier(ermod[[i]], 
                                            er.migdata[er.migdata$Start.year==year,],
                upper=final[i],lower=lower[i])
        TotalM[i]=sum(ermod[[i]]$y)*mult[i]
      }
      return(list(models=ermod,Total=Total,TotalM=TotalM,mult=mult))
   }
   else
      return(list(models=ermod,Total=Total))
   }
   else
   {
   i=0
   ern=NULL
   for (year in years)
   {
     i=i+1
     er=er.migdata[er.migdata$Start.year==year,]
     er$offset=log(er$effort)
     if(anchor)
     {
        if(lower[i] <= (er$begin[1]-1))
        {
           er=rbind(er[1,],er)
           er$time[1]=c(lower[i]+0.5)
           er$begin[1]=c(lower[i])
           er$end[1]=c(lower[i]+1)
           er$offset[1]=c(0)
           er$effort[1]=1
           er$nhat[1]=0
           er$vis[1]=1
           er$beaufort[1]=0                      
        }
        else
           lower[i]=floor(er$begin[1])
        if(er$end[nrow(er)]<(final[i]-1))
        {
           er=rbind(er,er[1,])
           er$time[nrow(er)]=final[i]-0.5
           er$begin[nrow(er)]=final[i]-1
           er$end[nrow(er)]=final[i]
           er$offset[nrow(er)]=0
           er$effort[nrow(er)]=1
           er$nhat[nrow(er)]=0
           er$vis[nrow(er)]=1
           er$beaufort[nrow(er)]=0           
        }
        else
           final[i]=ceiling(er$end[nrow(er)])
     }
     ern=rbind(ern,er)
     }
     ern$Start.year=factor(ern$Start.year)
     if(!is.null(sp))
     {
        ermod=gam(formula,data=ern,offset=offset,family=quasipoisson,sp=sp,...)
        ermod=gam(formula,data=ern,offset=offset,family=quasipoisson,start=ermod$coef,...)
     } else
        ermod=gam(formula,data=ern,offset=offset,family=quasipoisson,...)     
     if(anchor & show.anchor)
     {
       er=er.migdata
       er$offset=log(er$effort)
       er$Start.year=factor(er$Start.year)
       ermod.na=gam(formula,data=er,offset=offset,family=quasipoisson,...)
       Eppd.all.na=predict(ermod.na,newdata=ern,type="response")
     }
     
     ppd.list <- list()
     i=0
     Eppd.all=predict(ermod,type="response")
     for (year in years)
     {
        i=i+1
        er=ern[ern$Start.year==year,]
        ppd=tapply(er$nhat, 
                   floor(er$time),sum)/tapply(er$effort,floor(er$time),sum)
        ppd.list[[i]] <- ppd
        if(plotit)
        {
           Eppd=Eppd.all[ern$Start.year==year]
           ymax=max(c(Eppd,ppd))*1.05
           plot(er$time,Eppd, 
                ylim=c(0,ymax), type="l",
                main=paste(year,"/",year+1,sep=""),
                xlab="Days since 1 Dec", ylab=ylabel,xlim=c(0,100))
           points(as.numeric(names(ppd)),ppd)
           
           # add ggplot and save them:
           library(ggplot2)
           pred.plot <- ggplot(data.frame(time = er$time,
                                          gam.fit = Eppd)) +
             geom_line(aes(x = time, y = gam.fit)) +
             geom_point(data = data.frame(date = as.numeric(names(ppd)),
                                          ppd = ppd),
                        aes(x = date,
                            y = ppd)) + 
             xlab("Days since 1 Dec") +
             ylab("Whales per day") +
             labs(title = paste0(year,"/",year+1))
           
           ggsave(plot = pred.plot,
                  filename = paste0("figures/Laake_", 
                                    paste0(year,"-",year+1), ".png"),
                  device = "png",
                  dpi = 600)
           
           if(anchor & show.anchor)
           {
           lines((floor(lower[i]):(final[i]-1))+.5,
             Eppd.all.na[ern$Start.year==year],
             xlim=c(lower[i],final[i]),lty=2)
        }
       }
     }
     newdata=NULL
     for(i in 1:length(years))
     {
        ndata=data.frame(time=(lower[i]:(final[i]-1))+.5)
        ndata$Start.year=factor(rep(years[i],nrow(ndata)),levels=years)
        if(length(grep("vis",as.character(formula)))!=0)
          ndata$vis=1
        if(length(grep("beaufort",as.character(formula)))!=0)
          ndata$beaufort=0
        newdata=rbind(newdata,ndata)
     }
     Eppd.all=predict(ermod,newdata=newdata,type="response")
     Total=tapply(Eppd.all,newdata$Start.year,sum)  
     var.Total=NULL 
     sumN=NULL 
     for(i in 1:length(years))
     {
       newd=newdata[newdata$Start.year==years[i],]
       Xp <- predict(ermod,newd,type="lpmatrix") 
       Xs=Xp*as.vector(exp(Xp%*%coef(ermod)))
       var.Total=c(var.Total,sum(Xs%*%ermod$Vp%*%t(Xs))+sum(Eppd.all*ermod$scale))
      }   
     return(list(models=ermod,
                 Total=Total,
                 var.Total=var.Total,
                 pred=Eppd.all,
                 ppd = ppd.list))
   }
}

lnl.missed.pods <-
function(par,data,primary.only,years=NULL,formula,gsS,debug)
{
# Computes negative log-likelihood for true pod size distribution and
# detection probability parameters for double observer survey data
#
# Arguments:
#
#  par     - parameter values
#  data    - match data ordered in pairs with station=P then S where
#            P=Primary & S=Secondary
#  primary.only - data from observations of primary observer only on watch
#  years   - vector of years included in the analysis
#  formula - the formula for the logistic detection model
#  gsS     - pod size calibration matrix
#  debug   - if TRUE will show iteration values
#
# Value: negative log-likelihood value at specified parameters
################################################################################
#
# Extract parameter values for detection model and set up probability matrix
# prob which has a row for each data pair and a column for each possible
# true pod size.
nmax=ncol(gsS)
pbyyear=FALSE
if(is.null(years)) 
  npar=2
else
{
  pbyyear=TRUE
  npar=2*length(levels(factor(years)))
}
beta=par[(npar+1):length(par)]
n=nrow(data)
# Loop over each possible true pod size and compute the probability of observing
# the value of the capture history if the true size was i for each given
# station's recorded pod size if seen by the observer at the station.
# See pdf writeup.
prob=sapply(1:nmax,function(x) 
{
  data$True=x
  dmat=model.matrix(formula,data)
  p=plogis(dmat%*%beta)
  gs=(gsS[x,][data$podsize])^data$seen
  pp=p^data$seen*(1-p)^(1-data$seen)
  return( gs[seq(1,n,2)]*gs[seq(2,n,2)]*pp[seq(1,n,2)]*pp[seq(2,n,2)]/
               (1-(1-p[seq(1,n,2)])*(1-p[seq(2,n,2)])))
})
# Next compute the overall probability by summing across all true pod sizes
# weighted by the estimated probability that a pod is of that true size (fS).
# The negative sum of the log(total prob) is the negative log-likelihood.
# If pbyyear is TRUE then this is done by survey year because the fS values
# differ by year.
if(!pbyyear)
{
  fS=gammad(par[1:2],nmax)
  neglnl=-sum(log(prob%*%fS))
  if(!is.null(primary.only))neglnl=neglnl+ lnl.pods(par[1:2], formula=formula, dpar=beta, 
                              data=primary.only, gsS=gsS, debug=debug)
}
else
{
  neglnl=sum((sapply(1:length(years), function(i,years){
     year=years[i]
     fS=gammad(par[(2*(i-1)+1):(2*i)],nmax)
     neglnl=-sum(log(prob[as.character(data$Start.year[seq(1,n,2)])==year,]%*%fS))
     if(!is.null(primary.only))
        neglnl=neglnl+ lnl.pods(par[(2*(i-1)+1):(2*i)], formula=formula, dpar=beta, 
            data=primary.only[as.character(primary.only$Start.year)==year,], gsS=gsS, debug=debug)
     return(neglnl)
     }, years=years)))
}
# if debug, output iteration results
if(debug)
{
   cat("\npar=",paste("c(",paste(par,collapse=","),")"))
   cat("\nneglnl=",neglnl)
}
return(neglnl)
}

fit.missed.pods <-
function(data, primary, pbyyear=FALSE, formula=~1, par=NULL, gsS,
                         maxit=1000,refit=TRUE,debug=FALSE,hessian=FALSE,method="BFGS")
{
# Fits true pod size distribution and the detection probability parameters
# for a given detection model specified by formula.
#
# Arguments:
#
#  data    - match data ordered in pairs with station=P then S where
#            P=Primary & S=Secondary
#  primary - observations made during primary period
#  pbyyear - if TRUE, fits year-specific gamma parameters for pod size
#  formula - the formula for the logistic detection model
#  par     - initial parameter values
#  gsS     - pod size calibration matrix
#  maxit   - maximum iterations
#  refit   - if TRUE will continue to call optim to refit until convergence is achieved
#             regardless of value of maxit
#  debug   - if TRUE will show iteration values
#  hessian - if TRUE will return hessian
#
#  Value:
#
#   list with elements : par  - parameter estimates
#                        AIC  - AIC value for model
#                        model- output from final optim run
################################################################################

#print("fit.missed.pods starting")
#Sys.sleep(3)

nmax=ncol(gsS)
primary.only=primary[primary$only,]
if(nrow(primary.only)==0)primary.only=NULL
# add a dummy True pod size
data$True=1

# set number of pod size parameters depending on pbyyear and number of survey years
if(pbyyear){
   if(is.null(primary.only)){
     years = levels(data$Start.year)
   } else {
     years = sort(unique(c(levels(primary.only$Start.year),levels(data$Start.year))))
     # np = 2 *length(years) - this was here but it's wrong. 
   }
   np = 2 *length(years)
   #print(paste("years = ", years))
} else {
   years=NULL
   np=2
}

#print(paste("np = ", np))

# set number of parameters in detection probability model and create formula
# using podsize in place of True for io.glm
nc=ncol(model.matrix(formula,data))
if(length(grep("True",as.character(formula)))!=0){
   dformula=as.formula(paste("seen",paste(sub("True","podsize",as.character(formula)),collapse=""),sep=""))
} else {
   dformula=as.formula(paste("seen",paste(as.character(formula),collapse=""),sep=""))
}
# Depending on value of par argument create initial values for parameters
# Unless all are specified, a glm is used to specify the starting values for
# the detection parameters.  That is done using observed pod size (podsize) values
# rather than unknown true pod sizes.
if(is.null(par) | length(par)==np){
  if(length(par)==np){
    par=c(par,as.vector(coef(io.glm(data,dformula))))
  } else {
    dpar=as.vector(coef(io.glm(data,dformula)))
    if(pbyyear){
      psmod=podsize.computations(primary,gsS=gsS)
      spar=as.vector(sapply(psmod,function(x)x$par))
      sspar=NULL
      for(i in 1:length(years))
        sspar = c(sspar,
					fit.pods(spar[(2*(i-1)+1):(2*i)], 
								formula=formula, 
								dpar=dpar, 
								data=primary[primary$Start.year==years[i],],  
								gsS=gsS, 
								debug=FALSE, 
								hessian=FALSE)$par)
    } else {
       sspar=c(.4,-.4)
    }
    par=c(sspar,dpar) 
  }
} else {
  if(length(par)>nc+np){
    stop("initial par vector too long")
  } else {
    if(length(par) <nc+np) stop("initial vector too short")
	}
}

if(debug) cat("\nInitial values = ",par,"\n")
# Fit model with optim which use lnl.missed.pods for the negative log-likelihood
mod=optim(par,lnl.missed.pods,years=years,formula=formula,data=data,primary.only=primary.only,
              gsS=gsS,control=list(maxit=maxit,ndeps=rep(1e-5,length(par))),debug=debug,hessian=hessian,method=method)
# If the model didn't converge and refit is TRUE, continue to fit the model
# using final values from last fit until it converges.
if(refit)
   while(mod$convergence!=0)
      mod=optim(mod$par,lnl.missed.pods,years=years,formula=formula,data=data,primary.only=primary.only,
                            gsS=gsS,control=list(maxit=maxit,ndeps=rep(1e-5,length(par))),debug=debug,hessian=hessian,method=method)
# Add names to each parameter in the model and return the result list
if(pbyyear)
{
  cnames=paste(rep(years,each=2),rep(c("Gamma shape","Gamma rate"),length(years)),sep=":")
  cnames=c(cnames,colnames(model.matrix(dformula,data)))
} else
  cnames=c("Gamma shape","Gamma rate",colnames(model.matrix(dformula,data)))
parvec=mod$par
names(parvec)=cnames
AIC=mod$value*2+2*length(cnames)
return(list(par=parvec,AIC=AIC,model=mod))
}


fit.pods=function(par, formula, dpar, data,  gsS, debug, hessian)
{
# Optimize log-likelihood of pod size parameters with specfied detection model
  optim(par=par, lnl.pods,formula=formula,data=data,dpar=dpar,gsS=gsS,debug=debug,hessian=hessian)
}

lnl.pods <-
function(par, formula, dpar, data,  gsS, debug)
{
# Computes negative log-likelihood for true pod size distribution
# for a given detection probability and parameters for a single observer
#
# Arguments:
#
#  par     - parameter values for pod size distribution
#  formula - the formula for the logistic detection model
#  dpar    - parameter values for detection model
#  data    - single observer data
#  gsS     - pod size calibration matrix
#  debug   - if TRUE will show iteration values
#
# Value: negative log-likelihood value at specified parameters
################################################################################
#
# Extract parameter values for detection model and set up probability matrix
# prob which has a row for each data pair and a column for each possible
# true pod size.
nmax=ncol(gsS)
fS=gammad(par,nmax)
probs=function(j)
{
  data$True=j
  dmat=model.matrix(formula,data)
  plogis(dmat%*%dpar)
}
ps=sapply(1:nmax, probs)
ps=ps/as.vector(ps%*%fS)
prob=rowSums(ps*t(fS*gsS[,data$podsize]))
neglnl=-sum(log(prob))
if(debug)
{
   cat("\npar=",par)
   cat("\nneglnl=",neglnl)
}
return(neglnl)
}

fit.poisson.re=function(x,re.factor,formula,sigma.formula=~1,initial,
                         maxit=1000,method="BFGS",hessian=FALSE,z=5,nmax=20)
{
# Fits a Poisson mixed-effects model to the pod size calibration data
# with random effects defined by re.factor.
#
# Arguments:
#
#  x             - dataframe
#  re.factor     - list of one or more vectors of factor variables to split x into
#                  a list of dataframes. Length of each factor variable vector must
#                  match number of rows in x.
#  formula       - formula for lambda of Poisson
#  sigma.formula - formula for sigma of normal error
#  initial       - required vector of initial parameter values
#  maxit         - maximum number of iterations (passed to optim)
#  method        - method used by optim for optimization
#  hessian       - passed to optim, if TRUE returns hessian for v-c matrix
#  z             -  defines integration bounds for normal (default is 5 for integration range of -5*sigma,5*sigma)
#  nmax          -  maximum pod size
#
#  Value: list returned with optim results
#
   mm=lapply(split(x,re.factor,drop=TRUE),function(x) model.matrix(formula,x))
   ss=lapply(split(x,re.factor,drop=TRUE),function(x) model.matrix(sigma.formula,x))
   observed.size=split(x$Estimate,re.factor,drop=TRUE)
   return(optim(initial,repois.lnl,method=method,control=list(maxit=maxit),hessian=hessian,
                 mm=mm,ss=ss,observed.size=observed.size,z=z,nmax=nmax))
}

repois.lnl=function(par,mm,ss,observed.size,z=5,nmax=20)
{
#
# Compute negative log-likelihood for a random effects Poisson model for
# the pod size calibration data. A normal 0,sigma error distribution is used
# for the random effect.  It is used for the log(lambda) where lambda is the
# Poisson intensity.
#
#  Arguments:
#   par           -  parameter values
#   mm            -  list of model matrices for lambda; list elements are split by random effect factor levels
#   ss            -  list of model matrices for sigma; list elements are split by random effect factor levels
#   observed.size -  list of vectors of observed sizes; list elements are split by random effect factor levels
#   z             -  defines integration bounds for normal (default is 5 for integration range of -5*sigma,5*sigma)
#   nmax          -  maximum pod size
#
#  Note that the mm,ss, observed.size lists need to be defined to split the data up
#  by the intersection of the factors used to define mm and ss such that for each grouping
#  there is only a single sigma.
#
#  Value:  negative log-likelihood
#
  pois.reint=function(x,obs,xbeta,sigma,nmax=20)
  {
#   Internal function to be used by integrate function for integration over normal random effect
#
#   Arguments:
#
#   x    - vector of values passed by integrate
#   obs  - observed size values
#   xbeta- x%*%beta where beta is parameter vector for lambda
#   sigma- value of sigma for normal
#   nmax - maximum pod size
#
    lambda=exp(outer(xbeta,x,"+"))
    return(dnorm(x,sd=sigma)*apply(lambda,2,function(lambda)prod(pois.d(lambda,x=obs,nmax=nmax))))
  }
# The parameter vector contains sigma parameters and then lambda parameters
# The number of columns of the design matrix ss determines number of sigma parameters
  npar=ncol(ss[[1]])
  neglnl=0
# Loop over each random effect
  for(i in 1:length(mm))
  {
#    Compute xbeta (x%*%beta) where x is design matrix and beta are parameters for lambda
     xbeta=as.vector(mm[[i]]%*%par[(npar+1):length(par)])
#    Compute sigma
     sigma=unique(exp(as.vector(ss[[i]]%*%par[1:npar])))
     if(length(sigma)>1)stop("\nInvalid sigma values\n")
     obs=observed.size[[i]]
#    Compute negative log likelihood by integrating this portion (subset of data defined by
#    random effects) of the likelihood over the normal distribution for the random effect
#    and accumulate the value.
     neglnl=neglnl-log(integrate(pois.reint,-z*sigma,z*sigma,obs=obs,xbeta=xbeta,sigma=sigma,nmax=nmax,stop.on.error=FALSE)$value)
  }
  cat("\npar = ",par)
  cat("\n neglnl= ",neglnl)
  return(neglnl)
}

fit.poisson=function(x,formula,initial,maxit=1000,method="BFGS",hessian=FALSE,nmax=20)
{
# Fits a Poisson fixed-effects model to the pod size calibration data.
#
# Arguments:
#
#  x             - dataframe
#  formula       - formula for lambda of Poisson
#  initial       - required vector of initial parameter values
#  maxit         - maximum number of iterations (passed to optim)
#  method        - method used by optim for optimization
#  hessian       - passed to optim, if TRUE returns hessian for v-c matrix
#  nmax         -  maximum pod size
#
#  Value: list returned with optim results
#
   xmat=model.matrix(formula,x)
   return(optim(initial,pois.lnl,method=method,control=list(maxit=maxit),
                hessian=hessian,xmat=xmat,observed.size=x$Estimate,nmax=nmax))
}

pois.lnl=function(par,xmat,observed.size,nmax=20)
{
# Compute negative log-likelihood for a fixed effects Poisson model for
# the pod size calibration data.
#
#  Arguments:
#   par           -  parameter values
#   xmat          -  model matrix for lambda
#   observed.size -  vector of observed sizes
#   nmax          -  maximum pod size
#
#  Value:  negative log-likelihood
#
  lambda=exp(xmat%*%par)
  neglnl=-sum(log(pois.d(lambda,x=observed.size,nmax=nmax)))
  cat("\npar = ",par)
  cat("\n neglnl= ",neglnl)
return(neglnl)
}

pois.d=function(lambda,x,nmax)
{
# Computes the Poisson probabilities for a vector of lambdas matched with observed values (x).
# Because these are pod size estimates (1,2,...) the observed x is shifted (x-1) such that it
# matches a Poisson with values 0,1,2,...  Also, to avoid dealing with an infinite number of
# of possible values the maximum pod size is set (nmax) and the distribution is renormalized
# over the range 0,1,2,...,(nmax-1).
#
# Arguments:
#   lambda        - vector of lambdas
#   x             - vector of observed pod sizes
#   nmax          -  maximum pod size
#
   if(length(x)!=length(lambda))stop("\n****invalid lengths****/n")
   num=apply(matrix(ppois(c(x-1,x-2),lambda=rep(lambda,2),lower=FALSE),ncol=2),1,diff)
   denom=ppois(rep(nmax-1,length(x)),lambda)
   denom[num<1e-16]=1
   return(num/denom)
}

io.glm <-
function(datavec, fitformula, eps = 0.00001, iterlimit = 500)
{
# ---------------------------------------------------------------
#  This is the code that uses the iterative offset glm or gam
#  approach; iteration is done until parameters are within a
#  certain epsilon (eps) or iteration limit (iterlimit) exceeded.
#
#  Note: David used offset in formula and I've put it as an
#  argument to glm and gam functions.
#
# Input : datavec = dataframe
#         fitformula = formula
#         eps = convergence criterion - fixed
#         iterlimit = maximum number of iterations allowed - fixed
#
# Output: list with
#       Note: modified to return glm object only with
#               class("ioglm","glm","lm")
#               class("ioglm","gam")
#
# $glmobj:  glm model
# $offsetvalue: final offsetvalues from iterative fit
# $plotobj: gam plot object (if GAM & gamplot==TRUE, else NULL)
# ----------------------------------------------------------------
#
  done <- FALSE
  i <- 1
  plotobj <- NULL
  if(is.null(datavec$offsetvalue))datavec$offsetvalue=0
  while(i <= iterlimit & !done) {
#  fit the glm or gam
    offsetvalue=datavec$offsetvalue
    ioglm <- glm(formula = fitformula, family = binomial, data = datavec, offset=offsetvalue)
    coeff <- ioglm$coeff
    fittedp <- ioglm$fitted.values

    if(i == 1) {
      oldmodel <- ioglm
      oldcoeff <- coeff
      oldp <- fittedp
    }else{
#    calculate differences between previous and present set of model outputs
      reldiff <- max(abs(plogis(coeff) - plogis(oldcoeff))/plogis(oldcoeff))

      if(is.na(reldiff)) {
        print("Can't calculate regression coefficients - model has not converged")
        print(" - last fit used for estimation" )
        ioglm <- oldmodel
        done <- TRUE
      }

      if(reldiff < eps & !done) {
        done <- TRUE
      }else{
        oldmodel <- ioglm
        oldcoeff <- coeff
        oldp <- fittedp
      }
    }
    if(!done){
      oldoff <- datavec$offsetvalue
        off <-  - log(plogis(predict(ioglm) - datavec$offsetvalue))
      datavec$offsetvalue[datavec$station == "S"] <- off[
        datavec$station == "P"]
      datavec$offsetvalue[datavec$station == "P"] <- off[
        datavec$station == "S"]
    }
    i <- i + 1
  }

  if(!done){
    datavec$offsetvalue <- oldoff
    warning("Iteration limit exceeded - last fit used for estimation")
  }

  class(ioglm)=c("ioglm",class(ioglm))

  return(ioglm)
}

################################################################################
# Creates dbf files for ERAbund with Recent survey data
################################################################################
CreateERAbundFiles=function(directory="c:/gw",DBDirectory="",years=c(87,92,93,95,97,0,1,6))
{
library(RODBC)
library(foreign)
# Create function to write out the ERAbund dbf file
write.ERAbundFile=function(year,directory=directory)
{
   extract=ERSurveyData[substr(ERSurveyData$SEQUENCE,1,4)==year,2:24]
   extract$DATE=as.Date(as.character(extract$DATE))
   extract$ETIME=as.single(extract$ETIME)
   extract$NTIME=as.single(extract$NTIME)
   extract$STIME=as.single(extract$STIME)
   extract$SRET=as.single(extract$SRET)
   extract$NRET=as.single(extract$NRET)
   extract$C_CPAIR=as.character(extract$C_CPAIR)
   extract=extract[extract$EFLAG!="",]
   write.dbf(extract,paste(directory,"\\/ERSW",substr(year,3,4),substr(year+1,3,4),".dbf",sep=""))
}
# Get recent survey data but first check that it can be found
if(DBDirectory != "")
{
   fdir = file.path(DBDirectory,"GrayWhaleSurveyData.accdb")
} else
{
   fdir = "GrayWhaleSurveyData.accdb"
}
if(!file_test("-f", fdir)) stop(paste("Cannot find file: ",fdir))
con=odbcConnectAccess2007(fdir)
ERSurveyData=sqlFetch(con,"AllRecentData")
# Next check to make sure that all of the ERAbund directories are found and write out the files
for (year in years)
{
  if(directory != "")
  {
   fdir = file.path(directory,paste("ERSW",formatC(year,width=2,flag="0"),formatC(year+1,width=2,flag="0"),sep=""))
  } else
  {
   fdir = file.path(paste("ERSW",formatC(year,width=2,flag="0"),formatC(year+1,width=2,flag="0"),sep=""))
  }
  if(!file_test("-d", fdir)) stop(paste("Cannot find directory: ",fdir))
  fyear=year+1900
  if(fyear<1960)fyear=fyear+100
  write.ERAbundFile(fyear,fdir)
}
# Next do the same for the Secondary subdirectory
for (year in years)
{
  if(directory != "")
  {
   fdir = file.path(directory,paste("Secondary/ERSW",formatC(year,width=2,flag="0"),formatC(year+1,width=2,flag="0"),sep=""))
  } else
  {
   fdir = file.path(paste("Secondary/ERSW",formatC(year,width=2,flag="0"),formatC(year+1,width=2,flag="0"),sep=""))
  }
  if(!file_test("-d", fdir)) stop(paste("Cannot find directory: ",fdir))
  fyear=year+1900
  if(fyear<1960)fyear=fyear+100
  write.ERAbundFile(fyear,fdir)
}
close(con)
return(NULL)
}
# Process ER survey data and create dataframes:
# Primary,PrimaryOff,Secondary,Match,PrimaryEffort,and PrimarySightings
# This function is called by store.ERAbundData which stores the constructed
process.ERdata=function(fdir="c:/gw",maxbeauf=4,maxvis=4)
{
# Get survey data from package directory
  EarlyEffort=NULL
  EarlySightings=NULL
  ERSurveyData=NULL
  data(EarlyEffort,envir=environment()) #,package="ERAnalysis",envir=environment())
  data(EarlySightings,envir=environment()) #,package="ERAnalysis",envir=environment())
  data(ERSurveyData,envir=environment()) #,package="ERAnalysis",envir=environment())
  card2=EarlyEffort
  card3=EarlySightings
  card3$original.watch=card3$Watchperiod
  card3$watch=card3$Assignedwatchperiod
# Read in analysis output files from ERAbund using ReadERAbundFiles
  ERAbundFileList=ReadERAbundFiles(fdir)
  Primary=ERAbundFileList$Primary
  Secondary=ERAbundFileList$Secondary
  Match=ERAbundFileList$Match
  er.migdata=subset(ERAbundFileList$ERNorm,select=c("Start.year","period","begin","end","npods","nwhales","effort","vis","beaufort"))
  er.migdata.sec=subset(ERAbundFileList$ERNormSec,select=c("Start.year","period","begin","end","npods","nwhales","effort","vis","beaufort"))
# Convert some of the fields; distance converted to km
# station is a 5 character string with only last character needed
# watch field created based on crossing time and made into a factor
  Match$distance=Match$distance/1000
  Match$station=factor(substr(as.character(Match$station),5,5))
  Match$original.watch=Match$watch
  Match$watch=1
  Match$watch[Match$t241>=10.5&Match$t241<13.5]=2
  Match$watch[Match$t241>=13.5]=3
  Match$watch=factor(Match$watch)
  Primary$distance=Primary$distance/1000                  # from meters to km
  Primary$original.watch=Primary$watch
  names(Primary)[names(Primary)=="pods.per.hour"]="pphr"
  Primary$etime=NULL
  Primary$station=NULL
  Primary$seen=NULL
  Primary$date=NULL
  Primary$effort.period=paste(Primary$Start.year,"_",formatC(Primary$effort.period,digits=3,flag=0),sep="")
  Primary$watch=1
  Primary$watch[Primary$t241>=10.5]=2
  Primary$watch[Primary$t241>=13.5]=3
  Primary$watch=factor(Primary$watch)
  x=Match[Match$station=="P"&Match$seen==1,]
  x$date=as.Date(paste(x$Start.year,"-12-01",sep=""))+x$day-1
  x$skey=paste(x$date,x$t241)
  Primary$skey=paste(Primary$year,"-",formatC(Primary$month,flag=0,digits=1),
    "-",formatC(Primary$day,flag=0,digits=1)," ",Primary$t241,sep="")
  Primary$only=TRUE
  Primary$only[Primary$skey%in%x$skey]=FALSE
  Primary$skey=NULL
# Create common dataframe for migration curve data
# For early data, create the gwnorm data structure by selecting effort periods
# that meet vis and beaufort criterion
  merged.cards=merge(card2,card3,all.x=TRUE,by.x="key",by.y="key")
# set max beauf and vis values to select data; in data values for Use
# are 4 for both
  badsight=sapply(split(merged.cards,merged.cards$key), function(x) any(x$Visibility>maxvis | x$Beaufort>maxbeauf))
  badsight[is.na(badsight)]=FALSE
  bb=subset(card2,select=c("Block1beaufort","Block1visibility","Block2beaufort","Block2visibility","Block3beaufort","Block3visibility"))
  badbf=apply(bb[,c(1,3,5)],1,function(x)any(x>maxbeauf))
  badvis=apply(bb[,c(2,4,6)],1,function(x)any(x>maxvis))
  badbf[is.na(badbf) | is.nan(badbf)]=FALSE
  badvis[is.na(badvis)| is.nan(badvis)]=FALSE
  card2$Use="Yes"
  card2$Use[badvis | badbf | badsight]="No"
# create an average value for vis and beaufort using the sightings and effort values for the older data
  avgviseff=apply(subset(card2,select=c("Block1visibility","Block2visibility","Block3visibility")),1,mean,na.rm=TRUE)
  avgvis=with(merged.cards, tapply(Visibility,key,mean,na.rm=TRUE))
  avgvis[is.na(avgvis)]=avgviseff[is.na(avgvis)]
  avgbfeff=apply(subset(card2,select=c("Block1beaufort","Block2beaufort","Block3beaufort")),1,mean,na.rm=TRUE)
  avgbf=with(merged.cards, tapply(Beaufort,key,mean,na.rm=TRUE))
  avgbf[is.na(avgbf)]=avgbfeff[is.na(avgbf)]
  card2$vis=as.vector(avgvis)
  card2$beaufort=as.vector(avgbf)
# merge effort and sightings again with new Use values
  merged.cards=merge(card2,card3,all.x=TRUE,by.x="key",by.y="key")
# construct migration effort dataframe for early years
  Use=NULL
  gwnorm=subset(card2,subset=Use=="Yes",select=c("key","Start.year","Year","Month","Day","Begin.date.time","End.date.time","vis","beaufort","Observer"))
# construct a dataframe of pod and whale counts from observations
  npods=with(merged.cards[merged.cards$Use=="Yes" &merged.cards$Included,], tapply(Podsize,key,length))
  nwhales=with(merged.cards[merged.cards$Use=="Yes" &merged.cards$Included,], tapply(Podsize,key,sum))
  df=data.frame(key=names(npods),npods=as.vector(npods),nwhales=as.vector(nwhales))
# merge the counts with the effort data
  gwnorm=merge(gwnorm,df,all.x=TRUE,by="key")
  gwnorm$npods[is.na(gwnorm$npods)]=0
  gwnorm$nwhales[is.na(gwnorm$nwhales)]=0
# compute length of effort period in decimal days and compute begin and
# end values for period in decimal days from 1 Dec
  gwnorm$begin=gwnorm$Begin.date.time-strptime(paste(gwnorm$Start.year,"-12-01 00:00",sep=""),format="%Y-%m-%d %H:%M")
  gwnorm$end=gwnorm$End.date.time-strptime(paste(gwnorm$Start.year,"-12-01 00:00",sep=""),format="%Y-%m-%d %H:%M")
  gwnorm$effort=as.vector(gwnorm$end)-as.vector(gwnorm$begin)
  gwnorm=subset(gwnorm,select=c("Start.year","key","begin","end","npods","nwhales","effort","vis","beaufort","Observer"))
# Create single primary sightings file with data from all years
  PrimarySightings=merge(subset(card3,select=-Observer),gwnorm,by.x="key",by.y="key")
  PrimarySightings$only=TRUE
  PrimarySightings=PrimarySightings[PrimarySightings$Included,]
  PrimarySightings$distance=PrimarySightings$Distanceoffshore*1.852 # from nm to km
  PrimarySightings$distance[is.na(PrimarySightings$distance)]=mean(PrimarySightings$distance[!is.na(PrimarySightings$distance)])
  PrimarySightings$pphr=PrimarySightings$npods/(PrimarySightings$effort*24)
  names(PrimarySightings)[10]="podsize"
  names(PrimarySightings)[15]="vis"
  names(PrimarySightings)[31]="Start.year"
  T241c=as.character(PrimarySightings$Timeabeam)
  PrimarySightings$t241=as.numeric(substr(T241c,1,2))+as.numeric(substr(T241c,4,5))/60+as.numeric(substr(T241c,7,8))/3600
  PrimarySightings=subset(PrimarySightings,select=c("Day","Month","Year","watch","t241","distance",
             "podsize","Observer","vis","Beaufort","Winddirection","key","pphr","Start.year","original.watch","only"))
  names(PrimarySightings)=names(Primary)
  PrimarySightings$pphr=as.numeric(PrimarySightings$pphr)
  PrimarySightings=rbind(PrimarySightings,Primary)
  names(PrimarySightings)[names(PrimarySightings)=="effort.period"]="key"
  names(Primary)[names(Primary)=="effort.period"]="key"
  PrimarySightings$watch=factor(PrimarySightings$watch)
# merge newer effort with older effort in gwnorm
  names(er.migdata)[names(er.migdata)=="period"]="key"
  er.migdata$key=paste(er.migdata$Start.year,"_",formatC(er.migdata$key,digits=3,flag=0),sep="")
  gwnorm$key=as.character(gwnorm$key)
#  er.migdata=rbind(subset(gwnorm,select=-Observer),er.migdata)
  er.migdata$Observer=0
  er.migdata=rbind(gwnorm,er.migdata)
  er.migdata$time=as.vector((er.migdata$end+er.migdata$begin)/2)
  er.migdata$begin=as.numeric(er.migdata$begin)
  er.migdata$end=as.numeric(er.migdata$end)
  er.migdata$key=factor(er.migdata$key)
# define watch 1-3 to select data based on
# the presence of any beauf>maxbeauf or vis>maxvis as in earlier years
  er.migdata$watch=1
  er.migdata$watch[er.migdata$Start.year<1987][(er.migdata$time[er.migdata$Start.year<1987]-floor(er.migdata$time[er.migdata$Start.year<1987]))>=(12/24)]=2
  er.migdata$watch[er.migdata$Start.year>=1987][(er.migdata$time[er.migdata$Start.year>=1987]-floor(er.migdata$time[er.migdata$Start.year>=1987]))>=(10.5/24)]=2
  er.migdata$watch[er.migdata$Start.year>=1987][(er.migdata$time[er.migdata$Start.year>=1987]-floor(er.migdata$time[er.migdata$Start.year>=1987]))>=(13.5/24)]=3
  er.migdata$watch.key=paste(er.migdata$Start.year,"_",floor(er.migdata$begin),er.migdata$watch,sep="")
  notUse=sapply(split(er.migdata,er.migdata$watch.key),function(x) any(x$vis>maxvis) | any(x$beaufort>maxbeauf))
  notUse[is.na(notUse)]=TRUE
  er.migdata=merge(er.migdata,data.frame(watch.key=names(notUse),Use=!notUse))
#
# Restrict all PrimarySightings & PrimaryEffort to be within maxbeauf and maxvis
  PrimaryEffort=er.migdata[(er.migdata$vis<=maxvis| is.na(er.migdata$vis))&(is.na(er.migdata$beaufort) | er.migdata$beaufort<=maxbeauf),]
  PrimarySightings=PrimarySightings[(is.na(PrimarySightings$vis)|PrimarySightings$vis<=maxvis)&(is.na(PrimarySightings$beaufort)|PrimarySightings$beaufort<=maxbeauf),]
# Create observer experience dataframe and add the experience field in hours to Match and
# PrimarySightings
  ObserverExp=CreateObserverExperience()
  data(Observer,envir=environment()) #,package="ERAnalysis",envir=environment())
  PrimarySightings$Date=paste(PrimarySightings$year,"-",formatC(PrimarySightings$month,width=2,flag=0),"-",formatC(PrimarySightings$day,width=2,flag=0),sep="")
  PrimarySightings$observer=substr(PrimarySightings$observer,nchar(PrimarySightings$observer)-2,nchar(PrimarySightings$observer))
  Observer$key=as.character(Observer$Initials)
  Observer$key[is.na(Observer$Initials)]=Observer$Observer[is.na(Observer$Initials)]
  PrimarySightings=merge(PrimarySightings,ObserverExp,by="Date",all.x=TRUE)
  PrimarySightings$hours=PrimarySightings[cbind(1:nrow(PrimarySightings),match(PrimarySightings$observer,names(PrimarySightings)[18:76])+17)]
  PrimarySightings=PrimarySightings[,c(1:17,77)]
  PrimarySightings$seq=1:nrow(PrimarySightings)
  PrimarySightings=merge(PrimarySightings, subset(Observer, select = c("key", "Sex")), by.x ="observer",by.y = "key", all.x = TRUE)
  PrimarySightings=PrimarySightings[order(PrimarySightings$seq),]
  PrimarySightings$seq=NULL
  PrimarySightings$hours=as.numeric(PrimarySightings$hours)
  PrimaryEffort$Date=substr(as.character(strptime(paste(PrimaryEffort$Start.year,"-12-01",sep=""),format="%Y-%m-%d")+3600*24*(floor(PrimaryEffort$begin)-1)),1,10)
  PrimarySightings$observer[PrimarySightings$observer=="33"]="CDV"
  PrimarySightings$observer[PrimarySightings$observer=="35"]="DJR"
  PrimarySightings$Observer=factor(PrimarySightings$observer)
  PrimarySightings$observer=NULL
  Match$Date=as.character(strptime(paste(Match$Start.year,"-12-01",sep=""),format="%Y-%m-%d")+3600*24*(Match$day-1))
  Match=merge(Match,ObserverExp,by="Date",all.x=TRUE)
  Match$observer=as.character(Match$observer)
  Match$observer=substr(Match$observer,nchar(Match$observer)-2,nchar(Match$observer))
  Match$hours=Match[cbind(1:nrow(Match),match(Match$observer,names(Match)[18:76])+17)]
  Match=Match[,c(1:17,77)]
  Match$hours=as.numeric(Match$hours)
  Match$seq=1:nrow(Match)
  Match=merge(Match, subset(Observer, select = c("key", "Sex")), by.x ="observer",by.y = "key", all.x = TRUE)
  Match=Match[order(Match$seq),]
  Match$seq=NULL
  Match$observer[Match$observer=="33"]="CDV"
  Match$observer[Match$observer=="35"]="DJR"
  Match$Observer=factor(Match$observer)
  Match$observer=NULL
# Create Secondary effort file
  names(er.migdata.sec)[names(er.migdata.sec)=="period"]="key"
  er.migdata.sec$key=paste(er.migdata.sec$Start.year,"_",formatC(er.migdata.sec$key,digits=3,flag=0),sep="")
  er.migdata.sec$time=as.vector((er.migdata.sec$end+er.migdata.sec$begin)/2)
  er.migdata.sec$begin=as.numeric(er.migdata.sec$begin)
  er.migdata.sec$end=as.numeric(er.migdata.sec$end)
  er.migdata.sec$key=factor(er.migdata.sec$key)
# define watch 1-3 to select data based on
# the presence of any beauf>maxbeauf or vis>maxvis as in earlier years
  er.migdata.sec$watch=1
  er.migdata.sec$watch[er.migdata.sec$Start.year<1987][(er.migdata.sec$time[er.migdata.sec$Start.year<1987]-floor(er.migdata.sec$time[er.migdata.sec$Start.year<1987]))>=(12/24)]=2
  er.migdata.sec$watch[er.migdata.sec$Start.year>=1987][(er.migdata.sec$time[er.migdata.sec$Start.year>=1987]-floor(er.migdata.sec$time[er.migdata.sec$Start.year>=1987]))>=(10.5/24)]=2
  er.migdata.sec$watch[er.migdata.sec$Start.year>=1987][(er.migdata.sec$time[er.migdata.sec$Start.year>=1987]-floor(er.migdata.sec$time[er.migdata.sec$Start.year>=1987]))>=(13.5/24)]=3
  er.migdata.sec$watch.key=paste(er.migdata.sec$Start.year,"_",floor(er.migdata.sec$begin),er.migdata.sec$watch,sep="")
  notUse=sapply(split(er.migdata.sec,er.migdata.sec$watch.key),function(x) any(x$vis>maxvis) | any(x$beaufort>maxbeauf))
  er.migdata.sec=merge(er.migdata.sec,data.frame(watch.key=names(notUse),Use=!notUse))
#
# Restrict all SecondaryEffort to be within maxbeauf and maxvis
  SecondaryEffort=er.migdata.sec[(er.migdata.sec$vis<=maxvis| is.na(er.migdata.sec$vis))&(is.na(er.migdata.sec$beaufort) | er.migdata.sec$beaufort<=maxbeauf),]
  SecondaryEffort$Date=substr(as.character(strptime(paste(SecondaryEffort$Start.year,"-12-01",sep=""),format="%Y-%m-%d")+3600*24*(floor(SecondaryEffort$begin)-1)),1,10)
  convert.decimal.time=function(x)
  {
  xx=(x-floor(x))*60
  zz=trunc((xx-floor(xx))*60+.5)
  return(paste(formatC(floor(x),width=2,flag=0),formatC(floor(xx),width=2,flag=0),formatC(floor(zz),width=2,flag=0),sep=":"))
  }
  data(ERSurveyData,envir=environment()) #,package="ERAnalysis",envir=environment())
  nn=names(PrimaryEffort)
# 1987 surveys and beyond -- adding Observer to Primary & Secondary Effort data
# Primary Effort EFLAG=1,2 records
  EFLAG=NULL
  EXPERIMENT=NULL
  LOCATION=NULL
  Start.year=NULL
  xx=subset(ERSurveyData,subset=EFLAG<=2&(EXPERIMENT==1 | (EXPERIMENT==2 & ((LOCATION=="S" & !Start.year %in% 2000:2001)|(Start.year %in% 2000:2001 & LOCATION=="N")) )))
  xx=xx[xx$VISCODE<=maxvis&xx$WINDFORCE<=maxbeauf,]
  xx$Start.year=as.numeric(xx$Start.year)
  xx$key=paste(xx$DATE,"_",substr(convert.decimal.time(xx$ETIME),1,2),sep="")
  PrimaryEffort$seq=1:nrow(PrimaryEffort)
  PrimaryEffort$wkey=paste(PrimaryEffort$Date,substr(convert.decimal.time(24*(PrimaryEffort$begin-floor(PrimaryEffort$begin)+.5/3600)),1,2),sep="_")
  zz=merge(PrimaryEffort,subset(xx,select=c("key","OBSERVER","ETIME")),by.x="wkey",by.y="key",all.x=TRUE)
  zz$tt=24*(zz$begin-floor(zz$begin))
  zz=zz[abs(zz$tt-zz$ETIME)<=1.5/3600,]
  zz=zz[!is.na(zz$wkey),]
  zz=zz[order(zz$seq),]
  PrimaryEffort$Observer=as.character(PrimaryEffort$Observer)
  PrimaryEffort$Observer[PrimaryEffort$Observer==0]=as.character(zz$OBSERVER)
  PrimaryEffort$Observer=factor(PrimaryEffort$Observer)
  PrimaryEffort=subset(PrimaryEffort,select=nn)
# Secondary Effort EFLAG=1,2 records
  xx=subset(ERSurveyData,subset=EFLAG<=2& (EXPERIMENT==2 & ((LOCATION=="S" & Start.year %in% 2000:2001)|(!Start.year %in% 2000:2001 & LOCATION=="N")) ))
  xx=xx[xx$VISCODE<=maxvis&xx$WINDFORCE<=maxbeauf,]
  xx$Start.year=as.numeric(xx$Start.year)
  xx$key=paste(xx$DATE,"_",substr(convert.decimal.time(xx$ETIME),1,2),sep="")
  SecondaryEffort$seq=1:nrow(SecondaryEffort)
  SecondaryEffort$wkey=paste(SecondaryEffort$Date,substr(convert.decimal.time(24*(SecondaryEffort$begin-floor(SecondaryEffort$begin)+.5/3600)),1,2),sep="_")
  zz=merge(SecondaryEffort,subset(xx,select=c("key","OBSERVER","ETIME")),by.x="wkey",by.y="key",all.x=TRUE)
  zz$tt=24*(zz$begin-floor(zz$begin))
  zz=zz[abs(zz$tt-zz$ETIME)<=1.5/3600,]
  zz$Observer=zz$OBSERVER
  SecondaryEffort=zz
  SecondaryEffort=subset(SecondaryEffort,select=nn)
  Primary$Observer=factor(Primary$observer)
  Primary$observer=NULL
  return(list(Primary=Primary,PrimaryOff=ERAbundFileList$PrimaryOff,Secondary=Secondary,Match=Match,PrimarySightings=PrimarySightings,
               PrimaryEffort=PrimaryEffort,SecondaryEffort=SecondaryEffort))
}
################################################################################
# Reads in output data files created by ERAbund with recent survey data
# The function returns a list of 4 dataframes:
#  1) Primary   - on-effort sightings from primary platform
#  2) PrimaryOff  - off-effort sightings from primary platform
#  3) Secondary - on and off-effort sightings from secondary platform
#  4) Match     - matches for double-observer analysis
#  5) ERNorm    - data structure for GWNorm program
################################################################################
ReadERAbundFiles=function(directory="c:/gw",years=c(87,92,93,95,97,0,1,6))
{
Primary=NULL
PrimaryOff=NULL
Secondary=NULL
Match=NULL
ERNorm=NULL
ERNormSec=NULL
for (year in years)
{
  yearc=paste(formatC(year,width=2,flag="0"),formatC(year+1,width=2,flag="0"),sep="")
# Read in primary observation file
  if(directory != "")
  {
   fdir = file.path(directory,paste("ERSW",yearc,sep=""),paste("erprim",yearc,".dat",sep=""))
  } else
  {
   fdir = file.path(paste("ERSW",yearc,sep=""),paste("erprim",yearc,".dat",sep=""))
  }
  if(!file_test("-f", fdir)) stop(paste("Cannot find directory: ",fdir))
  Primaryx=read.fwf(fdir,widths=c(5,5,5,9,4,11,10,4,5,4,5,5,8,7,5,5),skip=1)
  names(Primaryx)=c("day","month","year","etime","watch","t241","distance","podsize",
       "observer","vis","beaufort","wind.direction","effort.period","pods.per.hour","station","seen")
  Primaryx$Start.year=year+1900
  Primaryx$Start.year[Primaryx$Start.year<1987]=Primaryx$Start.year[Primaryx$Start.year<1987]+100
  Primary=rbind(Primary,Primaryx)
# Read in primary observation 'not used' file
  if(directory != "")
  {
   fdir = file.path(directory,paste("ERSW",yearc,sep=""),paste("erprimUAE",yearc,".dat",sep=""))
  } else
  {
   fdir = file.path(paste("ERSW",yearc,sep=""),paste("erprimUAE",yearc,".dat",sep=""))
  }
  if(!file_test("-f", fdir)) stop(paste("Cannot find directory: ",fdir))
  Primaryx=read.fwf(fdir,widths=c(5,5,5,9,4,11,10,4,5,4,5,5,8,7,5,5),skip=1)
  names(Primaryx)=c("day","month","year","etime","watch","t241","distance","podsize",
       "observer","vis","beaufort","wind.direction","effort.period","pods.per.hour","station","seen")
  Primaryx$Start.year=year+1900
  Primaryx$Start.year[Primaryx$Start.year<1987]=Primaryx$Start.year[Primaryx$Start.year<1987]+100
  PrimaryOff=rbind(PrimaryOff,Primaryx)
# Read in secondary observation file
  if(directory != "")
  {
   fdir = file.path(directory,paste("ERSW",yearc,sep=""),paste("ersec",yearc,".dat",sep=""))
  } else
  {
   fdir = file.path(paste("ERSW",yearc,sep=""),paste("ersec",yearc,".dat",sep=""))
  }
  if(!file_test("-f", fdir)) stop(paste("Cannot find directory: ",fdir))
  Secondaryx=read.fwf(fdir,widths=c(5,5,5,9,11,12,4,4,3,5,10))
  names(Secondaryx)=c("day","month","year","etime","t241","distance","podsize",
       "vis","beaufort","wind.direction","off")
  Secondaryx$Start.year=year+1900
  Secondaryx$Start.year[Secondaryx$Start.year<1987]=Secondaryx$Start.year[Secondaryx$Start.year<1987]+100
  Secondary=rbind(Secondary,Secondaryx)
# Read in match observation file
  if(directory != "")
  {
   fdir = file.path(directory,paste("ERSW",yearc,sep=""),paste("ermatchsp",yearc,".dat",sep=""))
  } else
  {
   fdir = file.path(paste("ERSW",yearc,sep=""),paste("ermatchsp",yearc,".dat",sep=""))
  }
  if(!file_test("-f", fdir)) stop(paste("Cannot find directory: ",fdir))
  Matchx=read.fwf(fdir,widths=c(5,5,5,4,9,10,6,4,4,5,5,8,9,6),skip=1)
  names(Matchx)=c("seen","station","day","watch","t241","distance","podsize","observer",
       "vis","beaufort","wind.direction","pphr","mscore","mcode")
  Matchx$Start.year=year+1900
  Matchx$Start.year[Matchx$Start.year<1987]=Matchx$Start.year[Matchx$Start.year<1987]+100
  Match=rbind(Match,Matchx)
# Read in ernorm file
  if(directory != "")
  {
   fdir = file.path(directory,paste("ERSW",yearc,sep=""),paste("ernorm",yearc,".dat",sep=""))
  } else
  {
   fdir = file.path(paste("ERSW",yearc,sep=""),paste("ernorm",yearc,".dat",sep=""))
  }
  if(!file_test("-f", fdir)) stop(paste("Cannot find directory: ",fdir))
  normx=read.fwf(fdir,widths=c(7,12,11,4,4,rep(3,13),4))
  names(normx)=c("period","begin","end","npods","nwhales",paste("pod",c(1:10,"11+"),sep=""),"vis","beaufort","wind.direction")
  normx$Start.year=year+1900
  normx$Start.year[normx$Start.year<1987]=normx$Start.year[normx$Start.year<1987]+100
  ERNorm=rbind(ERNorm,normx)
  # Read in secondary ernorm file
  if(directory != "")
  {
   fdir = file.path(directory,paste("Secondary/ERSW",yearc,sep=""),paste("ernorm",yearc,".dat",sep=""))
  } else
  {
   fdir = file.path(paste("Secondary/ERSW",yearc,sep=""),paste("ernorm",yearc,".dat",sep=""))
  }
  if(!file_test("-f", fdir)) stop(paste("Cannot find directory: ",fdir))
  normx=read.fwf(fdir,widths=c(7,12,11,4,4,rep(3,13),4))
  names(normx)=c("period","begin","end","npods","nwhales",paste("pod",c(1:10,"11+"),sep=""),"vis","beaufort","wind.direction")
  normx$Start.year=year+1900
  normx$Start.year[normx$Start.year<1987]=normx$Start.year[normx$Start.year<1987]+100
  ERNormSec=rbind(ERNormSec,normx)
}
Primary$wind.direction=factor(Primary$wind.direction)
Primary$year=Primary$year+1900
Primary$year[Primary$year<1987]=Primary$year[Primary$year<1987]+100
Primary$date=as.Date(paste(Primary$year,Primary$month,Primary$day,sep="-"))
Secondary$wind.direction=factor(Secondary$wind.direction)
Secondary$year=Secondary$year+1900
Secondary$year[Secondary$year<1987]=Secondary$year[Secondary$year<1987]+100
Secondary$date=as.Date(paste(Secondary$year,Secondary$month,Secondary$day,sep="-"))
Match$observer=factor(Match$observer)
Match$wind.direction=factor(Match$wind.direction)
ERNorm$effort=as.vector(ERNorm$end-ERNorm$begin)
ERNormSec$effort=as.vector(ERNormSec$end-ERNormSec$begin)
return(list(Primary=Primary,PrimaryOff=PrimaryOff,Secondary=Secondary,Match=Match,ERNorm=ERNorm,ERNormSec=ERNormSec))
}
store.ERdata=function(DBDirectory="",package.dir="C:/Users/Jeff Laake/Desktop/MyDocuments/R Development/ERAnalysis/ERAnalysis/data",nmax=20)
{
# Retrieves data from Gray Whale Access Database and creates dataframes
# EarlyEffort, EarlySightings for surveys from 1967-1985 and ERSurveyData for
# surveys from 1987-2006.  Some computed fields are added to the dataframes and then
# they are saved to the ERAnalysis package data directory.
  library(RODBC)
  # Set DBDirectory to point to Access Database. "" means in local directory
  DBDirectory=""
  if(DBDirectory != "")
  {
     fdir = file.path(DBDirectory,"GrayWhaleSurveyData.accdb")
  } else
  {
     fdir = "GrayWhaleSurveyData.accdb"
  }
  if(!file_test("-f", fdir)) stop(paste("Cannot find file: ",fdir))
  con=odbcConnectAccess2007(fdir)
# get effort and sightings tables and save them in package directory
  EarlyEffort=sqlFetch(con,"EffortPre1987")
  EarlySightings=sqlFetch(con,"ERSightingsPre1987")
  ERSurveyData=sqlFetch(con,"AllRecentData")
  EarlyEffort$Begin.date.time=as.POSIXct(with(EarlyEffort,strptime(paste(Year,"-",Month,"-",Day," ",formatC(Begintime,digits=3,flag="0"),sep=""),
                                  format="%Y-%m-%d %H%M")))
  EarlyEffort$End.date.time=as.POSIXct(with(EarlyEffort,strptime(paste(Year,"-",Month,"-",Day," ",formatC(Endtime,digits=3,flag="0"),sep=""),
                                  format="%Y-%m-%d %H%M")))
  EarlyEffort$Start.year=EarlyEffort$Year
  EarlyEffort$Start.year[EarlyEffort$Month<4]=EarlyEffort$Start.year[EarlyEffort$Month<4]-1
  EarlySightings$Start.year=EarlySightings$Year
  EarlySightings$Start.year[EarlySightings$Month<4]=EarlySightings$Start.year[EarlySightings$Month<4]-1
  ERSurveyData$Start.year=substr(ERSurveyData$SEQUENCE,1,4)
  ERSurveyData$year=substr(ERSurveyData$DATE,1,4)
  ERSurveyData$month=substr(ERSurveyData$DATE,6,7)
  ERSurveyData$day=substr(ERSurveyData$DATE,9,10)
  save(EarlyEffort,file=file.path(package.dir,"EarlyEffort.rda"))
  save(EarlySightings,file=file.path(package.dir,"EarlySightings.rda"))
  save(ERSurveyData,file=file.path(package.dir,"ERSurveyData.rda"))
  Observer=sqlFetch(con,"ObserverCodes")
  AddObs=subset(Observer,subset=Observer%in%c(33,35))
  AddObs$Initials=NA
  Observer=rbind(Observer,AddObs)
  save(Observer,file=file.path(package.dir,"Observer.rda"))
  close(con)
  # Read in calibration data from text file and create and store dataframes
  PodsizeCalibrationData=read.delim(file.path(package.dir,"PScaliball.txt"),header=TRUE)
  PodsizeCalibrationData$key=paste(PodsizeCalibrationData$Type,PodsizeCalibrationData$Year,PodsizeCalibrationData$ID,sep="_")
  PodsizeCalibrationTable=as.data.frame(table(PodsizeCalibrationData$key,factor(PodsizeCalibrationData$Estimate,levels=1:nmax))[,1:nmax])
  PodsizeCalibrationTable$key=row.names(PodsizeCalibrationTable)
  row.names(PodsizeCalibrationTable)=NULL
  PodsizeCalibrationTable=merge(unique(data.frame(key=PodsizeCalibrationData$key,
        True=PodsizeCalibrationData$True)),PodsizeCalibrationTable,by="key")
  save(PodsizeCalibrationTable,file=file.path(package.dir,"PodsizeCalibrationTable.rda"))
  save(PodsizeCalibrationData,file=file.path(package.dir,"PodsizeCalibrationData.rda"))
}
store.ERAbundData=function(fdir="c:/gw",package.dir="C:/Users/Jeff Laake/Desktop/MyDocuments/R Development/ERAnalysis/ERAnalysis/data",
                                 maxbeauf=4,maxvis=4)
{
   df.list=process.ERdata(fdir=fdir,maxbeauf=maxbeauf,maxvis=maxvis)
   SecondarySightings=df.list$Secondary
   Primary=df.list$Primary
   PrimarySightings=df.list$PrimarySightings
   PrimaryEffort=df.list$PrimaryEffort
   SecondaryEffort=df.list$SecondaryEffort
   Match=df.list$Match
   PrimaryOff=df.list$PrimaryOff
   save(Primary,file=file.path(package.dir,"Primary.rda"))
   save(PrimarySightings,file=file.path(package.dir,"PrimarySightings.rda"))
   save(PrimaryEffort,file=file.path(package.dir,"PrimaryEffort.rda"))
   save(SecondaryEffort,file=file.path(package.dir,"SecondaryEffort.rda"))
   save(Match,file=file.path(package.dir,"Match.rda"))
   save(SecondarySightings,file=file.path(package.dir,"SecondarySightings.rda"))
   save(PrimaryOff,file=file.path(package.dir,"PrimaryOff.rda"))
   df.list=process.ERdata(fdir=fdir,maxbeauf=10,maxvis=10)
   AllPrimaryEffort=df.list$PrimaryEffort
   AllSecondaryEffort=df.list$SecondaryEffort
   save(AllPrimaryEffort,file=file.path(package.dir,"AllPrimaryEffort.rda"))
   save(AllSecondaryEffort,file=file.path(package.dir,"AllSecondaryEffort.rda"))
}

##########################################################################
#  Function linklist links primary and secondary observer records in Match.list
#  creates revised dataframe with t241 and distance averaged and group size summed
#  three possible types of spatial criterion are included:
# "box" criterion is maximum of distances
# "diamond" criterion is sum of distances
# "ellipse" criterion is the euclidean distance default is ellipse.
#  Three specific criteria based on prior analysis are included:
#  "ERAbund" Uses the diamond spatial criterion with user specified parameters,
#         default is the criteria used by Breiwick et all 2006
# 0.04   !Time weight factor for linking
# 1.75   !Distance weight factor for linking
# 0.00   !Pod size weight factor for linking
#-0.05   !Maximum sum score for linking (merging pods)
## Rugh et al. 1993
# Rugh et al. 1990
#  RCH 6/12/2009
##################################################################
linklist <- function(mspdf,lcut=-0.05,twt=.18,dwt=3.95,crittype="ellipse")
{
# Arguments:
#
# mspdf  - dataframe of sightings
# t2dist - multiplier to convert time to distancein km/hr
# crittype - criterion for linking
# lcut   - linking parameter cutoff
# twt  -time weight for crossing time difference in minutes
# dwt   distance weight for ratio of difference to longer distance of shore at 241
#
mdf <- as.data.frame(mspdf)
Pdf <- mdf[(mdf$station=="P" & mdf$t241>1 ),]
Sdf <- mdf[(mdf$station=="S" & mdf$t241>1 ),]
lenp <- nrow(Pdf)
lens <- nrow(Sdf)
if (lenp>=2 & lcut>= 0) {
    Pdf <- linkdf(Pdf,lcut,twt,dwt,crittype)
    lenp <- nrow(Pdf)
    }
if (lens>=2 & lcut>= 0) {
    Sdf <- linkdf(Sdf,lcut,twt,dwt,crittype)
    lens <- nrow(Sdf)
    }
linklist<-rbind(Pdf,Sdf)
return(linklist)
}

##########################################################################
#  Function linkdf identifies potential links within observer record and
#  creates revised dataframe with t241 and dist averaged and group size summed
#  three possible types of spatial criterion are included:
# "box" criterion is maximum of distances
# "diamond" criterion is sum of distances
# "ellipse" criterion is the euclidean distance default is ellipse.
#  Three specific criteria based on prior analysis are included:
#  "ERAbund" Uses the diamond spatial criterion with user specified parameters,
#         default is the criteria used by Breiwick et all 2006
# 0.04   !Time weight factor for linking
# 1.75   !Distance weight factor for linking
# 0.00   !Pod size weight factor for linking
#-0.05   !Maximum sum score for linking (merging pods)
# Rugh et al. 1993
# Rugh et al. 1990
#  RCH 6/12/2009
####################################################################
linkdf <- function(mspdf,lcut=-0.05,twt=.18,dwt=3.95,crittype="ellipse")
{
# Arguments:
#
# mspdf  - dataframe of sightings
# t2dist - multiplier to convert time to distancein km/hr
# crittype - criterion for linking
# lcut   - linking parameter cutoff
# twt  -time weight for crossing time difference in minutes
# dwt   distance weight for ratio of difference to longer distance of shore at 241
#
mdf <- as.data.frame(mspdf)
lenm <- nrow(mdf)
if (lenm<=1) return(mdf)
tmat <- twt*60*abs(outer(mdf$t241, mdf$t241, "-"))
 fd <- function(x, y) (x-y)/max(x,y)
 dmat <- dwt*abs(outer(mdf$distance, mdf$distance, fd))
if (crittype == "box")   cmat<- pmax(tmat,dmat)
else
  if (crittype == "diamond")   cmat<- tmat + dmat
else
  if (crittype == "ERAbund")     cmat<- tmat + dmat
else   cmat<- sqrt(tmat^2 + dmat^2)
#generate link matrix
 c1mat <- 1*( cmat < lcut & cmat > 0)
#generate conversion matrix
 comat<-c1mat + diag(lenm)
#count entries byrow
crvec<-  rowSums(comat)
ccvec<-  colSums(comat)
#weight columns
conmat<-  comat * sqrt(outer(1/crvec,1/ccvec))
#sum weights by row
 smvec <-  rowSums(conmat)
#check sum if all rows are prperly weighted then finished
 while(length(smvec[smvec != 1])>0)
 {
#find improperly weighted rows
   snvec<- (smvec!=1)
#mask matrix to get subset of rows to be fixed
   fmat<-outer(snvec,snvec)  *cmat * c1mat
 # find and remove maximum link distance
   c1mat<-c1mat * (fmat!= max(fmat))
# create new conversion matrix
   comat<-c1mat + diag(lenm)
crvec<-  rowSums(comat)
ccvec<-  colSums(comat)
conmat<-  comat * sqrt(outer(1/crvec,1/ccvec))
  #check sum
     smvec <-  rowSums(conmat)
 }
#find rows with elements in lower triangle indicating duplicates
convec<- rowSums(lower.tri(conmat)*conmat)
# generate linked db
linkdf<-mdf
#t241 and dist averaged
linkdf$t241<- as.vector(conmat%*%(mdf$t241 ))
linkdf$distance<- as.vector(conmat%*%(mdf$distance ))
#ps summed
sumps<- (1*(conmat>0))
linkdf$podsize<- as.vector(sumps%*%(mdf$podsize))
if(!is.null(mdf$corrected.podsize))linkdf$corrected.podsize<- as.vector(sumps%*%(mdf$corrected.podsize))
#remove duplicates
linkdf<- linkdf[ convec==0,]
return(linkdf)
}

##########################################################################
#  Function matchdf identifies potential links within a linked observer record and
#  creates a revised dataframe with matches ordered by primary then secondary
#  column 1 is seen column 2 is 1 for primary and 2 for secondary
#  three possible types of spatial criterion are included:
# "box" criterion is maximum of distances
# "diamond" criterion is sum of distances
# "ellipse" criterion is the euclidean distance default is ellipse.
#  Three specific criteria based on prior analysis are included:
#  "ERAbund" Uses the diamond spatial criterion with user specified parameters,
#         default is the criteria used by Breiwick et all 200??
#0.04         !Time weight factor for matching     time difference in minutes
#1.75         !Distance weight factor for matching
#0.05         !Pod size weight factor for matching
# 1.0         !Maximum sum score for matching
# Rugh et al. 1993
# Rugh et al. 1990
#  RCH 7/21/2009
####################################################################

 matchdf<- function(mspdf,mcut=1,twt=.18,dwt=3.95,pwt=.05,crittype="ellipse")
{
# Arguments:
#
# mspdf  - dataframe of sightings
# t2dist - multiplier to convert time to distancein km/hr
# crittype - criterion for matching
# mcut   - matching parameter cutoff
# twt  -time weight for crossing time difference in minutes
# dwt   distance weight for ratio of difference to longer distance of shore at 241
# pwt   weight for difference in podsize
#
mdf <- as.data.frame(mspdf)
Pdf <- mdf[(mdf$station=="P"),]
Sdf <- mdf[(mdf$station=="S"),]
lenp <- nrow(Pdf)
lens <- nrow(Sdf)
if((lenp*lens)>0)  {
tmat <- twt*60*abs(outer(Pdf$t241, Sdf$t241, "-"))
 fd <- function(x, y) (x-y)/max(x,y)
 dmat <- dwt*abs(outer(Pdf$distance, Sdf$distance, fd))
 pmat <- pwt*abs(outer(Pdf$podsize, Sdf$podsize, "-"))
if (crittype == "box")     cmat<- pmax(tmat,dmat) + pmat
else
  if (crittype == "diamond")     cmat<- tmat + dmat + pmat
else
  if (crittype == "ERAbund")     cmat<- tmat + dmat + pmat
  else     cmat<- sqrt(tmat^2 + dmat^2) + pmat
# sort by best and eliminate sequentialy from best to worst
omat<- cbind(c(cmat),rep(1:lenp,lens),rep(1:lens,each=lenp))
      ii<-order(omat[,1],omat[,2])
     omat1<-rbind(omat[ii,],9999999,9999999)
     matched<-rbind(c(0,0,0),c(0,0,0))
     while(omat1[1,1] <= mcut){
     om2<-  omat1[1,2]
     om3<-  omat1[1,3]
     matched<-rbind(matched,c(1,omat1[1,2:3] ))
     omat1<-omat1[omat1[,2]!= om2,]
     omat1<-omat1[omat1[,3]!= om3,]
       } #end while
# nomatch vectors
crvec<- (1:lenp %in% matched[,2])                       
ccvec<- (1:lens %in% matched[,3])                       
xr<-cbind(crvec,2,1:lenp,0)
xc<-cbind(ccvec,3,0,1:lens)
#build matchdf
matches<-rbind(matched[ matched[,1]==1,],xr[xr[,1]==0,2:4],xc[xc[,1]==0,2:4])
} # end of if((lenp*lens)>0)
# generate matches set if no sightings by one observer
else {
if (lenp >0) matches<-cbind(2,1:lenp,0)
else matches<-cbind(3,0,1:lens)
}   #end else
#build matchdf
nmat <- nrow(matches)
# make starter df
seen<-1*t((matches>0)[,2:3] )
mstat<-cbind(c(seen),c(1,2))
for(x in 1:nmat) {
if (matches[x,1] ==1) mat<-rbind(Pdf[matches[x,2],],Sdf[matches[x,3],])
 else if (matches[x,1] ==2) mat<-rbind(Pdf[matches[x,2],],Pdf[matches[x,2],])
else   mat<- rbind(Sdf[matches[x,3],],Sdf[matches[x,3],])
if (x==1)   matdf<-mat
else matdf<-rbind(matdf,mat)
}  #end for
matchdf <- cbind(mstat,matdf )
names(matchdf)[1] <- "seen"
names(matchdf)[2] <- "PorS"
return(matchdf)
}

create.podsize.calibration.matrix=function(package.dir="C:/Users/Jeff Laake/Desktop/MyDocuments/R Development/ERAnalysis/ERAnalysis/data",save=TRUE)
{
################################################################################
# Estimate pod size calibration matrix (gsS) from pod size calibration data and fitted
# gamma.pod model (assumed to be in workspace) and create additive 
# correction factor vectors and store all as data for package
################################################################################
PodsizeCalibrationData=NULL
PodsizeCalibrationTable=NULL
data(PodsizeCalibrationData,envir=environment()) #,package="ERAnalysis",envir=environment())
data(PodsizeCalibrationTable,envir=environment()) #,package="ERAnalysis",envir=environment())
# Create Reilly additive correction factors
# add.cf.all pools all of the data
# add.cf.reilly only uses 1970s aerial data that was applied in Buckland et al 1992
# add.cf.laake only uses 1990s aerial data that was applied in all surveys from 1992 forward
nmax=20
psdf=PodsizeCalibrationTable
add.cf.all=with(PodsizeCalibrationData,tapply(True-Estimate,cut(Estimate,breaks=c(0,1,2,3,20)),mean))
add.cf.reilly=with(PodsizeCalibrationData[PodsizeCalibrationData$Year=="1978/1979",],tapply(True-Estimate,cut(Estimate,breaks=c(0,1,2,3,20)),mean))
add.cf.laake=with(PodsizeCalibrationData[PodsizeCalibrationData$Type=="Aerial"&PodsizeCalibrationData$Year!="1978/1979",],
       tapply(True-Estimate,cut(Estimate,breaks=c(0,1,2,3,20)),mean))
if(save)save(add.cf.all,file=file.path(package.dir,"add.cf.all.rda"))
if(save)save(add.cf.reilly,file=file.path(package.dir,"add.cf.reilly.rda"))
if(save)save(add.cf.laake,file=file.path(package.dir,"add.cf.laake.rda"))
# Create gsS matrix from pod random effects model with gamma distribution
gss=create.gsS(gamma.pod,nmax=20)
row.names(gsS)=NULL
if(save)save(gsS,file=file.path(package.dir,"gsS.rda"))
return(NULL)
}

create.gsS=function(ps.results,nmax=20,True=NULL)
{
  gam.reint=function(x,obs,xbeta,shape,sigma,nmax=20)
  { 
#   Internal function to be used by integrate function for integration over normal random effect
#  
#   Arguments:
#
#   x    - vector of values passed by integrate
#   obs  - observed size values
#   xbeta- x%*%beta where beta is parameter vector for rate of gamma
#   shape- shape parameters for gamma
#   sigma- value of sigma for normal
#   nmax - maximum pod size
#   
    rate=exp(outer(xbeta,x,"+"))
    return(dnorm(x,sd=sigma)*gam.d(shape,rate,x=rep(obs,length(x)),nmax=nmax))
  }
  gpar=ps.results$par
  sigma=exp(gpar[1])
  if(is.null(True))
    nprime=nmax
  else
    nprime=length(True)
  gsS=matrix(0,nrow=nprime,ncol=nmax)
  for (i in 1:nprime)
  {
     if(is.null(True))
        x=data.frame(size=cut(i,c(1,2,3,4,nmax+1),right=FALSE),plus=as.numeric(i>3),True=i)
     else
        x=data.frame(size=cut(True[i],c(1,2,3,4,nmax+1),right=FALSE),plus=as.numeric(True[i]>3),True=True[i])     
     shape=exp(model.matrix(~size,x)%*%gpar[7:10])
     xbeta=model.matrix(~size+True:plus,x)%*%gpar[2:6]
     for (j in 1:nmax)
       gsS[i,j]=integrate(gam.reint,-5*sigma,5*sigma,xbeta=xbeta,shape=shape,sigma=sigma,obs=j)$val
  }
  return(gsS)
}
################################################################################
reilly.cf=function(x,add.cf) x+add.cf[as.numeric(cut(x,breaks=c(0,1,2,3,20)))]
################################################################################

################################################################################
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
################################################################################

podsize.computations=function(x,gsS,nmax=20,plot=FALSE)
{
################################################################################
# Observed pod size distribution and estimation of true pod size distribution
# assuming equal visibility of all sizes; returns parameters for "true" pod
# size distribution which can be used as starting values for the more
# involved computations of the parameters of true podsize.
################################################################################
# Create observed pod size table and summary of number of pods used
ps.table=table(x$Start.year,factor(x$podsize,levels=1:nmax))
# Plot combined pod size distribution based on proportions
all.years=unique(x$Start.year)
if(plot)
{
  dev.new()
  barplot(ps.table/rowSums(ps.table),beside=TRUE,col=heat.colors(nrow(ps.table)),
         legend.text=all.years)
# Plot observed pod size distributions
  pdf("ObservedPodsizeDistributions.pdf")
  par(mfrow=c(3,2))
  for(i in 1:length(all.years))
     barplot(ps.table[i,],xlab="Pod size",xlim=c(1,nmax),main=all.years[i])
  dev.off()
}
# For each year fit the distribution and plot the fitted/observed distribution
# This assumes that all sizes are equally visible which is not the case
# The fitted models are stored in the list psmod
if(plot)
{
  pdf("fittedps.pdf")
  par(mfrow=c(4,3),mar=c(2,2,1,2))
}
psmod=vector("list",length=length(all.years))
fS=matrix(0,nrow=length(all.years),ncol=nmax)
i=0
for(year in all.years)
{
  i=i+1
  psmod[[i]]=optim(c(0,0),fit.fS,gsS=gsS,pstable=ps.table[i,])
  fS[i,]=gammad(psmod[[i]]$par,nmax)
  if(plot)barplot(rbind(colSums(gsS*fS[i,])*sum(ps.table[i,]),ps.table[i,]),beside=TRUE,legend.text=c("Exp","Obs"),main=paste("Year=",year,"/",year+1,sep=""))
}
# Next plot the estimated "true" distribution of pod size and the distribution of
# measured (observed) pod size.
if(plot)
{
  dev.off()
  pdf("fS.pdf")
  par(mfrow=c(4,3),mar=c(2,2,1,2))
  i=0
  for(year in all.years)
  {
    i=i+1
    barplot(rbind(fS[i,]*sum(ps.table[i,]),ps.table[i,]),beside=TRUE,legend.text=c("True","Estimate"),main=paste("Year=",year,"/",year+1,sep=""))
  }
  dev.off()
}
return(psmod)
}


psfit.gamma <-  function(par,psmat,True=NULL,shape=TRUE,nmax=20)
{
# This function computes the negative log-likelihood for fitting the negative
# binomial to the calibration data  ps=gammad(par,20,True,shape)
  ps=gammad(par,nmax,True,shape)
  if(is.null(True))
    -sum(psmat%*%log(ps))
  else
    -sum(psmat*log(ps))
}

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

negbin.d <-
function(par,nmax,True=NULL)
{
# This function provides normalized probabilities for negative
# binomial for the range from 1 to nmax
#
# Arguments:
#
#  par  - parameters for negative binomial
#  nmax - maximum pod size
#  True - values of true pod size for a range of pod sizes (eg 4+) in which
#          a common distribution is being fitted with a relationship on
#          size with True
#
# The parameters:
#  True = NULL then log(size)=par[1], logit(p)=par[2]
#  True = NULL then log(size)=par[1]+par[2]*True, logit(p)=par[3]
  if(is.null(True))
  {
     pp=diff(c(0,pnbinom(0:(nmax-1),size=exp(par[1]),prob=plogis(par[2]))))
     ps=pp/sum(pp)
  } else
  {
     pp=do.call("rbind",lapply(True, function(tt) pnbinom(0:(nmax-1),size=exp(par[1]+par[2]*tt),prob=plogis(par[3]))))
     ps=t(apply(cbind(rep(0,nrow(pp)),pp),1,diff))/pp[,nmax]
  }
  return(ps)
}

psfit.negbin <-
function(par,psmat,True=NULL,nmax=20)
{
  ps=negbin.d(par,nmax,True)
  if(is.null(True))
    -sum(psmat%*%log(ps))
  else
    -sum(psmat*log(ps))
}

psfit.poisson <-
function(par,psmat,True=NULL,nmax=20)
{
# This function computes the negative log-likelihood for fitting the Poisson
# to the calibration data
  ps=poisson.d(par,nmax,True)
  ps=log(ps)
  ps[is.infinite(ps)]=-200
  if(is.null(True))
    -sum(psmat%*%ps)
  else
    -sum(psmat*ps)
}

poisson.d <-
function(par,nmax,True=NULL)
{
# This function provides normalized probabilities for Poisson
# for the range from 1 to nmax
#
# Arguments:
#
#  par  - parameters for Poisson
#  nmax - maximum pod size
#  True - values of true pod size for a range of pod sizes (eg 4+) in which
#          a common distribution is being fitted with a relationship between
#          True and lambda
#
# The parameters:
#  True = NULL then log(lambda)=par[1]
#  True = NULL then log(lambda)=par[1]+par[2]*True

  if(is.null(True))
  {
     pp=diff(c(0,ppois(0:(nmax-1),lambda=exp(par))))
     ps=pp/sum(pp)
  }
  else
  {
     pp=do.call("rbind",lapply(True, function(tt) ppois(0:(nmax-1),lambda=exp(par[1]+par[2]*tt))))
     ps=t(apply(cbind(rep(0,nrow(pp)),pp),1,diff))/pp[,nmax]
  }
  return(ps)
}

##############################################################################
#  Function from geofunc.xla translated for R and changed to return Kilometers
# instead of NM
#############################################################################
RetDistK<- function (Height, RadPerReticle=0.00497, Reticles)  {
# Height in meters
# RaderReticle = Radians per reticle mark
# Reticles = number of reticles below horizon
# RetDist = distance in kilometers
    if(length(Height)==1)
      Height=rep(Height,length(Reticles))
    x <- sqrt(2 * 6366 * Height / 1000 + (Height / 1000) ^ 2)
    distance=x
    Angle=vector("numeric",length(x))
    whichNA=which(is.na(Reticles))
    Reticles[is.na(Reticles)]=0
    Angle[Reticles>0] <- atan(x[Reticles>0] / 6366)
    distance[Reticles>0] <- ((6366 + Height[Reticles>0] / 1000) * sin(Angle[Reticles>0] + 
         Reticles[Reticles>0]*RadPerReticle) - sqrt(6366 ^ 2 - ((6366 + Height[Reticles>0] / 1000) *
         cos(Angle[Reticles>0] + Reticles[Reticles>0] * RadPerReticle)) ^ 2))
    distance[whichNA]=NA
    return(distance)
}
##############################################################################
#  Function to calculate projected time to cross Azmuth at 241 degrees magnetic
#  from gray whale sighting uses swim speed, time at sighting, and recorded
# bino compass reticle and angle returns decimal hours
# note should add variance
# note speed and variance of speed are from  Swartz et al. 1986
#############################################################################
T241H <- function (stime, sangle, sdist , Speed=6.0)  {
    #stime  time of sighting in decimal hours
    #sangle magnetic angle from bino compass to sighting
    #dist  distance to sighting
    #Speed=6.0 km/hr
  return(stime + ((sin(pi*(sangle-241)/180)*sdist/Speed)))
}
###############################################################################
#  Function to calculate projected distance offshore at 241 degrees magnetic
#  from gray whale sighting uses recorded bino compass angle
# note should add variance
#############################################################################
D241 <- function (sangle, sdist)   {
    #sangle magnetic angle from bino compass to sighting
    #dist  distance to sighting
    #Speed=6.0 km/hr
 return(cos(pi*abs(sangle-241)/180)*sdist)
}




        
		


 
   
   
   