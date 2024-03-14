
#I used the ERAnalysis package (https://github.com/jlaake/ERAnalysis) to obtain code and data (1967/1968 - 2006/2007). Data for recent years (2009/2010 - 2022/2023) were arranged to the input format and analysis ran.

rm(list = ls())

# from LaakeAnalysis_NewData.R
library(tidyverse)
library(ERAnalysis)


# Function from Laake's code
conf.int=function(abundance, CV, alpha=0.05, digits=2, prt=FALSE){
  # Computes confidence intervals based on lognormal distr.
  # JMB / NMML / 11 Sep 2008
  
  if (alpha <0 || alpha > .999) stop("alpha must be in (0,1)")
  z = round(abs(qnorm(alpha/2)),2)
  if (prt) cat("N:",abundance,"  cv:",CV,"  alpha:",alpha,"  z:",z,"\n")
  C <- exp(z * sqrt(log(1 + CV^2)))
  SL <- round(abundance/C,digits)
  SU <- round(abundance * C,digits)
  data.frame(SL,SU)
}

# Observers need to be in integers, not their initials.
data("Observer")   # from ERAnalysis package.
#Observer$Set = "old"

new.observers <- read.csv(file = "Data/ObserverList2023.csv") %>%
  transmute(ID = ID,
            Initials = obs)

Observer %>%
  left_join(new.observers, by = "Initials") %>%
  filter(!is.na(ID.y)) %>%
  transmute(ID = ID.y,
            Initials = Initials) -> tmp

new.observers %>%
  anti_join(tmp, by = "ID") -> tmp.2

tmp.2$new.ID = 68:((68+nrow(tmp.2))-1)

tmp.2 %>%
  select(new.ID, Initials) %>%
  mutate(ID = new.ID,
         Initials = Initials) %>%
  select(-new.ID) -> tmp.3

all.observers <- rbind(Observer %>% select(ID, Initials), tmp.3) %>%
  mutate(Observer = Initials) %>%
  na.omit() %>%
  filter(ID != "21") %>%  # ID - 21 is ARV/AVS, which is not useful
  droplevels()

# add those initials back in
all.observers <- rbind(all.observers, 
                       data.frame(ID = c("21", "21"), 
                                  Initials = c("ARV", "AVS"), 
                                  Observer = c("ARV", "AVS")))

# sightings
# data(PrimarySightings)
# data("SecondarySightings")
# 
# data(PrimaryEffort)

Laake_PrimarySightings <- read.csv(file = "Data/Laake_PrimarySightings.csv") %>%
  left_join(all.observers, by = "Observer") %>%
  select(-Initials) %>%
  mutate(ID.1 = ID,
         ID = Observer) %>%
  select(-Observer) %>%
  left_join(all.observers, by = "ID") %>%
  mutate(ID.new = ifelse(is.na(ID.1), ID, ID.1)) %>%
  select(-c(ID, ID.1, Initials, Observer)) %>%
  mutate(ID = ID.new) %>%
  select(-ID.new) #-> tmp

Laake_SecondarySightings <- read.csv(file = "Data/Laake_SecondarySightings.csv") 

# effort
Laake_PrimaryEffort <- read.csv(file = "Data/Laake_PrimaryEffort.csv") %>%
  filter(Use) %>%
  left_join(all.observers, by = "Observer") %>%
  select(-Initials) %>%
  mutate(ID.1 = ID,
         ID = Observer) %>%
  select(-Observer) %>%
  left_join(all.observers, by = "ID") %>%
  mutate(ID.new = ifelse(is.na(ID.1), ID, ID.1)) %>%
  select(-c(ID, ID.1, Initials, Observer)) %>%
  mutate(ID = ID.new) %>%
  select(-ID.new) #-> tmp
  
Laake_SecondaryEffort <- read.csv(file = "Data/Laake_SecondaryEffort.csv") %>%
  filter(Use) %>%
  left_join(all.observers, by = "Observer") %>%
  select(-Initials) %>%
  mutate(ID.1 = ID,
         ID = Observer) %>%
  select(-Observer) %>%
  left_join(all.observers, by = "ID") %>%
  mutate(ID.new = ifelse(is.na(ID.1), ID, ID.1)) %>%
  select(-c(ID, ID.1, Initials, Observer)) %>%
  mutate(ID = ID.new) %>%
  select(-ID.new) #-> tmp

# Only 2010 and 2011 had two observation stations, which are needed for applying
# Laake's method beyond naive estimates
years <- c(2010, 2011, 2015, 2016, 2020, 2022, 2023)

# sightings and effort
sightings.list.primary <- effort.list.primary  <- list()
sightings.list.secondary <- effort.list.secondary  <- list()

for (k in 1:length(years)){
  # These raw data files contain repeated observations of all groups.
   
  tmp.sightings <- read.csv(paste0("Data/all_sightings_", 
                                   years[k], "_Tomo_v2.csv")) 
  
  tmp.sightings %>%
    mutate(Date1 = as.Date(Date, format = "%m/%d/%Y")) %>%
    #group_by(group) %>%
    transmute(Date = Date1,
              Time = Time_PST, 
              day = day(Date),
              month = month(Date),
              year = year(Date),
              watch = shift,
              t241 = difftime((paste(Date, Time)),
                              (paste0((years[k] - 1), 
                                             "-11-30 00:00:00"))) %>%
                as.numeric(),
              Group_ID = Group_ID,
              distance = Distance,
              podsize = nwhales,
              vis = vis,
              beaufort = beaufort,
              #wind.direction = NA,
              key = key,
              pphr = pphr,
              Start.year = years[k]-1,
              #original.watch = NA,
              only = TRUE,
              #hours = NA,
              #Sex = NA,
              Observer = toupper(observer),
              station = station) %>%
    arrange(Date, Group_ID) %>%
    filter(vis < 5, beaufort < 5) %>%
    select(-c(Group_ID, Time)) %>%
    na.omit() -> sightings.all

  sightings.all %>% 
    filter(station == "P") %>%
    select(-station) -> sightings.list.primary[[k]]
  
  sightings.all %>% 
    filter(station == "S") %>%
    select(-station) -> sightings.list.secondary[[k]]
  
  tmp.effort <- read.csv(paste0("Data/all_effort_", 
                                years[k], "_Tomo_v2.csv")) 
  
  tmp.effort %>%
    transmute(watch.key = watch.key,
              Start.year = years[k] - 1,
              key = key,
              begin = begin,
              end = end,
              npods = npods,
              nwhales = nwhales,
              effort = effort,
              vis = vis,
              beaufort = beaufort,
              Observer = toupper(observer),
              time = time,
              watch = shift,
              date.shift = date.shift,
              Use = T,
              Date = as.Date(Date, format = "%m/%d/%Y"),
              station = station) %>%
    filter(vis < 5, beaufort < 5) %>%
    na.omit() -> effort.all
  
  effort.all %>%
    filter(station == "P") -> effort.list.primary[[k]]
  
  effort.all %>%
    filter(station == "S") -> effort.list.secondary[[k]]
}

# sightings
sightings.primary <- do.call("rbind", sightings.list.primary) %>%
  na.omit() %>%
  transmute(X = row_number(),
            Date = as.character(Date),
            day = day,
            month = month,
            year = year,
            watch = watch,
            t241 = t241,
            distance = distance,
            podsize = podsize,
            vis = vis,
            beaufort = beaufort,
            wind.direction = NA,
            key = key,
            pphr = pphr,
            Start.year = Start.year,
            original.watch = NA,
            only = TRUE,
            hours = NA,
            Sex = NA,
            ID = Observer)
              
# sightings.primary %>% 
#   left_join(all.observers, by = "Observer") %>%
#   select(-c(Observer, Initials)) %>%
#   rename(Observer = ID) -> sightings.primary

# combine with Laake's
sightings.primary.all <- rbind(Laake_PrimarySightings, sightings.primary)

sightings.secondary <- do.call("rbind", sightings.list.secondary) %>%
  na.omit()%>%
  transmute(X = row_number(),
            #Date = Date,
            day = day,
            month = month,
            year = year,
            etime = NA,
            #watch = watch,
            t241 = t241,
            distance = distance,
            podsize = podsize,
            vis = vis,
            beaufort = beaufort,
            wind.direction = NA,
            off = NA,
            Start.year = Start.year,
            date = as.character(Date))

# sightings.secondary %>% 
#   left_join(all.observers, by = "Observer") %>%
#   select(-c(Observer, Initials)) %>%
#   rename(Observer = ID) -> sightings.secondary

# combine with Laake's
sightings.secondary.all <- rbind(Laake_SecondarySightings, sightings.secondary)

# Effort
effort.primary <- do.call("rbind", effort.list.primary)  %>%
  na.omit()

effort.primary %>% 
  left_join(all.observers, by = "Observer") %>%
  #select(-c(Observer, Initials)) %>%
  #rename(Observer = ID) %>%
  transmute(X = row_number(),
            watch.key = watch.key,
            Start.year = Start.year,
            key = key,
            begin = begin,
            end = end,
            npods = npods,
            nwhales = nwhales,
            effort = effort,
            vis = vis,
            beaufort = beaufort,
            time = time,
            watch = watch,
            Use = Use,
            Date = as.character(Date), 
            ID = ID) -> effort.primary

effort.primary.all <- rbind(Laake_PrimaryEffort, effort.primary)

effort.secondary <- do.call("rbind", effort.list.secondary)  %>%
  na.omit()

effort.secondary %>% 
  left_join(all.observers, by = "Observer") %>%
  transmute(X = row_number(),
            watch.key = watch.key,
            Start.year = Start.year,
            key = key,
            begin = begin,
            end = end,
            npods = npods,
            nwhales = nwhales,
            effort = effort,
            vis = vis,
            beaufort = beaufort,
            time = time,
            watch = watch,
            Use = Use,
            Date = as.character(Date),
            ID = ID) -> effort.secondary

effort.secondary.all <- rbind(Laake_SecondaryEffort, effort.secondary)

# gsS: nmax x nmax pod size calibration matrix; each row is a true pod size 
# from 1 to nmax and the value for each column is the probability that a pod of 
# a true size S is recorded as a size s (1..nmax columns)
# 

# Compute naive estimates of abundance for all surveys.  These use the uncorrected
# counts of whales from the primary observer during watches in which neither Beaufort nor
# vis exceeded 4.  For each year a gam with a smooth over time is fitted and this is
# used to predict total abundance throughout the migration from the counts of whales
# during the sampled periods.  There is no correction for missed pods or for measurement
# error in podsize. Each fitted migration gam is plotted with the observed values and
# saved in the file NaiveMigration.pdf.

final.time <- sapply(tapply(floor(effort.primary.all$time),
                            effort.primary.all$Start.year,max), function(x) ifelse(x>90,100,90))

lower.time <- rep(0,length(final.time))

all.start.years <- unique(sightings.primary.all$Start.year)

naive.abundance.models <- vector("list", length(all.start.years))
i=0
for (year in all.start.years){
  i=i+1
  primary=sightings.primary.all[sightings.primary.all$Start.year==year,]
  primary$Start.year=factor(primary$Start.year)
  ern=subset(effort.primary.all,
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

# Next compute the series of abundance estimates for the most recent 8 years (1987/1988 - 2006/2007) by
# fitting and selecting the best detection model.  From those 8 estimates and the
# naive estimates, compute an average ratio and apply it to generate the estimates
# for the first 15 surveys prior to 1987. Note with hessian=TRUE, the analysis can
# take about 30-60 minutes to complete. (TE: Takes about 6 minutes now. But added
# the if-else. 2023-08-31) This was run separately and results saved.

abundance.estimates <- readRDS("RData/Laake_abundance_estimates.rds")

ratio <- abundance.estimates$ratio
ratio.SE <- 0.03  # from Laake et al 2012, Table 8

# Compute series of estimates for before 1987 without nighttime correction factor (eqn 24)  
W.hat.1 <- c(sapply(naive.abundance.models[1:15], 
                    function(x)x$Total)*ratio)

# Apply nighttime correction factor (eqn 29)
fn = 1.0817
SE.fn <- 0.0338

# Need to add CI or SE and var-cov matrix. 
# Bring in the results from Laake_example_code.R
abundance.vc <- read_rds(file = "RData/abundance.vc.rds")

# Var(W.hat) for year < 1987
W.tilde.1 <- c(sapply(naive.abundance.models[1:15], function(x) x$Total))  # naive abundance
var.W.tilde.1 <- c(sapply(naive.abundance.models[1:15], function(x) x$var.Total))
var.W.hat.1 <- W.tilde.1^2 * ratio.SE^2 * 9 + ratio^2 * var.W.tilde.1   # eqn 27

N.hat.1 <- W.hat.1 * fn

# SE values are a little different from what I calcualted above (SE.Nhat.1) but not much
SE.Nhat.1 <- abundance.vc$se[1:length(N.hat.1)]

# var(W.hat) for year > 1985 eqn. 25
# From Table 8 in Laake et al. 2012
W.hat.2 <- setNames(c(24883, 14571, 18585, 19362, 19539, 15133, 14822, 17682),
                    c("1987", "1992", "1993", "1995", "1997", "2000", "2001", "2006"))

#W.hat <- c(W.hat.1, W.hat.2)
N.hat.2 <- abundance.estimates$summary.df$Nhat
SE.Nhat.2 <- abundance.vc$se[(length(N.hat.1)+1):length(abundance.vc$se)]

# The same approach for year < 1987 will be used for years 2009 - 2022
# Although the values didn't match exactly, they were quite close. So, 
# I'm not going to worry too much about it.
W.tilde.3 <- c(sapply(naive.abundance.models, function(x) x$Total))
var.W.tilde.3 <- c(sapply(naive.abundance.models, function(x) x$var.Total))
W.hat.3 <- c(sapply(naive.abundance.models, 
                    function(x)x$Total)*ratio)

var.W.hat.3 <- W.tilde.3^2 * ratio.SE^2 * 9 + ratio^2 * var.W.tilde.3   # eqn 27
N.hat.3 <- W.hat.3 * fn

# Fix the following three lines according to what I find on lines 721-724
var.Nhat.3 <- (fn * W.hat.3)^2 * ((SE.fn/fn)^2 + (var.W.hat.3/((W.hat.3)^2)))  # eqn 30
SE.Nhat.3 <- sqrt(var.Nhat.3)

all.estimates.Laake <- data.frame(Year = all.start.years, 
                            Nhat = N.hat.3,
                            SE = SE.Nhat.3) %>%
  mutate(CV = SE/Nhat,
         Season = paste0(Year, "/", (Year + 1)))

CI <- conf.int(all.estimates.Laake$Nhat, all.estimates.Laake$CV)

all.estimates.Laake$LCL <- CI$SL
all.estimates.Laake$UCL <- CI$SU
all.estimates.Laake$Method <- "Laake"

