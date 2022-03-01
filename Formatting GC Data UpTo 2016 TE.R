
rm(list = ls())
#library(dplyr)
library(tidyverse)
library(lubridate)


##################################
#This works for years up to 2016
##################################

YEAR <- 2016
FILES <- list.files(paste0("Data/", YEAR, "/"))
#FILES <- list.files(paste0(getwd(),"/2016 Edited for JWD"))
#FILES <- list.files(paste0(getwd(),"/GRANITE CANYON 2019_2020_RAW AND EDITED DATA FILES (for Josh Stewart)/Granite Canyon 2020 Visual Data-EDITED"))

ff <- 1
for(ff in 1:length(FILES)){ 
  
  # Read in the file, fill blanks with NA, skip the first line (ONLY SKIP IF FILES HAVE 'EDITED FOR' LINE TO BEGIN)
  #data <- read.table(paste0(getwd(),"/2016 Edited for JWD/",FILES[ff]), fill=T, skip=1, na.strings = "", stringsAsFactors = F)
  # data <- read.table(paste0(getwd(),
  #                           "/GRANITE CANYON 2019_2020_RAW AND EDITED DATA FILES (for Josh Stewart)/Granite Canyon 2020 Visual Data-EDITED/",
  #                           FILES[ff]), fill=T, 
  #                    na.strings = "", stringsAsFactors = F)
  
  # TE: the first line to read in starts with "001". So, skip lines until it sees "001"
  
  
  # It does not like to separate the first two fields: line number and event code,
  # which are separated by a space, not tab.
  data <- read.table(paste0("Data/", YEAR, "/", FILES[ff]), 
                     #sep = " ", 
                     fill=T, 
                     skip = 1,  # eliminates "EIDTED for JWD"
                     na.strings = "", 
                     stringsAsFactors = F)
  
  Shifts <- which(data$V2 %in% c('P','E')) #start/end of all shifts
  
  #Shifts <- grep("[PE]", data$V1) #start/end of all shifts
  i <- 1
  for(i in 1:(length(Shifts)-1)){
    #Only use the first observer for model random effect
    Observer <- data[Shifts[i],5] 
    # Days since Nov 30th (TE: why 2015? - because the study period starts in the previous year)
    BeginDay <- mdy(data[Shifts[i],3]) - mdy(paste0("11/30/", (YEAR-1)))
    
    # Decimal hour of shift start time
    BeginHr <- (hour(hms(data[Shifts[i],4])) + 
                  (minute(hms(data[Shifts[i],4]))/60)) 
    
    # Decimal hour of next shift start time
    NextBeginHr <- (hour(hms(data[Shifts[i+1],4])) + 
                      (minute(hms(data[Shifts[i+1],4]))/60)) 
    # End time is just before next start time (replicating J Durban's calculations)
    EndHr <- NextBeginHr - 0.00001 
    # Beginning time as a decimal day
    Begin <- BeginDay + BeginHr/24 
    #End time as a decimal day
    End <- BeginDay + (EndHr/24)
    #Beaufort (maximum from watch period)
    BF <- max(data[Shifts[i]:(Shifts[i+1]-1),12],na.rm=T) 
    #Visibility (maximum from watch period)
    VS <- max(data[Shifts[i]:(Shifts[i+1]-1),13],na.rm=T) 
    
    if(BF ==-Inf){BF <- data[Shifts[i]+1,5]}
    if(VS ==-Inf){VS <- data[Shifts[i]+1,6]}
    #First set of code only works for shifts prior to the final shift
    if(i < (length(Shifts)-1)){ 
      # Group numbers from this watch period
      GroupsThisWatch <- data[Shifts[i]:(Shifts[i+1]-1),] %>% 
        filter(V2=="S") %>%
        distinct(V5) %>%
        pull()
      
      # Group numbers from next watch period
      GroupsNextWatch <- data[Shifts[i+1]:(Shifts[i+2]-1),] %>% 
        filter(V2=="S") %>%
        distinct(V5) %>%
        pull()
      # Which groups from watch i were also observed in watch i+1? 
      # They should be excluded from i and counted in i+1
      Spillover <- GroupsThisWatch[GroupsThisWatch %in% GroupsNextWatch] 
      
      #if there are groups that spill over into following watch, 
      if(length(Spillover>0)){ 
        # Calculate the number of whales in the watch period
        N <- data[(Shifts[i]):(Shifts[i+1]-1),] %>% 
          filter(V2=="S", !(V5 %in% Spillover), V14!="North") %>% #select sightings only (code S)
          group_by(V5) %>% #group by the whale group number
          summarize(Count = last(V9)) %>% #take the last count of the group (repeated counts tend to increase with more observations)
          summarize(N = sum(as.numeric(Count))) #sum all final group counts
        #if there aren't groups that spill over into following watch, 
        # the spillover argument is removed from filter() to avoid errors  
      }else{ 
        # Calculate the number of whales in the watch period
        N <- data[(Shifts[i]):(Shifts[i+1]-1),] %>% 
          filter(V2=="S", V14!="North") %>% #select sightings only (code S)
          group_by(V5) %>% #group by the whale group number
          summarize(Count = last(V9)) %>% #take the last count of the group (repeated counts tend to increase with more observations)
          summarize(N = sum(as.numeric(Count))) #sum all max group counts
      }#ifelse
    }else{#if i > length(Shifts)-1 (In the last shift, there's obviously no spillover to the next shift. And the code creates errors if you try to do that. So this is a separate snippet for just the last shift)
      N <- data[(Shifts[i]):(Shifts[i+1]-1),] %>% # Calculate the number of whales in the watch period
        filter(V2=="S", V14!="North") %>% #select sightings only (code S)
        group_by(V5) %>% #group by the whale group number
        summarize(Count = last(V9)) %>% #take the last count of the group (repeated counts tend to increase with more observations)
        summarize(N = sum(as.numeric(Count))) #sum all max group counts
    }#ifelse2
    
    if(i == 1){
      Output <- data.frame(begin=as.numeric(Begin),
                           end=as.numeric(End),
                           bf=as.numeric(BF),
                           vs=as.numeric(VS),
                           n=as.numeric(N),
                           obs=as.character(Observer))
    }else{
      Output <- rbind(Output,
                      data.frame(begin=as.numeric(Begin),
                                 end=as.numeric(End),
                                 bf=as.numeric(BF),
                                 vs=as.numeric(VS),
                                 n=as.numeric(N),
                                 obs=as.character(Observer)))
      #Remove any watches where beaufort or vis was > 4
      Output <- filter(Output, bf!=5, vs!=5) 
    }#else output
    
  }#i
  
  if(ff==1){
    Data_Out <- Output
  }else{
    Data_Out <- rbind(Data_Out,Output)  
  }  
  
}#ff (files loop)         


Data_Out[which(Data_Out$end-Data_Out$begin < 0.059),]

dplyr::filter(data, V2=='P')
is.na(as.numeric(data$V5))


