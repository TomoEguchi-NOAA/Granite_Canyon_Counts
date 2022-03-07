

# this script has been splitted into two files: Extract_Data_Upto2016.R and 
# Extract_Data_From2019.R

#library(dplyr)
library(tidyverse)
library(lubridate)



##################################
#This works for years up to 2016
##################################

YEAR <- 2016 #Enter the year of the data files
FILES <- list.files(paste0("Data/", YEAR, "/"))

#FILES <- list.files(paste0(getwd(),"/2016 Edited for JWD"))
#FILES <- list.files(paste0(getwd(),"/GRANITE CANYON 2019_2020_RAW AND EDITED DATA FILES (for Josh Stewart)/Granite Canyon 2020 Visual Data-EDITED"))

for(ff in 1:length(FILES)){ 
  
  # Read in the file, fill blanks with NA, skip the first line (ONLY SKIP IF FILES HAVE 'EDITED FOR' LINE TO BEGIN)
  #data <- read.table(paste0(getwd(),"/2016 Edited for JWD/",FILES[ff]), fill=T, skip=1, na.strings = "", stringsAsFactors = F)
  # data <- read.table(paste0(getwd(),
  #                           "/GRANITE CANYON 2019_2020_RAW AND EDITED DATA FILES (for Josh Stewart)/Granite Canyon 2020 Visual Data-EDITED/",
  #                           FILES[ff]), fill=T, 
  #                    na.strings = "", stringsAsFactors = F)
  
  # TE: the first line to read in starts with "001". So, skip lines until it sees "001"
  all.lines <- read_lines(file = paste0("Data/", YEAR, "/", FILES[ff]))
  idx.line1 <- grep("001", all.lines)
  
  # It does not like to separate the first two fields: line number and event code,
  # which are separated by a space, not tab.
  data <- read.table(paste0("2016 Edited for JWD/", FILES[ff]), 
                     #sep = " ", 
                     fill=T, 
                     skip = (idx.lin1-1),  # eliminates lines up to "001"
                     na.strings = "", 
                     stringsAsFactors = F)
  
  Shifts <- which(data$V2 %in% c('P','E')) #start/end of all shifts
  
  #Shifts <- grep("[PE]", data$V1) #start/end of all shifts
  i <- 1
  for(i in 1:(length(Shifts)-1)){
    #Only use the first observer for model random effect
    Observer <- data[Shifts[i],5] 
    
    # Days since Nov 30th (TE: why 2015?) - each season starts in the previous year
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


##################################
#Modifications for 2019 onwards
##################################
rm(list = ls())
# FILES <- list.files(paste0(getwd(),
#                            "/GRANITE CANYON 2019_2020_RAW AND EDITED DATA FILES (for Josh Stewart)/Granite Canyon 2020 Visual Data-EDITED"))
YEAR <- 2020 #Enter the year of the data files
FILES <- list.files(paste0("Data/", YEAR, "/"))

for(ff in 1:length(FILES)){ 
  
  #Comment lines create formatting problems because each word is interpreted as a column
  #So, first read in individual lines from each file and skip the comment lines when creating a data frame with read.table
  # dataLines <- readLines(paste0(getwd(),
  #                               "/GRANITE CANYON 2019_2020_RAW AND EDITED DATA FILES (for Josh Stewart)/Granite Canyon 2020 Visual Data-EDITED/",FILES[ff]))
  dataLines <- readLines(paste0("Data/", YEAR, "/", FILES[ff]))
  COMMENTS <- which(substr(dataLines,5,5)=="C") #Comment lines, which are problematic for input / formatting
  
  if(length(COMMENTS)>0){
    data <- read.table(text = dataLines[-COMMENTS], 
                       fill=T, 
                       na.strings = "", 
                       stringsAsFactors = F, 
                       col.names = paste0("V",1:16))
  }else{
    data <- read.table(text = dataLines, 
                       fill=T, 
                       na.strings = "", 
                       stringsAsFactors = F, 
                       col.names = paste0("V",1:16))
  }
  # Read in the file, fill blanks with NA, skip the first line (ONLY SKIP IF FILES HAVE 'EDITED FOR' LINE TO BEGIN)
  #data <- read.table(paste0(getwd(),"/2016 Edited for JWD/",FILES[ff]), fill=T, skip=1, na.strings = "", stringsAsFactors = F)
  #data <- read.table(paste0(getwd(),"/GRANITE CANYON 2019_2020_RAW AND EDITED DATA FILES (for Josh Stewart)/Granite Canyon 2020 Visual Data-EDITED/",FILES[ff]), fill=T, na.strings = "", stringsAsFactors = F)
  
  
  #In some cases, an observer accidentally ends a watch and quickly restarts it after recognizing their error
  #To this formatting code, it would look like multiple short shifts, which would get thrown out
  #So we want to identify shift stops and starts within some threshold of each other 
  #I've picked 90 seconds, since it seems like a reasonable amount of time to correct an error
  ### NOTE: THIS COULD GET SQUIRRELY WITH UPDATED DATA / MORE OR LESS ERRORS. CHECK HERE FIRST IF THINGS AREN'T WORKING WITH DATA UPDATE
  Starts <- which(data$V2=="B") #Find all start times
  Ends <- which(data$V2=="E") #Find all end times
  
  if(length(Starts)>0 & length(Ends)>0){
    Diffs <- matrix(NA,ncol=length(Ends),nrow=length(Starts)) #Make an array to hold the time differences of starts and ends
    for(t in 1:length(Starts)){
      Diffs[t,] <- seconds(hms(data[Starts[t],4]))-seconds(hms(data[Ends,4])) #Subtract all end times from each start time (if an end time is less than 90 seconds before a start time, that's probably an error)
    }
    Oops <- which(Diffs >=0 & Diffs<91,arr.ind = T) #Select the differences that are likely errors (End, <90 seconds, Start. Oops!)
    
    if(length(Oops)>0){ #If there are any suspect errors, remove them from the data file
      data <- data[-c(Starts[Oops[1]],(Starts[Oops[1]]+1),Ends[Oops[2]]),]
    }
  }
  #Then carry on as normal:
  
  Shifts <- which(data$V2 %in% c('P','E')) #start/end of all shifts
  
  for(i in 1:(length(Shifts)-1)){
    Observer <- data[Shifts[i],5] #Only use the first observer for model random effect
    BeginDay <- mdy(data[Shifts[i],3]) - mdy(paste0("11/30/",(YEAR-1))) # Days since Nov 30th
    BeginHr <- (hour(hms(data[Shifts[i],4])) + (minute(hms(data[Shifts[i],4]))/60)) # Decimal hour of shift start time
    NextBeginHr <- (hour(hms(data[Shifts[i+1],4])) + (minute(hms(data[Shifts[i+1],4]))/60)) # Decimal hour of next shift start time
    EndHr <- NextBeginHr - 0.00001 # End time is just before next start time (replicating J Durban's calculations)
    Begin <- BeginDay + BeginHr/24 # Beginning time as a decimal day
    End <- BeginDay + (EndHr/24)#End time as a decimal day
    BF <- max(data[Shifts[i]:(Shifts[i+1]-1),12],na.rm=T) #Beaufort (maximum from watch period)
    VS <- max(data[Shifts[i]:(Shifts[i+1]-1),13],na.rm=T) #Visibility (maximum from watch period)
    
    #If there are no sightings, then columns 12/13 won't have BF/VS. So instead, take the VS 
    #record at the start of the shift:
    try(
      {if(!is.numeric(BF)){BF <- data[Shifts[i]+1,5]}
        if(is.na(BF)){BF <- data[Shifts[i]+1,5]}
        if(BF ==-Inf){BF <- data[Shifts[i]+1,5]}
        if(!is.numeric(VS)){VS <- data[Shifts[i]+1,5]}
        if(is.na(VS)){VS <- data[Shifts[i]+1,6]}
        if(VS ==-Inf){VS <- data[Shifts[i]+1,6]}},
      silent=T)
    
    
    
    if(i < (length(Shifts)-1)){ #First set of code only works for shifts prior to the final shift
      GroupsThisWatch <- data[Shifts[i]:(Shifts[i+1]-1),] %>% # Group numbers from this watch period
        filter(V2=="S") %>%
        distinct(V5) %>%
        pull()
      GroupsNextWatch <- data[Shifts[i+1]:(Shifts[i+2]-1),] %>% # Group numbers from next watch period
        filter(V2=="S") %>%
        distinct(V5) %>%
        pull()
      Spillover <- GroupsThisWatch[GroupsThisWatch %in% GroupsNextWatch] # Which groups from watch i were also observed in watch i+1? They should be excluded from i and counted in i+1
      
      if(length(Spillover>0)){ #if there are groups that spill over into following watch, 
        try(
          {N <- data[(Shifts[i]):(Shifts[i+1]-1),] %>% # Calculate the number of whales in the watch period
            filter(V2=="S", !(V5 %in% Spillover), V14!="North") %>% #select sightings only (code S)
            group_by(V5) %>% #group by the whale group number
            summarize(Count = last(V9)) %>% #take the last count of the group (repeated counts tend to increase with more observations)
            summarize(N = sum(as.numeric(Count)))}, #sum all final group counts
          silent = T)
      }else{ #if there aren't groups that spill over into following watch, the spillover argument is removed from filter() to avoid errors
        try(
          {N <- data[(Shifts[i]):(Shifts[i+1]-1),] %>% # Calculate the number of whales in the watch period
            filter(V2=="S", V14!="North") %>% #select sightings only (code S)
            group_by(V5) %>% #group by the whale group number
            summarize(Count = last(V9)) %>% #take the last count of the group (repeated counts tend to increase with more observations)
            summarize(N = sum(as.numeric(Count)))}, #sum all max group counts
          silent = T)
      }#ifelse
    }else{#if i > length(Shifts)-1 (In the last shift, there's obviously no spillover to the next shift. And the code creates errors if you try to do that. So this is a separate snippet for just the last shift)
      try( #If there were no sightings, this code won't work. So wrap it in try() to avoid an error stoppage
        {N <- data[(Shifts[i]):(Shifts[i+1]-1),] %>% # Calculate the number of whales in the watch period
          filter(V2=="S", V14!="North") %>% #select sightings only (code S)
          group_by(V5) %>% #group by the whale group number
          summarize(Count = last(V9)) %>% #take the last count of the group (repeated counts tend to increase with more observations)
          summarize(N = sum(as.numeric(Count)))}, #sum all max group counts
        silent = T)
    }#ifelse2
    
    if(exists("N")==F){N <- 0} #If there were no sightings, the above code will not assign an N. If N doesn't exist, set it to 0
    
    if(i == 1){
      Output <- data.frame(begin=as.numeric(Begin),
                           end=as.numeric(End),
                           bf=as.numeric(BF),
                           vs=as.numeric(VS),
                           n=as.numeric(N),
                           obs=as.character(Observer),
                           ff=ff, #ff and i for easy file reference
                           i=i,
                           BeginHr = BeginHr,
                           BeginDay = BeginDay)
    }else{
      Output <- rbind(Output,
                      data.frame(begin=as.numeric(Begin),
                                 end=as.numeric(End),
                                 bf=as.numeric(BF),
                                 vs=as.numeric(VS),
                                 n=as.numeric(N),
                                 obs=as.character(Observer),
                                 ff=ff,
                                 i=i,
                                 BeginHr=BeginHr,
                                 BeginDay = BeginDay))
      
      #Output <- filter(Output, bf!=5, vs!=5) #Remove any watches where beaufort or vis was > 4
    }#else output
    
  }#i
  
  if(ff==1){
    Data_Out <- Output
  }else{
    Data_Out <- rbind(Data_Out,Output)  
  }  
  
}#ff (files loop)         


View(Data_Out[which(Data_Out$end-Data_Out$begin < 0.059),]) #Entries that are less than 1.5hrs (5 minute grace period)

View(Data_Out[which(Data_Out$end-Data_Out$begin > 0.059),]) #Entries more than 1.5hrs (5 minute grace period)
#0.0625 is the exact watch period, but I've given some leeway. If it's less than 5 minutes short, I'm counting it
#If it's less than 10 minutes over the 1.5hrs, I'm also counting it (guessing that they forgot to log off or something)

View(Data_Out[which(Data_Out$end-Data_Out$begin > 0.0695),]) #Entries more than 1.5 hrs + 10 min

View(Data_Out[which(Data_Out$end-Data_Out$begin > 0.066),])#Entries more than 1.5hrs + 5 min
View(Data_Out[which(Data_Out$end-Data_Out$begin > 0.0625),])#Entries more than 1.5hrs


# Final data
# Remove watches that were less than 85 minutes or greater than 95 minutes:

CorrectLength <- Data_Out[which(Data_Out$end-Data_Out$begin > 0.059 & Data_Out$end-Data_Out$begin < 0.066),]
Chaff <- Data_Out[which(Data_Out$end-Data_Out$begin < 0.059),]


FinalData <- CorrectLength %>%
  filter(bf < 5, vs < 5)

WhalesDays <- FinalData %>%
  group_by(BeginDay) %>%
  mutate(PropDay = end-begin) %>%
  summarize(TotalWatch = sum(PropDay), TotalWhales=sum(n))

plot(x=WhalesDays$BeginDay,y=WhalesDays$TotalWhales/WhalesDays$TotalWatch)

ShiftsPerDay <- FinalData %>%
  group_by(BeginDay) %>%
  summarize(Watches = n())


#Spot Checks:

#Check shifts that passed muster to confirm the compiled data is correct
set.seed(1199)
View(FinalData[sample(1:196,20,replace=F),])

#Check shifts that were thrown out to make sure they deserved it
set.seed(1200)
View(Chaff[sample(1:54,10,replace=F),])


#Summary Stats for Report
Complete_Data <- Data_Out[complete.cases(Data_Out),]
Complete_Data$Eff <- Complete_Data$end-Complete_Data$begin

TotalHrs <- sum(Complete_Data$Eff)*24
TotalDays <- length(unique(floor(Complete_Data$begin)))
TotalObservers <- length(unique(Complete_Data$obs))
TotalWhales <- sum(Complete_Data$n)


WPH <- Complete_Data %>%
  group_by(floor(Complete_Data$begin)) %>%
  summarize(TotalWhales = sum(n), TotalEffort = sum(Eff), WPH = sum(n)/(sum(Eff)*24)) 

out.obj <- list(WPH = WPH,
                FinalData = FinalData,
                Data_Out = Data_Out,
                WhalesDays = WhalesDays,
                Complete_Data = Complete_Data)
saveRDS(out.obj, file = paste0("RData/out_", YEAR, "_Joshs.rds"))

