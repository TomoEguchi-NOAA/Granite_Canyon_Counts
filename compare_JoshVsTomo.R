
rm(list=ls())

get.data <- function(dir, YEAR, ff){
  FILES <- list.files(paste0(dir, "/", YEAR))
  all.lines <- read_lines(file = paste0(dir, "/", YEAR, "/", FILES[ff]))
  event.code <- str_sub(all.lines, start = 5, end = 5)
  COMMENTS <- which(event.code == "C")
  
  if(length(COMMENTS)>0){
    data <- read.table(text = all.lines[-COMMENTS],
                       fill=T,
                       na.strings = "",
                       stringsAsFactors = F,
                       col.names = paste0("V",1:16))
  }else{
    data <- read.table(text = all.lines,
                       fill=T,
                       na.strings = "",
                       stringsAsFactors = F,
                       col.names = paste0("V",1:16))
    
  }
  
  # Files with extensive comments in sightings are problematic because they get split into multiple lines. 
  # We need to pull out lines that contain only numeric V1 (when they are converted into numeric)
  
  data <- data[!is.na(as.numeric(data$V1)),]
  
  
  Starts <- which(data$V2=="B") #Find all start times
  Ends <- which(data$V2=="E") #Find all end times
  
  if(length(Starts)>0 & length(Ends)>0){
    #Make an array to hold the time differences of starts and ends
    Diffs <- matrix(NA, ncol=length(Ends), nrow=length(Starts)) 
    
    for(t in 1:length(Starts)){
      #Subtract all end times from each start time (if an end time is less than 90 seconds 
      # before a start time, that's probably an error)
      # TE: I added [t] to Ends in the following line. I think it's needed. NO... 
      # Ends does not need the subscript. 
      Diffs[t,] <- seconds(hms(data[Starts[t],4])) - seconds(hms(data[Ends,4])) 
    }
    
    #Select the differences that are likely errors (End, <90 seconds, Start. Oops!)
    Oops <- which(Diffs >=0 & Diffs<91,arr.ind = T) 
    
    if(length(Oops)>0){ #If there are any suspect errors, remove them from the data file
      data <- data[-c(Starts[Oops[1]], (Starts[Oops[1]]+1), Ends[Oops[2]]),]
    }
  }
  
  BeginDay <- mdy(data$V3) - mdy(paste0("11/30/", (YEAR - 1)))
  # Decimal hour of shift start time
  BeginHr <- hour(hms(data$V4)) + minute(hms(data$V4))/60
  
  data %>% mutate(begin = as.numeric(BeginDay) + BeginHr/24) -> data
  return(data)
}

get.shift <- function(YEAR, data, i){
  Shifts.begin <- which(data$V2 %in% "P") 
  Shifts.end <- c((Shifts.begin[2:length(Shifts.begin)] - 1), which(data$V2 %in% "E"))
  #Shifts.end <- which(data$V2 %in% "E")

  max.shifts <- length(Shifts.begin)
  #Only use the first observer for model random effect
  Observer <- data[Shifts.begin[i], 5] 
  # Days since Nov 30th of the previous year
  BeginDay <- mdy(data[Shifts.begin[i], 3]) - mdy(paste0("11/30/", (YEAR - 1)))
  # Decimal hour of shift start time
  BeginHr <- (hour(hms(data[Shifts.begin[i], 4])) + (minute(hms(data[Shifts.begin[i], 4]))/60)) 
  # Decimal hour of next shift start time
  if (i < max.shifts){
    NextBeginHr <- (hour(hms(data[Shifts.begin[i+1], 4])) + (minute(hms(data[Shifts.begin[i+1], 4]))/60)) 
  } else {
    NextBeginHr <- (hour(hms(data[which(data$V2 %in% "E"), 4])) + (minute(hms(data[which(data$V2 %in% "E"), 4]))/60)) + 0.00001
  }
  # End time is just before next start time (replicating J Durban's calculations)
  EndHr <- NextBeginHr - 0.00001 
  # Beginning time as a decimal day
  Begin <- as.numeric(BeginDay) + (BeginHr/24) 
  #End time as a decimal day
  End <- as.numeric(BeginDay) + (EndHr/24)
  BFs <- as.numeric(data[Shifts.begin[i]:Shifts.end[i], 12])
  if (sum(!is.na(BFs)) == 0){
    BF <- NA
  } else {
    BF <- max(BFs, na.rm=T)
  }
  
  #Visibility (maximum from watch period)
  VSs <- as.numeric(data[Shifts.begin[i]:Shifts.end[i], 13])
  if (sum(!is.na(VSs)) == 0){
    VS <- NA 
  } else {
    VS <- max(VSs, na.rm=T) 
  }
  
  # if still NA
  if (is.na(BF)) {BF <- data[Shifts.begin[i]+1, 5]}
  if (is.na(VS)) {VS <- data[Shifts.begin[i]+1, 6]}

  Spillover <- vector(length = 0)
  if (i < max.shifts){
    # Groups = Observers. Only the first (primary) observer is considered (V5)
    GroupsThisWatch <- data[Shifts.begin[i]:Shifts.end[i],] %>% # Group numbers from this watch period
      filter(V2 == "S") %>%
      distinct(V5) %>%
      pull()
    
    GroupsNextWatch <- data[Shifts.begin[i+1]:Shifts.end[i+1],] %>% # Group numbers from next watch period
      filter(V2 == "S") %>%
      distinct(V5) %>%
      pull()
    
    # Which groups from watch i were also observed in watch i+1? They should be excluded from i and counted in i+1
    Spillover <- GroupsThisWatch[GroupsThisWatch %in% GroupsNextWatch] 
    
  }
  
    
  if(length(Spillover > 0)){ #if there are groups that spill over into following watch, 
    # figure out if there were any sightings that need to be considered:
    sub.data <- data[(Shifts.begin[i]):(Shifts.end[i]),] %>% 
      filter(V2 == "S", !(V5 %in% Spillover), V14 != "North")
    
    if (nrow(sub.data) > 0){
      N <- sub.data %>%
        group_by(V5) %>% #group by the whale group number
        select(V5, V9) %>%
        summarize(N = max(as.numeric(V9), na.rm = T)) %>% 
        select(N)  %>% sum()
    } else {
      N <- 0
    }
    
  } else {   # if there were no spillover
    sub.data <- data[Shifts.begin[i]:Shifts.end[i],]  %>%  
      filter(V2 == "S", V14 != "North")
      
    if (nrow(sub.data) > 0){
      N <- sub.data %>%
        group_by(V5) %>% #group by the whale group number
        select(V5, V9) %>%
        summarize(N = max(as.numeric(V9), na.rm = T)) %>% 
        select(N)  %>% sum()
    } else {
      N <- 0
    }
    
  }
  
  out.list <- list(out.df = data.frame(begin=as.numeric(Begin),
                                       end=as.numeric(End),
                                       dur = as.numeric(End) - as.numeric(Begin),
                                       bf=as.numeric(BF),
                                       vs=as.numeric(VS),
                                       n = N,
                                       obs=as.character(Observer),
                                       i=i,
                                       BeginHr = BeginHr,
                                       BeginDay = BeginDay),
                   data = sub.data)
  return( out.list )
}


Tomo.out <- readRDS("RData/out_2020_Tomos.rds")
Josh.out <- readRDS("RData/out_2020_Joshs.rds")
FinalData.Josh <- Josh.out$FinalData %>% mutate(v = "J")
FinalData.Tomo <- Tomo.out$FinalData %>% mutate(v = "T")
FinalData.Both <- rbind(FinalData.Josh, FinalData.Tomo)

dim(FinalData.Josh)
dim(FinalData.Tomo)

FinalData.Josh %>% 
  group_by(ff) %>% 
  summarize(nrow = n()) -> Josh.summary

FinalData.Tomo %>% 
  group_by(ff) %>% 
  summarize(nrow = n()) -> Tomo.summary

Tomo.summary %>% left_join(Josh.summary, by = "ff") %>%
  mutate(dif = nrow.x - nrow.y) -> TomoVsJosh_n

#%>% filter(abs(dif) > 0) 
comp.df <- data.frame(nrow = nrow(Tomo.out$FinalData), ncol = 6) 
#begin = NA, end = NA, max.bf = NA, max.vs = NA, total.n = NA)
Fs <- unique(FinalData.Tomo$ff)
for (f in 1:length(Fs)){
  J1 <- Josh.out$FinalData %>% filter(ff == Fs[f])
  T1 <- Tomo.out$FinalData %>% filter(ff == Fs[f])
  
  comp.df[f,1] <- J1$begin[1] - T1$begin[1]
  comp.df[f,2] <- J1$end[1] - T1$end[1]
  comp.df[f,3] <- max(J1$bf) - max(T1$bf)
  comp.df[f,4] <- max(J1$vs) - max(T1$vs)
  comp.df[f,5] <- sum(J1$n) - sum(T1$n)
  comp.df[f,6] <- Fs[f]

}

# The following discrepancies for 2020 survey need to be examined

###########################################################################
FinalData.Both %>% filter(ff == 23) %>% filter(v == "J")
FinalData.Both %>% filter(ff == 23) %>% filter(v == "T")
# i = 5 contains 2 more whales (16) in Tomo's than in Josh's (14)

# pull out the data, examine the spillover groups
data.23 <- get.data("Data/", 2020, ff = 23)
output.23.5 <- get.shift(2020, data.23, i = 5)
output.23.6 <- get.shift(2020, data.23, i = 6)

output.23.5$data %>% distinct(V5)

data.23 %>% filter(V2 == "S") %>%
  filter(begin > 39.562 & begin < 39.625) %>% 
  distinct(V5)
# Group 30 has been moved to i = 6

output.23.5$data %>% 
  filter(begin > 39.562 & begin < 39.625) %>% 
  group_by(V5) %>%
  summarise(max.n = max(V9),
            last.n = last(V9),
            d.max.last = max(V9) - last(V9))
  
data.23 %>% filter(V2 == "S") %>%
  filter(begin > 39.562 & begin < 39.625) %>% 
  group_by(V5) %>%
  filter(V5 == "25")

# Group 25 changed the number from 2, to 7, to 5. This explains the discrepancy of 2.

#########################################################################
FinalData.Both %>% filter(ff == 29) %>% filter(v == "J")
FinalData.Both %>% filter(ff == 29) %>% filter(v == "T")
# i = 3 contains 1 more whale (37) in Tomos than in Josh's (36)

# pull out the data, examine the spillover groups
data.29 <- get.data("Data/", 2020, ff = 29)

output.29.3 <- get.shift(2020, data.29, i = 3)
output.29.4 <- get.shift(2020, data.23, i = 4)

ID.in <- output.29.3$data %>% distinct(V5)

ID.all <- data.29 %>% filter(V2 == "S") %>%
  filter(begin > 48.4 & begin < 48.5) %>%
  distinct(V5)

# These were moved to i = 4
ID.spillover <- ID.all$V5[!(ID.all$V5 %in% ID.in$V5)]

output.29.3$data %>% 
  filter(V2 == "S") %>%
  filter(begin > 48.4 & begin < 48.5) %>%
  group_by(V5) %>%
  summarise(max.n = max(V9),
            last.n = last(V9),
            d.max.last = max(V9) - last(V9))

data.29 %>% filter(V2 == "S") %>%
  filter(begin > 48.4 & begin < 48.5) %>% 
  group_by(V5) %>%
  filter(V5 == "26")

# Group 26 changed the number from 3, to 1, to 2. This explains the discrepancy of 1.

###########################################################################
FinalData.Both %>% filter(ff == 34) 
# Josh's missed these two sightings (i = 1 and 2)
# Need to figure out why these were dropped in Josh's analysis.

# Josh's CorrectLength contains two lines but vs is NA and 8!
Josh.out$CorrectLength %>% filter(ff == 34)


data.34 <- get.data("Data/", 2020, ff = 34)

output.34.1 <- get.shift(2020, data.34, i = 1)
output.34.2 <- get.shift(2020, data.34, i = 2)
output.34.3 <- get.shift(2020, data.34, i = 3)

ID.in <- output.34.1$data %>% distinct(V5)

ID.all <- data.34 %>% 
  filter(V2 == "S") %>% 
  filter(begin > 55.3 & begin < 55.38) %>%
  distinct(V5)

# this was moved to i = 2
ID.spillover <- ID.all$V5[!(ID.all$V5 %in% ID.in$V5)]

output.34.1$data %>% 
  filter(V2 == "S") %>%
  filter(begin > 55.3 & begin < 55.38) %>%
  group_by(V5) %>%
  summarise(max.n = max(V9),
            last.n = last(V9),
            d.max.last = max(V9) - last(V9))

data.34 %>% filter(V2 == "S") %>%
  filter(begin > 55.3 & begin < 55.38) %>% 
  group_by(V5) %>%
  filter(V5 == "4")

# Group 26 changed the number from 3, to 1, to 2. This explains the discrepancy of 1.


FinalData.Both %>% filter(ff == 39) 
# Josh's missed i = 5 (n = 8)

FinalData.Both %>% filter(ff == 43) %>% filter(v == "J")
FinalData.Both %>% filter(ff == 43) %>% filter(v == "T")
# i = 1 contains 2 more whales (11) in Tomos than in Josh's (9)


FinalData.Both %>% filter(ff == 44) %>% filter(v == "J")
FinalData.Both %>% filter(ff == 44) %>% filter(v == "T")
# i = 4 contains 1 more whales (12) in Tomos than in Josh's (11)


FinalData.Both %>% filter(ff == 47) %>% filter(v == "J")
FinalData.Both %>% filter(ff == 47) %>% filter(v == "T")
# i = 1 contains 1 more whales (18) in Tomos than in Josh's (17)
# i = 6 contains 1 more whales (27) in Tomos than in Josh's (26)


ggplot(data = FinalData.Both) +
  geom_point(aes(x = begin, y = end, color = v))

delta.begin <- FinalData.Josh$begin - FinalData.Tomo$begin
sum(abs(delta.begin))
diff.begin.Josh <- FinalData.Josh[abs(delta.begin) > 0,]
diff.begin.Tomo <- FinalData.Tomo[abs(delta.begin) > 0,]
row.idx.diff.begin <- row.names(diff.begin)
file.idx.diff.begin <- unique(diff.begin$ff)
file.idx.diff.begin

sum(abs(FinalData.Josh$end - FinalData.Tomo$end))

sum(abs(FinalData.Josh$bf - FinalData.Tomo$bf))
FinalData.Josh[abs(FinalData.Josh$bf - FinalData.Tomo$bf) > 0,]


sum(abs(FinalData.Josh$vs - FinalData.Tomo$vs))


