
# define some functions


# Multiple plot function
# from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



Richards_fcn <- function(d, S1, S2, K, P, min, max){
  K <- abs(K)
  if (S1 > 0) S1 <- -S1
  if (S2 < 0) S2 <- -S2
  
  M1 <- (1 + (2 * exp(K) - 1) * exp((1/S1) * (P - d))) ^ (-1/exp(K))
  M2 <- (1 + (2 * exp(K) - 1) * exp((1/S2) * (P - d))) ^ (-1/exp(K))
  N <- min + (max - min) * (M1 * M2)
  return(N)
}




# A function to get one data file from selected directory
# Inputs are data directory name, year of survey (2021/2022 is 2022), and
# which file to be extracted (sequential number from 1 to length(files)).
get.data <- function(dir, YEAR, ff){
  FILES <- list.files(paste0(dir, "/", YEAR))
  all.lines <- read_lines(file = paste0(dir, "/", YEAR, "/", FILES[ff]))
  
  # look at all event code
  event.code <- str_sub(all.lines, start = 5, end = 5)
  
  # are there any comments?
  COMMENTS <- which(event.code == "C")
  
  # if there were comments, remove them
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
  
  #data <- data[, colSums(is.na(data)) < nrow(data)]
  
  # Files with extensive comments in sightings are problematic because they get split into multiple lines. 
  # We need to pull out lines that contain only numeric V1 (when they are converted into numeric)
  
  data <- data[!is.na(as.numeric(data$V1)),]
  
  Starts <- which(data$V2=="B") #Find all start times
  Ends <- which(data$V2=="E") #Find all end times
  
  # if there is no "E" at the end, Add "E" with time equal to 5 seconds after
  # the last entry.
  if (length(Ends) == 0 | max(Ends) != nrow(data)){
    row.num <- as.numeric(data[nrow(data), 1])
    row.num.char <- ifelse(row.num < 100, 
                           paste0("0", as.character(row.num+1)),
                           as.character(row.num + 1))
    tmp <- hour(hms(data[nrow(data),4])) + 
      minute(hms(data[nrow(data),4]))/60 + 
      (second(hms(data[nrow(data),4])) + 5)/3600 
    
    h <- trunc(tmp)
    m <- trunc((tmp - h) * 60)
    s <- (((tmp - h) * 60) - m) * 60
    data <- rbind(data, c(row.num.char, "E", data[nrow(data), 3],
                          paste(h,m,s, sep = ":"),
                          rep(NA, times = 12)))
  }
  
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
  BeginHr <- hour(hms(data$V4)) + minute(hms(data$V4))/60 + second(hms(data$V4))/3600
  
  data %>% 
    mutate(begin = as.numeric(BeginDay) + BeginHr/24,
           shift = cumsum(V2=="P")) -> data
  return(data)
}

# 
# A function to extract one shift from a data file. Use get.data first and
# use the output of get.data in this function. 
get.shift <- function(YEAR, data, ff, i){
  # Each shift always begins with "P"
  Shifts.begin <- which(data$V2 %in% "P")
  Shifts.begin.df <- data.frame(event = "P",
                                shift = 1:length(Shifts.begin),
                                row = Shifts.begin)
  # But each shift does not always have an explicit end. 
  # The end of data file does not always contain "E" either.
  Shifts.end <- which(data$V2 %in% "E")
  Shifts.end.df <- data.frame(event = "E",
                              shift = NA,
                              row = Shifts.end)
  
  Shifts.df <- arrange(rbind(Shifts.begin.df, Shifts.end.df), row)
  
  max.shifts <- length(Shifts.begin)
  #Only use the first observer for model random effect
  Observer <- data[Shifts.begin[i], 5] 
  # Days since Nov 30th of the previous year
  BeginDay <- mdy(data[Shifts.begin[i], 3]) - mdy(paste0("11/30/", (YEAR - 1)))
  # Decimal hour of shift start time - need to add seconds because sometimes
  # the last sighting and next shift starts within one minute. This happened
  # in a 2020 data file (file 41, 2020-02-04)
  
  # Beginning hr of the shift.
  BeginHr <- (hour(hms(data[Shifts.begin[i], 4])) + 
                (minute(hms(data[Shifts.begin[i], 4]))/60) 
              + (second(hms(data[Shifts.begin[i], 4]))/3600))
  
  # Decimal hour of next shift start time
  if (i < max.shifts){
    event.idx <- which(Shifts.df$shift %in% i)
    next.event <- Shifts.df[event.idx + 1,]
    if (next.event$event == "P"){
      NextBeginHr <- (hour(hms(data[next.event$row, 4])) + 
                        (minute(hms(data[next.event$row, 4]))/60)
                      + (second(hms(data[next.event$row, 4]))/3600))
      EndHr <- NextBeginHr - 0.00001
    } else {  # if the event is "E"
      next.P <- Shifts.df[event.idx+2,]
      NextBeginHr <- (hour(hms(data[next.P$row, 4])) + 
                        (minute(hms(data[next.P$row, 4]))/60) 
                      + (second(hms(data[next.P$row, 4]))/3600))
      
      EndHr <- (hour(hms(data[next.event$row, 4])) + 
                  (minute(hms(data[next.event$row, 4]))/60) + 
                  (second(hms(data[next.event$row, 4]))/3600))
    }

    # Find next end hr to find the next shift to figure out spillovers
    event.idx2 <- which(Shifts.df$shift %in% (i+1))
    next.event2 <- Shifts.df[event.idx2+1,]
    
    if (next.event2$event == "P"){   # if the next event is also "P"
       NextEndHr <- (hour(hms(data[next.event2$row, 4])) + 
                       (minute(hms(data[next.event2$row, 4]))/60) + 
                       (second(hms(data[next.event2$row, 4]))/3600) ) - 0.00001
       # 
    } else {  # if the event is "E"
       NextEndHr <- (hour(hms(data[next.event2$row, 4])) + 
                       (minute(hms(data[next.event2$row, 4]))/60) 
                     + (second(hms(data[next.event2$row, 4]))/3600))
    }
    

  } else {    # for the last shift
    event.idx <- which(Shifts.df$shift %in% max.shifts)
    next.event <- Shifts.df[event.idx + 1,]  # This has to be E
    if (length(next.event) == 0){
      end.row <- nrow(data)
    } else {
      end.row <- next.event$row
    }
    EndHr <-  (hour(hms(data[end.row, 4])) + 
                 (minute(hms(data[end.row, 4]))/60) + 
                 (second(hms(data[end.row, 4]))/3600))
  }
  
  
  # End time is just before next start time (replicating J Durban's calculations)
  # TE: This is incorrect. If there was "E", we should use it. 
  #EndHr <- NextBeginHr - 0.00001 
  # Beginning time as a decimal day - NextBeginHr needs to include seconds for the
  # rare occasions when a sighting happens within a minute of the start of a
  # shift. 
  Begin <- as.numeric(BeginDay) + (BeginHr/24) 
  #End time as a decimal day
  End <- as.numeric(BeginDay) + (EndHr/24)
  
  data.shift <- data %>% filter(begin >= Begin & 
                                  begin <= End)
  
  if (i < max.shifts){
    # when there are multiple Es in one file: Take the first of positive values
    if (length(NextBeginHr) > 1){
      dif.BeginHr <- NextBeginHr - BeginHr
      NextBeginHr <- NextBeginHr[dif.BeginHr>0] %>% first()
    }
    
    # following shift for finding out spillovers.
    data.shift2 <- data %>% filter(begin >= as.numeric(BeginDay) + NextBeginHr/24 & 
                                     begin <= as.numeric(BeginDay) + NextEndHr/24)
    
  } else {
    data.shift2 <- NA
  }

    # pull out the data for this shift
  # if (i < max.shifts){
  #   data.shift <- data[(Shifts.begin[i]):(Shifts.begin[i+1]-1),]    
  # } else {
  #   data.shift <- data %>% filter(begin >= data[Shifts.begin[i], "begin"] & 
  #                                   begin <= End)
  # }

  # This really removes information from "V" entries because there is no
  # V12 when "V" is entered (i.e., NA). I think this is better because
  # there was a case (1/12/2015 9:00 - 10:30) where visibility changed from 4
  # to 5 at the end of a period (30 sec from the end). I think those sightings
  # should be included, rather than excluded. The same is true for VS below. 
  #
  # However... there were shifts (e.g., 2/5/2015 9:00-10:30) where BF changed to
  # 5 in the middle of the shift (1 hr in), so that shift should be removed. But
  # by using just "S" entries, it was not removed. So, I need to fix that. 2022-03-30
  BFs.S <- data.shift %>% 
    filter(V2 == "S") %>% 
    select(V12, begin) %>%
    transmute(BF = as.numeric(V12),
              time = begin)
    #pull() %>% as.numeric()
  
  BFs.V <- data.shift %>% 
    filter(V2 == "V") %>% 
    select(V5, begin) %>%
    transmute(BF = as.numeric(V5),
              time = begin)
  
    #pull()%>% as.numeric()
  
  BFs.dt <- rbind(BFs.S, BFs.V) %>%
    arrange(time) %>%
    mutate(dt = (time - min(time)) * (24*60))  # dt in minutes
  
  # if BF changed to 5 within the last 5 minutes, keep the entire period (max(BF) < 5)
  # but if BF changed to 5 before then, make the max BF = 5.
  BFs.dt %>% 
    filter(BF > 4) -> High.BFs
  
  if (nrow(High.BFs) == 0) {
    BFs <- BFs.dt$BF
  } else {
    if (min(High.BFs$dt) > 85){   # when the first change happened within 5 min of the shift change
      BFs <- BFs.dt$BF[BFs.dt$BF < 5]
    } else {
      BFs <- BFs.dt$BF
    }
  }
  # if (i < max.shifts){
  #   BFs <- as.numeric(data[Shifts.begin[i]:(Shifts.begin[i+1]-1), 12])    
  # } else {
  #   BFs <- data.shift %>% 
  #     select(V12) %>% pull()
  # }
  
  
  if (sum(!is.na(BFs)) == 0){
    BF <- NA
  } else {
    BF <- max(BFs, na.rm=T)
  }
  
  VSs.S <- data.shift %>% 
    filter(V2 == "S") %>% 
    select(V13, begin) %>%
    transmute(VS = as.numeric(V13),
              time = begin)
  #pull() %>% as.numeric()
  
  VSs.V <- data.shift %>% 
    filter(V2 == "V") %>% 
    select(V6, begin) %>%
    transmute(VS = as.numeric(V6),
              time = begin)
  
  #pull()%>% as.numeric()
  
  VSs.dt <- rbind(VSs.S, VSs.V) %>%
    arrange(time) %>%
    mutate(dt = (time - min(time)) * (24*60))  # dt in minutes
  
  # if VS changed to 5 within the last 5 minutes, keep the entire period (max(VS) < 5)
  # but if VS changed to 5 before then, make the max VS = 5.
  VSs.dt %>% 
    filter(VS > 4) -> High.VSs
  
  if (nrow(High.VSs) == 0) {
    VSs <- VSs.dt$VS
  } else {
    if (min(High.VSs$dt) > 85){   # when the first change happened within 5 min of the shift change
      VSs <- VSs.dt$VS[VSs.dt$VS < 5]
    } else {
      VSs <- VSs.dt$VS
    }
  }
  
  if (sum(!is.na(VSs)) == 0){
    VS <- NA
  } else {
    VS <- max(VSs, na.rm=T)
  }

  # if still NA, take the first "V" entry
  if (is.na(BF)) {BF <- data[Shifts.begin[i]+1, 5]}
  if (is.na(VS)) {VS <- data[Shifts.begin[i]+1, 6]}
  
  Spillover <- vector(length = 0)
  # No spillover for the last shift
  if (i < max.shifts){
    # Groups = Observers. Only the first (primary) observer is considered (V5)
    # Group numbers from this watch period
    GroupsThisWatch <- data.shift %>% #[Shifts.begin[i]:(Shifts.begin[i+1]-1),] %>% 
      filter(V2 == "S") %>%
      distinct(V5) %>%
      pull()
    # Group numbers from next watch period
    GroupsNextWatch <- data.shift2 %>% #[Shifts.begin[i+1]:(Shifts.begin[i+2]-1),] %>% 
      filter(V2 == "S") %>%
      distinct(V5) %>%
      pull()
    
    # Which groups from watch i were also observed in watch i+1? They should be excluded from i and counted in i+1
    Spillover <- GroupsThisWatch[GroupsThisWatch %in% GroupsNextWatch] 
    
  }
  
  # Changed V14 != "North" to tolower(V14) != "north" because of entries using
  # uppercase "NORTH" in 2015 (file 11), which was not filtered out correctly in V2. 
  if(length(Spillover > 0)){ #if there are groups that spill over into following watch, 
    # figure out if there were any sightings that need to be considered:
    sub.data <- data.shift %>% #data[(Shifts.begin[i]):(Shifts.begin[i+1]-1),] %>% 
      filter(V2 == "S", !(V5 %in% Spillover), tolower(V14) != "north")
    
    if (nrow(sub.data) > 0){
      N <- sub.data %>%
        group_by(V5) %>% #group by the whale group number
        select(V5, V9) %>%
        #summarize(N = max(as.numeric(V9), na.rm = T)) %>% 
        summarize(N = last(as.numeric(V9))) %>% 
        select(N)  %>% sum()
    } else {
      N <- 0
    }
    
  } else {   # if there were no spillover
    sub.data <- data.shift %>% #data[Shifts.begin[i]:(Shifts.begin[i+1]-1),]  %>%  
      filter(V2 == "S", tolower(V14) != "north")
    
    if (nrow(sub.data) > 0){
      N <- sub.data %>%
        group_by(V5) %>% #group by the whale group number
        select(V5, V9) %>%
        #summarize(N = max(as.numeric(V9), na.rm = T)) %>% 
        summarize(N = last(as.numeric(V9))) %>% 
        select(N)  %>% sum()
    } else {
      N <- 0
    }
    
  }
  
  out.list <- list(out.df = data.frame(begin = as.numeric(Begin),
                                       end = as.numeric(End),
                                       dur = as.numeric(End) - as.numeric(Begin),
                                       bf = as.numeric(BF),
                                       vs = as.numeric(VS),
                                       n = N,
                                       obs = as.character(Observer),
                                       ff = ff,
                                       i = i,
                                       BeginHr = BeginHr,
                                       BeginDay = BeginDay),
                   data = sub.data,
                   data.shift = data.shift,
                   data.next.shift = data.shift2)
  return( out.list )
}

# converts a fractional date into YMD and hms
fractional_Day2YMDhms <- function(x, YEAR){
  n.days <- floor(x)
  dec.hr <- (x - n.days) * 24
  hr <- floor(dec.hr)
  dec.min <- (dec.hr - hr) * 60
  m <- floor(dec.min)
  s <- floor((dec.min - m) * 60)
  

  mdy <- n.days + as.Date(paste0((YEAR-1), "-11-30"))          

  return(list(YMD = mdy,
              hms = paste(ifelse(hr < 10, paste0("0", hr), hr), 
                          ifelse(m < 10, paste0("0", m), m), 
                          ifelse(s < 10, paste0("0", s), s), sep = ":")))  
}

# This function compares Ver1.0 and Ver2.0 data extraction code for using
# raw data files (starting the 2020 season). The raw data files should be 
# analyzed using Formatting GC Data TE.R for Ver1.0 (saves in a .rds file)
# and Extract_Data_All_v2.Rmd for Ver2.0 (saves in a .rds file).
compare.V0.V2.raw <- function(YEAR, obs.list){
  v0.out <- readRDS(paste0("RData/out_", YEAR, "_Joshs.rds"))
  v2.out <- readRDS(paste0("RData/out_", YEAR, "_Tomo_v2.rds"))
  
  FinalData.v2 <- v2.out$FinalData %>% mutate(v = "V2") %>% select(-dur)
  FinalData.v0 <- v0.out$FinalData %>% mutate(v = "V0") 
  
  FinalData.v2 <- v2.out$FinalData %>% 
    mutate(v = "V2") %>% 
    left_join(obs.list, by = "obs") 
  
  # find if there is NA in ID - not in the look up table  
  ID.NA <- filter(FinalData.v2, is.na(ID))
  
  unique.ID.NA <- unique(ID.NA$obs)

  if (length(unique.ID.NA) > 0){
    new.obs <- data.frame(obs = NA, ID = NA)
    for (k in 1:length(unique.ID.NA)){
      FinalData.v2[FinalData.v2$obs == unique.ID.NA[k], "ID"] <- max(obs.list$ID) + k
      new.obs[k,] <- c(unique.ID.NA[k], max(obs.list$ID)+k)
    }
    obs.list <- rbind(obs.list, new.obs)
    
  }
  
  
  # replace column names
  FinalData.v2 %>% select(-obs) %>%
    mutate(obs = ID) %>%
    select(-ID) -> FinalData.v2
  
  # rearrange the columns to match v0
  FinalData.v2 <- FinalData.v2[, names(FinalData.v0)]
  FinalData.Both <- rbind(FinalData.v2, FinalData.v0)
  
  min.begin <- min(floor(FinalData.Both$begin))
  max.begin <- max(ceiling(FinalData.Both$begin))
  
  time.steps <- min.begin:max.begin
  difs <- data.frame(begin = double(),
                     end = double(),
                     min.begin = double(), 
                     max.end = double(), 
                     n.periods = integer(), 
                     max.bf = integer(), 
                     max.vs = integer(), 
                     total.whales = integer(),
                     time.step = integer(),
                     stringsAsFactors = F)
  
  c <- k <- 1
  for (k in 1:(length(time.steps)-1)){
    tmp <- filter(FinalData.Both, begin >= time.steps[k] & begin < time.steps[k+1])
    if (nrow(tmp) > 0){
      tmp %>% filter(v == "V0") -> tmp.1
      tmp %>% filter(v == "V2") -> tmp.2
      
      difs[c,] <- c(min(tmp$begin), 
                    max(tmp$end),
                    min(tmp.1$begin) - min(tmp.2$begin), 
                    max(tmp.1$end) - max(tmp.2$end),
                    nrow(tmp.1) - nrow(tmp.2),
                    max(tmp.1$bf) - max(tmp.1$bf),
                    max(tmp.1$vs) - max(tmp.1$vs),
                    sum(tmp.1$n) - sum(tmp.2$n),
                    time.steps[k])
      c <- c + 1
      
    }
    
  }  
  
  difs %>% filter(n.periods != 0 | total.whales != 0) -> difs.1
  FinalData.Both %>% mutate(time.steps = floor(FinalData.Both$begin)) -> FinalData.Both
  
  v2.out$Data_Out %>% 
    mutate(time.steps = floor(v2.out$Data_Out$begin)) -> Data_Out.v2 
  
  v2.out$CorrectLength %>%
    mutate(time.steps = floor(v2.out$CorrectLength$begin)) -> CorrectLength.v2 
  
  return(out.list <- list(difs = difs,
                          difs.1 = difs.1,
                          FinalData.Both = FinalData.Both,
                          Data_Out.v2 = Data_Out.v2,
                          CorrectLength.v2 = CorrectLength.v2,
                          v0.out = v0.out,
                          v2.out = v2.out,
                          obs.list = obs.list))
}

# This function compares outputs from data extraction codes. It uses BUGS input
# data (for data before the 2020 season because the old version (Ver1.0) does not
# work for old files) and Ver2.0. The raw data files should be 
# analyzed using Extract_Data_All_v2.Rmd for Ver2.0 (saves in a .rds file).
compare.V0.V2.BUGSinput <- function(YEAR, idx.yr, periods, obs.list){
  
  #Watch start times, as fraction of a day - stored in a different file
  begin <- as.matrix(read.table("Data/begin.txt", 
                                header=T, 
                                nrows = max(periods)))
  
  #watch end times
  end <- as.matrix(read.table("Data/end.txt", 
                              header=T,
                              nrows = max(periods)))
  
  
  # this file contains all input data for WinBUGS.
  V0.out <- readRDS("RData/2006-2019_GC_Formatted_Data.RDS")
  
  # Pull out the information for 2015
  periods.2015 <- V0.out$periods[idx.yr]
  n.2015 <- V0.out$n[1:periods.2015,,idx.yr]
  n.com.2015 <- V0.out$n.com[1:periods.2015,,idx.yr]
  n.sp.2015 <- V0.out$n.sp[1:periods.2015,,idx.yr]
  obs.2015 <- V0.out$obs[1:periods.2015,,idx.yr]
  
  vs.2015 <- V0.out$vs[1:periods.2015,idx.yr]
  bf.2015 <- V0.out$bf[1:periods.2015,idx.yr]
  day.2015 <- V0.out$day[1:periods.2015,idx.yr]
  
  FinalData.V0 <- data.frame(begin = begin[1:periods[idx.yr], idx.yr],
                             end = end[1:periods[idx.yr], idx.yr],
                             bf = bf.2015,
                             vs = vs.2015,
                             n = n.2015[,1],
                             obs = obs.2015[,1],
                             BeginDay = day.2015,
                             v = "V0")
  
  # This contains the results from my version
  v2.out <- readRDS(paste0("RData/out_", YEAR, "_Tomo_v2.rds"))
  FinalData.v2 <- v2.out$FinalData %>% 
    mutate(v = "V2") %>% 
    left_join(obs.list, by = "obs") %>%
    select(-c(dur, ff, i, BeginHr)) 
  
  # find if there is NA in ID - not in the look up table  
  ID.NA <- filter(FinalData.v2, is.na(ID))
  
  unique.ID.NA <- unique(ID.NA$obs)
  
  if (length(unique.ID.NA) > 0){
    new.obs <- data.frame(obs = NA, ID = NA)
    
    for (k in 1:length(unique.ID.NA)){
      FinalData.v2[FinalData.v2$obs == unique.ID.NA[k], "ID"] <- max(obs.list$ID) + k
      new.obs[k,] <- c(unique.ID.NA[k], as.numeric(max(obs.list$ID)+k))
    }
    obs.list <- rbind(obs.list, new.obs)
    
  }
  
  # replace column names
  FinalData.v2 %>% 
    select(-obs) %>%
    mutate(obs = ID) %>%
    select(-ID) -> FinalData.v2
  
  # rearrange the columns to match V0
  FinalData.v2 <- FinalData.v2[, names(FinalData.V0)]
  FinalData.Both <- rbind(FinalData.v2, FinalData.V0)
  
  v2.out$Data_Out %>% 
    mutate(time.steps = floor(v2.out$Data_Out$begin)) -> Data_Out.v2 
  
  v2.out$CorrectLength %>%
    mutate(time.steps = floor(v2.out$CorrectLength$begin)) -> CorrectLength.v2 
  
  min.begin <- min(floor(FinalData.Both$begin))
  max.begin <- max(ceiling(FinalData.Both$begin))
  
  time.steps <- min.begin:max.begin
  difs <- data.frame(begin = double(),
                     end = double(),
                     min.begin = double(), 
                     max.end = double(), 
                     n.periods = integer(), 
                     max.bf = integer(), 
                     max.vs = integer(), 
                     total.whales = integer(),
                     time.step = integer(),
                     stringsAsFactors = F)
  
  c <- k <- 1
  for (k in 1:(length(time.steps)-1)){
    tmp <- filter(FinalData.Both, begin >= time.steps[k] & begin < time.steps[k+1])
    if (nrow(tmp) > 0){
      tmp %>% filter(v == "V0") -> tmp.1
      tmp %>% filter(v == "V2") -> tmp.2
      
      difs[c,] <- c(min(tmp$begin), 
                    max(tmp$end),
                    min(tmp.1$begin) - min(tmp.2$begin), 
                    max(tmp.1$end) - max(tmp.2$end),
                    nrow(tmp.1) - nrow(tmp.2),
                    max(tmp.1$bf) - max(tmp.1$bf),
                    max(tmp.1$vs) - max(tmp.1$vs),
                    sum(tmp.1$n) - sum(tmp.2$n),
                    time.steps[k])
      c <- c + 1
      
    }
    
  }
  
  
  difs %>% filter(n.periods != 0 | total.whales != 0) -> difs.1
  FinalData.Both %>% mutate(time.steps = floor(FinalData.Both$begin)) -> FinalData.Both
  
  return(out.list <- list(difs = difs,
                          difs.1 = difs.1,
                          FinalData.Both = FinalData.Both,
                          FinalData.V0 = FinalData.V0,
                          Data_Out.v2 = Data_Out.v2,
                          CorrectLength.v2 = CorrectLength.v2,
                          v0.out = V0.out,
                          v2.out = v2.out,
                          obs.list = obs.list))
  
}

# This function compares whale counts, Beaufort, and visibility for
# all shifts that were recorded by Ver1.0 and Ver2.0 and creates
# a table that is sorted by the beginning time of each shift. 
# The table is returned. 2022-04-01
n.comparison <- function(FinalData.Both, difs.1, idx, YEAR){
  FinalData.Both %>% 
    filter(time.steps == difs.1[idx, "time.step"]) %>%
    mutate(begin.time = fractional_Day2YMDhms(begin, YEAR)$hms,
           end.time = fractional_Day2YMDhms(end, YEAR)$hms) %>% 
    select(begin.time, end.time, n, bf, vs, v) -> tmp
  
  tmp %>%
    filter(v == "V0") -> tmp.0
  tmp %>%
    filter(v == "V2") -> tmp.2
  
  tmp.0 %>% 
    full_join(tmp.2, by = "begin.time") %>%
    arrange(begin.time) %>%
    transmute(begin.time = begin.time,
              end.time.Ver1.0 = end.time.x,
              end.time.Ver2.0 = end.time.y,
              n.Ver1.0 = n.x,
              n.Ver2.0 = n.y,
              vs.Ver1.0 = vs.x,
              vs.Ver2.0 = vs.y,
              Bf.Ver1.0 = bf.x,
              Bf.Ver2.0 = bf.y) -> tmp.0.2
  
  total <- c(NA, NA, "Total", 
           sum(tmp.0.2$n.Ver1.0, na.rm = T), 
           sum(tmp.0.2$n.Ver2.0, na.rm = T), 
           NA, NA, NA, NA)
  
  return(rbind(tmp.0.2, total))  
}


