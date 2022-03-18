


# define some functions
# A function to get one data file from selected directory
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
  
  # if there is no "E" at the end
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
  BeginHr <- hour(hms(data$V4)) + minute(hms(data$V4))/60
  
  data %>% mutate(begin = as.numeric(BeginDay) + BeginHr/24) -> data
  return(data)
}

# 
# A function to extract one shift from a data file.
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
  # the last sighting and next shift starts within one minute. 
  
  START HERE NEXT 2022-03-17
  
  BeginHr <- (hour(hms(data[Shifts.begin[i], 4])) + 
                (minute(hms(data[Shifts.begin[i], 4]))/60)) 
  #  + (second(hms(data[Shifts.begin[i], 4]))/3600)
  
  # Decimal hour of next shift start time
  if (i < max.shifts){
    event.idx <- which(Shifts.df$shift %in% i)
    next.event <- Shifts.df[event.idx + 1,]
    if (next.event$event == "P"){
      NextBeginHr <- (hour(hms(data[next.event$row, 4])) + 
                        (minute(hms(data[next.event$row, 4]))/60))
      #  + (second(hms(data[next.event$row, 4]))/3600)
      EndHr <- NextBeginHr - 0.00001
    } else {  # if the event is "E"
      next.P <- Shifts.df[event.idx+2,]
      NextBeginHr <- (hour(hms(data[next.P$row, 4])) + 
                        (minute(hms(data[next.P$row, 4]))/60) )
      # + (second(hms(data[next.P$row, 4]))/3600)
      EndHr <- (hour(hms(data[next.event$row, 4])) + 
                  (minute(hms(data[next.event$row, 4]))/60))
    }

    # Find next end hr to find the next shift to figure out spillovers
    event.idx2 <- which(Shifts.df$shift %in% (i+1))
    next.event2 <- Shifts.df[event.idx2+1,]
    
    if (next.event2$event == "P"){   # if this is the end of the file
       NextEndHr <- (hour(hms(data[next.event2$row, 4])) + 
                       (minute(hms(data[next.event2$row, 4]))/60) ) - 0.00001
       # + (second(hms(data[next.event2$row, 4]))/3600)
    } else {  # if the event is "E"
       NextEndHr <- (hour(hms(data[next.event2$row, 4])) + 
                       (minute(hms(data[next.event2$row, 4]))/60) )
       # + (second(hms(data[next.event2$row, 4]))/3600)
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
                 (minute(hms(data[end.row, 4]))/60))
  }
  
  
  # End time is just before next start time (replicating J Durban's calculations)
  # TE: This is incorrect. If there was "E", we should use it. 
  #EndHr <- NextBeginHr - 0.00001 
  # Beginning time as a decimal day
  Begin <- as.numeric(BeginDay) + (BeginHr/24) 
  #End time as a decimal day
  End <- as.numeric(BeginDay) + (EndHr/24)
  
  data.shift <- data %>% filter(begin >= data[Shifts.begin[i], "begin"] & 
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

  
  BFs <- data.shift %>% 
    filter(V2 == "V" | V2 == "S") %>% 
    select(V12) %>% 
    pull()%>% as.numeric()
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
  
  #Visibility (maximum from watch period)
  VSs <- data.shift %>% 
    filter(V2 == "V" | V2 == "S") %>% 
    select(V13) %>%
    pull() %>% as.numeric()
  # if (i < max.shifts){
  #   VSs <- as.numeric(data[Shifts.begin[i]:(Shifts.begin[i+1]-1), 13])
  # } else {
  #   VSs <- data %>% filter(begin >= data[Shifts.begin[i], "begin"] & 
  #                            begin <= as.numeric(BeginDay + NextBeginHr/24)) %>%
  #     select(V13) %>% pull()
  # }
  
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
  
  if(length(Spillover > 0)){ #if there are groups that spill over into following watch, 
    # figure out if there were any sightings that need to be considered:
    sub.data <- data.shift %>% #data[(Shifts.begin[i]):(Shifts.begin[i+1]-1),] %>% 
      filter(V2 == "S", !(V5 %in% Spillover), V14 != "North")
    
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
      filter(V2 == "S", V14 != "North")
    
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
