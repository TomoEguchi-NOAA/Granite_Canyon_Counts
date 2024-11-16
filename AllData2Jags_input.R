#AllData2Jags_input.R

AllData2JagsInput <- function(min.dur){
  source("DataSince2006_Jags_input.R")
  source("LaakeData2Jags.R")
  
  load("Data/PrimaryEffort.rda")
  Laake.jags.data <- LaakeData2JagsInput()
  Laake.start.year <- unique(PrimaryEffort$Start.year)
  
  Jags.input.2006<- data2Jags_input(min.dur = min.dur)
  .data <- Jags.input.2006$jags.data
  .start.year <- lapply(Jags.input.2006$seasons, 
                        FUN = str_split, pattern = "/") %>% 
    lapply(FUN = function(x) {tmp <- unlist(x); tmp[1]}) %>%
    unlist()
  
  all.start.year <- c(Laake.start.year, .start.year)
  all.start.year <- all.start.year[!duplicated(all.start.year)]
  
  # Both datasets contain 2006. Take it out from the recent one.
  bf <- vs <- all.n <- all.obs <- array(dim = c(max(dim(Laake.jags.data$n)[1], dim(.data$n)[1]), 
                                                2, (dim(Laake.jags.data$n)[3] + dim(.data$n)[3]-1)))
  
  day <- watch.prop <- array(dim = c(max(dim(Laake.jags.data$n)[1], dim(.data$n)[1]), 
                                     2, (dim(Laake.jags.data$n)[3] + dim(.data$n)[3]-1)))
  
  c2 <- 2
  c1 <- c <- 1
  for (k in 1:length(all.start.year)){
    if (all.start.year[k] < 2007){
      all.n[1:dim(Laake.jags.data$n)[1], 1, c] <- Laake.jags.data$n[, 1, c1]
      all.n[1:dim(Laake.jags.data$n)[1], 2, c] <- Laake.jags.data$n[, 2, c1]
      
      all.obs[1:dim(Laake.jags.data$obs)[1], 1, c] <- Laake.jags.data$obs[, 1, c1]
      all.obs[1:dim(Laake.jags.data$obs)[1], 2, c] <- Laake.jags.data$obs[, 2, c1]
      
      watch.prop[1:dim(Laake.jags.data$watch.prop)[1], 1, c] <- Laake.jags.data$watch.prop[, 1, c1]
      watch.prop[1:dim(Laake.jags.data$watch.prop)[1], 2, c] <- Laake.jags.data$watch.prop[, 2, c1]
      
      day[1:dim(Laake.jags.data$day)[1], 1, c] <- Laake.jags.data$day[, 1, c1]
      day[1:dim(Laake.jags.data$day)[1], 2, c] <- Laake.jags.data$day[, 2, c1]
      
      bf[1:dim(Laake.jags.data$bf)[1], 1, c] <- Laake.jags.data$bf[, 1, c1]
      bf[1:dim(Laake.jags.data$bf)[1], 2, c] <- Laake.jags.data$bf[, 2, c1]
      
      vs[1:dim(Laake.jags.data$vs)[1], 1, c] <- Laake.jags.data$vs[, 1, c1]
      vs[1:dim(Laake.jags.data$vs)[1], 2, c] <- Laake.jags.data$vs[, 2, c1]
      
      c <- c + 1
      c1 <- c1 + 1
    } else {
      all.n[1:1:dim(.data$n)[1], 1, c] <- .data$n[, 1, c2]
      all.n[1:1:dim(.data$n)[1], 2, c] <- .data$n[, 2, c2]
      all.obs[1:1:dim(.data$obs)[1], 1, c] <- .data$obs[, 1, c2]
      all.obs[1:1:dim(.data$obs)[1], 2, c] <- .data$obs[, 2, c2]
      watch.prop[1:1:dim(.data$watch.prop)[1], 1, c] <- .data$watch.prop[, 1, c2]
      watch.prop[1:1:dim(.data$watch.prop)[1], 2, c] <- .data$watch.prop[, 2, c2]
      day[1:1:dim(.data$day)[1], 1, c] <- .data$day[, 1, c2]
      day[1:1:dim(.data$day)[1], 2, c] <- .data$day[, 2, c2]
      bf[1:1:dim(.data$bf)[1], 1, c] <- .data$bf[, 1, c2]
      bf[1:1:dim(.data$bf)[1], 2, c] <- .data$bf[, 2, c2]
      vs[1:1:dim(.data$vs)[1], 1, c] <- .data$vs[, 1, c2]
      vs[1:1:dim(.data$vs)[1], 2, c] <- .data$vs[, 2, c2]
      c <- c + 1
      c2 <- c2 + 1
    }
    
  }
  
  n.station <- c(Laake.jags.data$n.station, .data$n.station[2:length(.data$n.station)])
  n.year <- length(all.start.year)
  n.obs <- length(unique(as.vector(all.obs))) 
  
  periods <- rbind(Laake.jags.data$periods, .data$periods[2:length(.data$n.station),])
  
  jags.data <- list(n = all.n,
                    n.station = n.station,
                    n.year = n.year,
                    n.obs = n.obs,
                    periods = periods,
                    obs = all.obs,
                    vs = vs,
                    bf = bf,
                    watch.prop = watch.prop,
                    day = day,
                    n.days = 94) 
  return(jags.data)
}

