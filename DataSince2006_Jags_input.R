#data2Jags_input.R

# Creates Jags inputs from data files

data2Jags_input <- function(min.dur = 85, 
                            seasons = c("2006/2007", "2007/2008", "2009/2010", "2010/2011",
                                        "2014/2015", "2015/2016", "2019/2020", "2021/2022",
                                        "2022/2023", "2023/2024"), 
                            WinBUGS.out.file = "RData/WinBUGS_10yr_v2_min85.rds",
                            years = c("2010", "2011", "2015", 
                                      "2016", "2020", "2022", 
                                      "2023", "2024"),
                            n.stations = c(1, 1, 2, 2, rep(1, times = 6)),
                            data.dir = "RData/V2.1_Nov2024"){
  min.dur <- min.dur
  seasons <- seasons
  
  # v2 refers to v2 data extraction. 
  WinBUGS.out <- readRDS(WinBUGS.out.file)
  #data.WinBUGS <- data.v2$BUGS.data
  
  # New as of 2024-11-14
  # Use data2WinBUGS_input_fcn.R
  source("data2WinBUGS_input_fcn.R")
  WinBUGS.inputs <- data2WinBUGS_input(data.dir = data.dir,
                                       years = years,
                                       min.duration = min.dur)
  
  data.WinBUGS <- WinBUGS.inputs$data
  # watch lengths are assumed equal between primary and secondary stations in
  # WinBUGS code. But not in Jags. So, I duplicate the secondary watch effort
  
  bf <- vs <- array(dim = c(max(data.WinBUGS$periods),
                            2, 
                            length(data.WinBUGS$periods)))
  
  watch.prop <- day <- array(dim = c(dim(data.WinBUGS$Watch.Length)[1],
                                     2, 
                                     dim(data.WinBUGS$Watch.Length)[2]))
  
  bf[,1,] <- data.WinBUGS$bf
  bf[,2,] <- data.WinBUGS$bf
  
  vs[,1,] <- data.WinBUGS$vs
  vs[,2,] <- data.WinBUGS$vs
  
  day[,1,] <- data.WinBUGS$day
  day[,2,] <- data.WinBUGS$day
  
  watch.prop[,1,] <- data.WinBUGS$Watch.Length
  watch.prop[,2,] <- data.WinBUGS$Watch.Length
  
  jags.data <- list(  n = data.WinBUGS$n,
                      n.station = n.stations,
                      n.year = length(seasons),
                      n.obs = data.WinBUGS$n.obs,
                      #Daily.N = daily.N,
                      periods = cbind(data.WinBUGS$periods,
                                      data.WinBUGS$periods),
                      n.days = 90,
                      #first.day = unlist(as.vector(first.day)),
                      obs = data.WinBUGS$obs,
                      vs = vs,
                      bf = bf,
                      watch.prop = watch.prop,
                      #watch.prop = (data.WinBUGS$Watch.Length*24*60)/540,
                      day = day)
  
  out.list <- list(jags.data = jags.data,
                   min.dur = min.dur, 
                   seasons = seasons, 
                   WinBUGS.out.file = WinBUGS.out.file,
                   years = years,
                   data.dir = data.dir)
  return(out.list)
}

