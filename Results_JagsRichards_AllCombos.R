# Results_JagsRichards_AllCombo.R
# Show results of Run_JagsRichards_AllCombo.R
# 
# 

rm(list=ls())
# Run the Run_ script to get all input information
source("Run_JagsRichards_AllCombo.R")

library(loo)
library(tidyverse)

# out.file.name contains all out put file names.
# Comparisons need to be made among models (v1, v3, v4, and v5), among minimum
# watch durations (10, 30, 85), and among different data creation processes 
# (NoBUGS, AllYears, Since2006)
# 
# Metrics to compare are convergence (maximum Rhat) and model fit (LOOIC extreme
# values). Then compare estimated abundance and their precisions.

#file.name <- out.file.name[[1]]
get.results <- function(file.name){
  
  out <- readRDS(paste0("RData/", file.name))
  max.Rhats <- lapply(out$jm$Rhat, max, na.rm = T)
  
  # To compute LOOIC, need to turn zeros into NAs when there were no second station:
  jags.data <- out$jags.input$jags.data
  data.array <- jags.data$n
  data.array[,2,which(jags.data$n.station == 1)] <- NA
  
  LOOIC.n <- compute.LOOIC(loglik.array = out$jm$sims.list$log.lkhd,
                           data.array = data.array,
                           MCMC.params = out$MCMC.params)
  
  bad.Pareto <- length(which(LOOIC.n$loo.out$diagnostics$pareto_k>0.7))/length(LOOIC.n$loo.out$diagnostics$pareto_k)
  max.Pareto <- max(LOOIC.n$loo.out$diagnostics$pareto_k)
  
  # model version
  filename.parts <- unlist(strsplit(out$jags.model, "_"))
  model.v <- filename.parts[length(filename.parts)] %>% strsplit(".txt") %>% unlist()
  
  # watch duration - I forgot to extract minimum watch duration for Laake data
  # analysis... It has been fixed but rather than re-running the code, I extract
  # this information from the file name
  watch.dur <- out$jags.input$min.dur
  if (is.null(watch.dur)){
    file.name.parts <- unlist(strsplit(file.name, "min"))
    watch.dur <- unlist(strsplit(file.name.parts[2], "_"))[1] %>% as.numeric()
  }
  
  # Find out estimates - these are different among which datasets were used
  if (stringr::str_detect(file.name, paste0("min", watch.dur, "_NoBUGS_"))){
    all.start.years <- c(out$jags.input$jags.input.Laake$all.start.year,
                         out$jags.input$jags.input.new$start.years)
    data.set <- "_NoBUGS_"
    
  } else if (stringr::str_detect(file.name, "AllYears")){
    all.start.years <- out$jags.input$start.years  
    data.set <- "AllYears"
  } else if (stringr::str_detect(file.name, "Since2006")){
    all.start.years <- out$jags.input$start.years
    data.set <- "Since2006"
  } else if (stringr::str_detect(file.name, "LaakeData")){
    all.start.years <- out$jags.input$all.start.year
    data.set <- "LaakeData"
  } else if (stringr::str_detect(file.name, "_Since2010_NoBUGS")){
    all.start.years <- out$jags.input$start.years
    data.set <- "Since2010_NoBUGS"
  }
  
  model.fit <- data.frame(model = model.v,
                          data.set = data.set,
                          watch.dur = watch.dur,
                          bad.Pareto = bad.Pareto,
                          max.Pareto = max.Pareto,
                          max.Rhats = max(unlist(max.Rhats), na.rm = T))
  
  Nhats.df <- data.frame(start.year = all.start.years,
                         Mean = out$jm$mean$Corrected.Est,
                         SE = out$jm$sd$Corrected.Est,
                         LCL = out$jm$q2.5$Corrected.Est,
                         UCL = out$jm$q97.5$Corrected.Est,
                         model = model.v,
                         min.watch = watch.dur,
                         data.set = data.set)
  
  return(out.list <- list(model.fit = model.fit,
                          Nhats = Nhats.df))
  
}

all.results <- lapply(out.file.name, FUN = get.results)

all.model.fit <- lapply(all.results, 
                        function(x) x$model.fit)

all.model.fit.df <-  do.call(rbind, all.model.fit)

library(ggplot2)
p.model.fit <- ggplot(all.model.fit.df) +
  geom_point(aes(x = watch.dur, 
                 y = bad.Pareto, 
                 color = model,
                 shape = data.set))

# ggsave("figures/model_fit.png", p.model.fit, device = "png", dpi = 600)

p.max.pareto <- ggplot(all.model.fit.df) +
  geom_point(aes(x = watch.dur, 
                 y = max.Pareto, 
                 color = model,
                 shape = data.set))

# ggsave("figures/max_pareto.png", p.max.pareto, device = "png", dpi = 600)

p.Rhat <- ggplot(all.model.fit.df) +
  geom_point(aes(x = watch.dur, 
                 y = max.Rhats, 
                 color = model,
                 shape = data.set))

# ggsave("figures/max_Rhat.png", p.Rhat, device = "png", dpi = 600)

jags.Nhats <- lapply(all.results, function(x) x$Nhats)
jags.Nhats.df <- do.call(rbind, jags.Nhats) %>%
  mutate(Method = "Eguchi")

# WinBUGS results:
BUGS.all.file.names <- list.files(path = "RData/", pattern = "WinBUGS_2007to2024_v2_min") 
BUGS.file.names <- BUGS.all.file.names[str_detect(BUGS.all.file.names, "_85000_")]

BUGS.file.name <- BUGS.file.names[1]
get.results.BUGS <- function(BUGS.file.name){
  out <- readRDS(paste0("RData/", BUGS.file.name))
  out$BUGS.out$summary %>%
    as.data.frame() %>%
    rownames_to_column(var = "parameter") %>%
    filter(grepl("Corrected", parameter)) -> summary.Nhat

  file.name.parts <- strsplit(BUGS.file.name, "_") %>% unlist()
  watch.dur <- strsplit(file.name.parts[4], "min") %>% unlist() %>% as.numeric()
  
  WinBUGS.Nhats.df <- data.frame(start.year = out$BUGS.input$all.years - 1,
                                 Mean = summary.Nhat$mean,
                                 SE = summary.Nhat$sd,
                                 LCL = summary.Nhat$`2.5%`,
                                 UCL = summary.Nhat$`97.5%`,
                                 model = "BUGS",
                                 min.watch = watch.dur[2],
                                 Method = "Durban",
                                 data.set = "2007to2024")
 return(WinBUGS.Nhats.df) 
  
}

BUGS.Nhats <- lapply(BUGS.file.names, FUN = get.results.BUGS)
BUGS.Nhats.df <- do.call(rbind, BUGS.Nhats)

published.Nhats <- read.csv(file = "Data/all_estimates_2024.csv") %>%
  select(Year, Nhat, LCL, UCL, Method) %>%
  mutate(start.year = Year - 1,
         min.watch = 85) %>%
  mutate(model = ifelse(Method == "Laake", "GAM", "BUGS_pub"),
         SE = NA,
         data.set = "1967to2006") %>%
  select(start.year, Nhat, SE, LCL, UCL, model, min.watch, Method, data.set) %>%
  rename(Mean = Nhat)

all.Nhats.df <- rbind(jags.Nhats.df, 
                      BUGS.Nhats.df, 
                      published.Nhats) %>%
  mutate(CV = SE/Mean,
         model.f = as.factor(model),
         Method.f = as.factor(Method))

p.CV <- ggplot(all.Nhats.df) +
  geom_point(aes(x = min.watch,
                 y = CV,
                 color = model.f))

p.Nhats <- ggplot(all.Nhats.df) +
  geom_point(aes(x = start.year, 
                 y = Mean, 
                 color = model.f, 
                 size = min.watch),
             alpha = 0.5) +
  geom_errorbar(aes(x = start.year, 
                    ymin = LCL, 
                    ymax = UCL, 
                    color = model.f)) +
  facet_wrap(~ data.set)

# v4 and v5 are not doing so well... 
all.Nhats.df %>% 
  filter(model != "v4" & model != "v5") -> all.Nhats.df.2

p.2.Nhats <- ggplot(all.Nhats.df.2) +
  geom_point(aes(x = start.year, 
                 y = Mean, 
                 color = model.f, 
                 size = min.watch),
             alpha = 0.5) +
  geom_errorbar(aes(x = start.year, 
                    ymin = LCL, 
                    ymax = UCL, 
                    color = model.f))

all.Nhats.df.2 %>%
  filter(model != "v1") -> all.Nhats.df.3

p.3.Nhats <- ggplot(all.Nhats.df.3 ) +
  geom_point(aes(x = start.year, 
                 y = Mean, 
                 color = model.f,
                 size = min.watch),
             alpha = 0.5) +
  geom_errorbar(aes(x = start.year, 
                    ymin = LCL, 
                    ymax = UCL, 
                    color = model.f))

all.Nhats.df.3 %>% 
  filter(!(Method == "Eguchi" & min.watch > 20)) -> all.Nhats.df.4

p.4.Nhats <- ggplot(all.Nhats.df.4 ) +
  geom_point(aes(x = start.year, 
                 y = Mean, 
                 color = model.f,
                 size = min.watch),
             alpha = 0.5) +
  geom_errorbar(aes(x = start.year, 
                    ymin = LCL, 
                    ymax = UCL, 
                    color = model.f))

