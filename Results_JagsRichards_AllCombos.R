# Results_JagsRichards_AllCombo.R
# Show results of Run_JagsRichards_AllCombo.R
# 
# 

rm(list=ls())
# Run the Run_ script to get all input information
source("Run_JagsRichards_AllCombo.R")

library(loo)
library(tidyverse)

options(mc.cores = 5)
# out.file.name contains all out put file names.
# Comparisons need to be made among models (v1, v3, v4, and v5), among minimum
# watch durations (10, 30, 85), and among different data creation processes 
# (NoBUGS, AllYears, Since2006)
# 
# Metrics to compare are convergence (maximum Rhat) and model fit (LOOIC extreme
# values). Then compare estimated abundance and their precisions.

#file.name <- out.file.name[[1]]

all.results <- lapply(out.file.name, FUN = get.results.jags)

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

#BUGS.file.name <- BUGS.file.names[1]

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
                    color = model.f))+
  facet_wrap(~ data.set)

# v3 seems to be closer to Laake's and Durban's estimates than v1
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
                    color = model.f)) +
  facet_wrap(~ data.set)

# Using just 10 min watch duration for v3
all.Nhats.df.3 %>% 
  filter(!(Method == "Eguchi" & min.watch > 40)) -> all.Nhats.df.4

p.4.Nhats <- ggplot(all.Nhats.df.4 ) +
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

# For v3, watch duration doesn't seem to make any difference.
ggplot(jags.Nhats.df %>% filter(model == "v3")) +
  geom_point(aes(x = start.year, y = Mean,
                 size = min.watch),
             alpha = 0.5) +
  geom_errorbar(aes(x = start.year, 
                    ymin = LCL, 
                    ymax = UCL))+
  facet_wrap(~ data.set)

# NoBUGS and LaakeData result in about the same Nhats from 1967 to 2005. 
# AllYears results in a different pattern. Why is that?
# Look at AllYears.R for comparing data that went into the analysis. 


p.2.CV <- ggplot(all.Nhats.df.4) +
  geom_point(aes(x = min.watch, y = CV,
                 color = model.f))


# START HERE 2025-02-10
# Compare estiamates among Laake's, noBUGS and WinBUGS. Leave out "AllData" option


