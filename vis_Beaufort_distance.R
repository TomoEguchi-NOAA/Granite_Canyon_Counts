# A script to look at effects of visibility and Beaufort sea state on sighting distances

library(tidyverse)
library(readr)

Laake.data <- read.csv("Data/Laake_PrimarySightings.csv") %>%
  mutate(data = "Laake",
         date = as.Date(Date, format = "%Y-%m-%d")) %>%
  select(date, beaufort, vis, distance, podsize, data)

new.years <- c(2010, 2011, 2015, 2016, 2020, 2022, 2023, 2024)
new.data <- list(length(new.years))
for (y in 1:length(new.years)){
  new.data[[y]] <- read.csv(paste0("Data/all_sightings_", new.years[y], "_Tomo_v2.csv"))
}

new.data.df <- do.call("rbind", new.data) %>%
  mutate(date = as.Date(Date, format = "%m/%d/%Y"),
         data = "New",
         podsize = nwhales,
         distance = Distance) %>%
  select(date, beaufort, vis, distance, podsize, data)

all.data <- rbind(Laake.data, new.data.df) %>%
  na.omit() %>%
  filter(beaufort < 5.1) %>%
  mutate(BF.factor = as.factor(beaufort),
         VS.factor = as.factor(vis))

all.data %>%
  group_by(BF.factor, VS.factor) %>%
  summarise(min.distance = min(distance),
            max.distance = max(distance),
            mean.distance = mean(distance),
            n = n()) -> summary.distance

ggplot(summary.distance) +
  geom_point(aes(x = BF.factor, 
                 y = VS.factor, 
                 #                 color = n,
                 size = mean.distance))

all.data %>%
  group_by(BF.factor, VS.factor) %>%
  summarise(min.pod = min(podsize),
            max.pod = max(podsize),
            mean.pod = mean(podsize),
            n = n()) -> summary.podsize

ggplot(summary.podsize) +
  geom_point(aes(x = BF.factor, 
                 y = VS.factor, 
                 #                 color = n,
                 size = max.pod))

ggplot(all.data) +
  geom_point(aes(x = podsize, y = VS.factor,
                 color = BF.factor,
                 size = distance))

ggplot(all.data) +
  geom_point(aes(x = podsize, y = BF.factor,
                 color = VS.factor,
                 size = distance))

ggplot(all.data) +
  geom_point(aes(x = distance, y = VS.factor,
                 color = BF.factor,
                 size = podsize))

all.data %>% filter(VS.factor == "1") -> data.vis.1
all.data %>% filter(VS.factor == "2") -> data.vis.2
all.data %>% filter(VS.factor == "3") -> data.vis.3
all.data %>% filter(VS.factor == "4") -> data.vis.4

ggplot(all.data) +
  geom_point(aes(x = distance, y = BF.factor,
                 color = VS.factor,
                 size = podsize))

all.data %>% filter(BF.factor == "0") -> data.Bft.0
all.data %>% filter(BF.factor == "1") -> data.Bft.1
all.data %>% filter(BF.factor == "2") -> data.Bft.2
all.data %>% filter(BF.factor == "3") -> data.Bft.3
all.data %>% filter(BF.factor == "4") -> data.Bft.4

ggplot(all.data) + 
  geom_point(aes(x = distance, y = podsize))
