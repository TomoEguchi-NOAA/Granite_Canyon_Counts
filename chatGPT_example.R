dat <- read.csv("data/nmixture_simulated_data.csv")


# Build spline basis on day of year
K <- 4
Bday <- ns(day_vec, df = K, intercept = TRUE)

jags_data <- list(
  yobs = yobs_vec,
  x.p = x.p_vec,
  year = as.integer(as.factor(year_vec)),  # convert to 1:nyears
  day = day_vec,
  nyears = length(unique(year_vec)),
  max_day = max(day_vec),
  Bday_indexed = array(NA, dim = c(nyears, max_day, K)),
  nobs = length(yobs_vec),
  K = K
)

# Fill Bday_indexed with spline basis by year and day
for (i in 1:jags_data$nobs) {
  y <- jags_data$year[i]
  d <- jags_data$day[i]
  jags_data$Bday_indexed[y, d, ] <- Bday[i, ]
}

