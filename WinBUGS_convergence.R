# Calculates WinBUGS convergence statistics
# 
# Bring in the results from a previous run:

rm(list = ls())
YEAR <- 2026
min.dur <- 60#10 #

# It was run in the GrayWhaleAbundance project. 
WinBUGS.Run.Date <- "2026-05-06"  

WinBugs.out <- readRDS(file = paste0("RData/WinBUGS_2007to", YEAR, "_v2_min", 
                                     min.dur,
                                     "_100000_",
                                     WinBUGS.Run.Date, ".rds"))

grp_est  <- "^Corrected|^Raw"
grp_cont <- "^beta|^mean\\.beta|^mean\\.prob|^sigma\\.Obs|^OBS\\.RF|^BF\\.Fixed|^VS\\.Fixed|^b\\.sp|^sd\\.b\\.sp|^corr\\.factor"
grp_disc <- "^z\\[|^N\\[|^N\\.com|^N\\.sp|Switch"

grp_deriv <- "^lambda|^com|^sp\\[|^Daily\\.Est|^selected|^Common|^Specific|^obs\\.prob"

BUGS.samples <- WinBugs.out$BUGS.out$sims.array

BUGS.params <- grp_est
convergence.check <- function(BUGS.samples, BUGS.params, group = NA) {
  nm <- grep(BUGS.params, dimnames(BUGS.samples)[[3]], value = TRUE, perl = TRUE)
  if (length(nm) == 0) stop("no parameters matched: ", BUGS.params)
  
  sub <- BUGS.samples[, , nm, drop = FALSE]          # keep 3-D
  d   <- posterior::as_draws_array(sub)              # iter x chain x var already
  
  s <- posterior::summarise_draws(d, "mean", "sd", "rhat", "ess_bulk", "ess_tail")
  s$group    <- group
  s$constant <- is.na(s$rhat)                        # no variation across draws
  s
}

group.summary <- function(s) {
  data.frame(
    group           = s$group[1],
    n_param         = nrow(s),
    n_constant      = sum(s$constant),
    n_rhat_over_101 = sum(s$rhat > 1.01, na.rm = TRUE),
    max_rhat        = max(s$rhat, na.rm = TRUE),
    min_ess_bulk    = min(s$ess_bulk, na.rm = TRUE),
    median_ess_bulk = median(s$ess_bulk, na.rm = TRUE)
  )
}


conv.grp_est <- convergence.check(BUGS.samples = BUGS.samples,
                              BUGS.params = grp_est, group = "Estimand")

conv.grp_cont<- convergence.check(BUGS.samples = BUGS.samples,
                                  BUGS.params = grp_cont, group = "Continuous")

conv.grp_disc<- convergence.check(BUGS.samples = BUGS.samples,
                                  BUGS.params = grp_disc, group = "Discrete")

conv.grp_deriv<- convergence.check(BUGS.samples = BUGS.samples,
                                  BUGS.params = grp_deriv, group = "Derived")

do.call(rbind, lapply(list(conv.grp_est, conv.grp_cont, 
                           conv.grp_disc, conv.grp_deriv), group.summary))

all_nm  <- dimnames(BUGS.samples)[[3]]
matched <- unlist(lapply(list(grp_est, grp_cont, grp_disc, grp_deriv),
                         function(p) grep(p, all_nm, value = TRUE, perl = TRUE)))
unique(gsub("\\[.*", "", setdiff(all_nm, matched)))
