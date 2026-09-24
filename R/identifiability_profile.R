#!/usr/bin/env Rscript
# ============================================================================
# Numerical identifiability profile — Hantavirus QSP Model
# ============================================================================
# The analytic ridge (R/identifiability_and_humoral.R) exhibits one direction
# in parameter space along which every calibration-visible quantity is
# invariant. This script asks the blunter question: taking the calibration
# data at face value, which parameters does a single mortality proportion
# actually pin?
#
# The calibration datum is 7 deaths in 75 placebo patients (Huggins 1991).
# Its exact binomial 95% interval is wide, so a range of parameter values is
# statistically consistent with it. For each parameter this script walks a
# grid, records the placebo mortality and the day-1 combination relative risk
# reduction, and marks the grid points whose mortality falls inside that
# interval. The span of the predicted treatment effect ACROSS THE CONSISTENT
# POINTS is the quantity of interest: it is what the trial cannot resolve.
#
# No re-fitting is attempted. An earlier version of this analysis re-fitted the
# infection-rate constant at each grid point, which is ill-posed: placebo
# mortality is flat in beta once infection establishes (it moves 0.14
# percentage points across a 33-fold range), so the root is not unique and a
# bracketed search returns whichever end of the plateau it lands on.
#
# Outputs:
#   outputs/identifiability_profile.csv   grid, mortality, RRR, consistency flag
#   outputs/identifiability_profile.png   RRR against each parameter
#   outputs/identifiability_profile.txt   ranked summary
#
# Usage (from project root): Rscript R/identifiability_profile.R
# ============================================================================

suppressPackageStartupMessages({
  library(deSolve); library(ggplot2); library(dplyr)
})

source("model/parameters.R")
source("model/hantavirus_qsp.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/simulate_trial.R")

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
pars0 <- get_parameters()

# The calibration datum and its exact binomial interval.
OBS_DEATHS <- 7
OBS_N      <- 75
ci  <- stats::binom.test(OBS_DEATHS, OBS_N)$conf.int
CI_LO <- 100 * ci[1]
CI_HI <- 100 * ci[2]

endpoints <- function(p, arm, t_start) {
  s <- simulate_patient(p, arm, t_start = t_start, t_end = 21)
  if (is.null(s)) return(NULL)
  extract_endpoints_from_simulation(s, p, "HFRS")
}

#' Placebo mortality (%) and day-1 combination RRR (%) for a parameter set.
evaluate <- function(p) {
  ep <- endpoints(p, "placebo", 0)
  ec <- endpoints(p, "combination", 1)
  if (is.null(ep) || is.null(ec)) return(c(NA_real_, NA_real_))
  m <- ep$mortality_prob
  if (!is.finite(m) || m <= 0) return(c(100 * m, NA_real_))
  c(100 * m, 100 * (m - ec$mortality_prob) / m)
}

# Mechanistic constants that carry the output. Endpoint-mapping constants are
# excluded so the profile probes mechanism rather than definition.
PROFILE <- list(
  beta  = list(label = "Infection rate (beta)",                mults = c(0.5, 0.6, 0.8, 1, 1.5, 2, 5, 10, 20)),
  p     = list(label = "Viral production rate (p)",            mults = c(0.1, 0.2, 0.5, 1, 2, 5, 10)),
  c     = list(label = "Viral clearance rate (c)",             mults = c(0.25, 0.5, 1, 2, 4)),
  k_PV  = list(label = "Viral drive on permeability (k_PV)",   mults = c(0.25, 0.5, 1, 2, 4)),
  k_KP  = list(label = "Permeability drive on kidney (k_KP)",  mults = c(0.25, 0.5, 1, 2, 4)),
  V_ref = list(label = "Permeability reference load (V_ref)",  mults = c(0.1, 0.2, 0.5, 1, 2, 5, 10))
)

rows <- list()
for (nm in names(PROFILE)) {
  spec <- PROFILE[[nm]]
  message(sprintf("Profiling %s ...", nm))
  for (mult in spec$mults) {
    q <- pars0
    q[[nm]] <- pars0[[nm]] * mult
    v <- evaluate(q)
    rows[[length(rows) + 1]] <- data.frame(
      parameter = nm, label = spec$label, multiplier = mult,
      value = q[[nm]], placebo_mortality = v[1], combo_RRR = v[2],
      consistent = !is.na(v[1]) && v[1] >= CI_LO && v[1] <= CI_HI,
      stringsAsFactors = FALSE)
  }
}
prof <- do.call(rbind, rows)
utils::write.csv(prof, "outputs/identifiability_profile.csv", row.names = FALSE)
message("  Saved: outputs/identifiability_profile.csv")

# --- Figure -----------------------------------------------------------------
plotd <- prof[!is.na(prof$combo_RRR), ]
g <- ggplot(plotd, aes(x = multiplier, y = combo_RRR)) +
  geom_line(linewidth = 0.8, colour = "grey55") +
  geom_point(aes(colour = consistent, shape = consistent), size = 2.6) +
  facet_wrap(~ label, scales = "free_x", ncol = 3) +
  scale_x_log10() +
  scale_colour_manual(values = c(`TRUE` = "#1f78b4", `FALSE` = "#bdbdbd"),
                      labels = c(`TRUE` = "yes", `FALSE` = "no"),
                      name = "Consistent with 7/75") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "yes", `FALSE` = "no"),
                     name = "Consistent with 7/75") +
  labs(x = "Parameter value, as a multiple of the calibrated value",
       y = "Day-1 combination relative risk reduction (%)",
       title = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave("outputs/identifiability_profile.png", g, width = 11, height = 7, dpi = 300)
message("  Saved: outputs/identifiability_profile.png")

# --- Summary ----------------------------------------------------------------
sink("outputs/identifiability_profile.txt")
cat("NUMERICAL IDENTIFIABILITY PROFILE\n")
cat("Generated by R/identifiability_profile.R\n")
cat("Date:", format(Sys.Date()), "\n")
cat(strrep("=", 78), "\n\n")
cat(sprintf("Calibration datum: %d deaths in %d placebo patients (%.2f%%)\n",
            OBS_DEATHS, OBS_N, 100 * OBS_DEATHS / OBS_N))
cat(sprintf("Exact binomial 95%% interval: %.2f%% to %.2f%%\n", CI_LO, CI_HI))
cat(sprintf("Day-1 combination RRR at the calibrated parameters: %.1f%%\n\n",
            evaluate(pars0)[2]))
cat("For each parameter: the range of the predicted treatment effect across\n")
cat("the grid points whose placebo mortality falls inside that interval, i.e.\n")
cat("across values the calibration datum cannot distinguish.\n\n")
cat(sprintf("%-8s %-38s %7s %9s %9s %9s\n",
            "param", "description", "n consis", "RRR min", "RRR max", "span"))
summ <- plotd %>% filter(consistent) %>% group_by(parameter, label) %>%
  summarise(n = dplyr::n(), lo = min(combo_RRR), hi = max(combo_RRR), .groups = "drop") %>%
  mutate(span = hi - lo) %>% arrange(desc(span))
for (i in seq_len(nrow(summ))) {
  r <- summ[i, ]
  cat(sprintf("%-8s %-38s %7d %9.1f %9.1f %9.1f\n",
              r$parameter, r$label, r$n, r$lo, r$hi, r$span))
}
cat("\nFull grid, including points the datum excludes:\n\n")
cat(sprintf("%-8s %9s %14s %11s %6s\n", "param", "multiple", "placebo mort%", "RRR%", "consis"))
for (i in seq_len(nrow(prof))) {
  r <- prof[i, ]
  cat(sprintf("%-8s %9.3g %14s %11s %6s\n", r$parameter, r$multiplier,
              if (is.na(r$placebo_mortality)) "-" else sprintf("%.3f", r$placebo_mortality),
              if (is.na(r$combo_RRR)) "-" else sprintf("%.1f", r$combo_RRR),
              if (r$consistent) "yes" else "no"))
}
sink()
message("  Saved: outputs/identifiability_profile.txt")
message("Done.")
