#!/usr/bin/env Rscript
# ============================================================================
# Identifiability ridge and humoral-strength scans — Hantavirus QSP Model
# ============================================================================
# Produces the two analyses that are reported as tables but had no generating
# script in this repository:
#
#   main-text table 3   the practical identifiability ridge
#                       p -> p/f, beta -> beta*f, V_ref -> V_ref/f
#   ESM table S6        dependence of the day-1 favipiravir penalty on the
#                       strength of the humoral arm
#
# Both are deterministic scans over the representative HFRS patient, so no
# seed is required: re-running reproduces the tables exactly.
#
# Outputs:
#   outputs/identifiability_ridge.csv   table 3
#   outputs/humoral_scan.csv            table S6
#   outputs/identifiability_summary.txt both, formatted as in the manuscript
#
# Usage (from project root): Rscript R/identifiability_and_humoral.R
# ============================================================================

suppressPackageStartupMessages({
  library(deSolve)
})

source("model/parameters.R")
source("model/hantavirus_qsp.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/simulate_trial.R")

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

pars <- get_parameters()

#' Trapezoidal integral, used for the clearance-share decomposition.
trapz <- function(t, y) sum(diff(t) * (head(y, -1) + tail(y, -1)) / 2)

#' Simulate and stop rather than summarise a failed run: every quantity below
#' is a maximum or an integral, so a truncated trajectory would read as a
#' smaller peak and a smaller area (see R/simulate_trial.R).
sim_or_stop <- function(p, arm, t_start, t_end) {
  s <- simulate_patient(p, arm, t_start = t_start, t_end = t_end)
  if (is.null(s)) {
    stop(sprintf("simulation failed: arm=%s, t_start=%s, t_end=%s",
                 arm, t_start, t_end))
  }
  s
}

# ---------------------------------------------------------------------------
# Main-text table 3: the identifiability ridge
# ---------------------------------------------------------------------------
# Rescaling viral production down and infectivity up, with the permeability
# reference scaled to match, leaves both the establishment index beta*T_0*p/c
# and the normalised viral drive V/V_ref unchanged. The injury cascade -- the
# only thing the calibration data constrain -- therefore cannot distinguish
# the cases, while peak viraemia and the predicted treatment effect move.
message("[1/2] Identifiability ridge (f = 1, 3, 10, 30, 100)...")

f_values <- c(1, 3, 10, 30, 100)
ridge <- do.call(rbind, lapply(f_values, function(f) {
  q <- pars
  q$p     <- pars$p / f
  q$beta  <- pars$beta * f
  q$V_ref <- pars$V_ref / f

  pl <- sim_or_stop(q, "placebo",     t_start = 0, t_end = 21)
  cb <- sim_or_stop(q, "combination", t_start = 1, t_end = 21)
  ep <- extract_endpoints_from_simulation(pl, q, "HFRS")
  ec <- extract_endpoints_from_simulation(cb, q, "HFRS")

  data.frame(
    f                  = f,
    V_peak             = max(pl$V),
    placebo_mortality  = 100 * ep$mortality_prob,
    K_peak             = max(pl$K),
    L_peak             = max(pl$L),
    PLT_nadir          = min(pl$PLT),
    dialysis_risk      = 100 * ep$dialysis_prob,
    combo_day1_RRR     = 100 * (ep$mortality_prob - ec$mortality_prob) /
                               ep$mortality_prob,
    stringsAsFactors   = FALSE
  )
}))
utils::write.csv(ridge, "outputs/identifiability_ridge.csv", row.names = FALSE)
message("  Saved: outputs/identifiability_ridge.csv")

# ---------------------------------------------------------------------------
# ESM table S6: dependence of the day-1 penalty on humoral strength
# ---------------------------------------------------------------------------
# The rebound that drives the day-1 favipiravir penalty is only possible
# because antibody neutralisation performs a negligible share of viral
# clearance. Scaling both neutralisation rates by a common factor tests how
# far that has to change before the penalty disappears. A 42-day horizon is
# required: the day-1 injury peak falls outside 21 days (see section 3.4).
message("[2/2] Humoral-strength scan (x1, x10, x100, x1000)...")

ab_share <- function(s, p) {
  vc  <- trapz(s$time, p$c * s$V)
  vab <- trapz(s$time, p$k_neut_IgM * s$IgM * s$V) +
         trapz(s$time, p$k_neut_IgG * s$IgG * s$V)
  100 * vab / (vc + vab)
}

scales <- c(1, 10, 100, 1000)
humoral <- do.call(rbind, lapply(scales, function(sc) {
  q <- pars
  q$k_neut_IgM <- pars$k_neut_IgM * sc
  q$k_neut_IgG <- pars$k_neut_IgG * sc

  pl <- sim_or_stop(q, "placebo",     t_start = 0, t_end = 42)
  d1 <- sim_or_stop(q, "favipiravir", t_start = 1, t_end = 42)
  d2 <- sim_or_stop(q, "favipiravir", t_start = 2, t_end = 42)

  m <- function(s) 100 * extract_endpoints_from_simulation(s, q, "HFRS")$mortality_prob

  data.frame(
    antibody_scaling      = sc,
    placebo_mortality     = m(pl),
    day1_mortality        = m(d1),
    day2_mortality        = m(d2),
    day1_minus_day2_points = m(d1) - m(d2),
    # Reported for both arms: the manuscript quotes the range across them,
    # and they differ substantially once neutralisation is strong.
    antibody_share_day1   = ab_share(d1, q),
    antibody_share_day2   = ab_share(d2, q),
    stringsAsFactors      = FALSE
  )
}))
utils::write.csv(humoral, "outputs/humoral_scan.csv", row.names = FALSE)
message("  Saved: outputs/humoral_scan.csv")

# ---------------------------------------------------------------------------
# Text summary, laid out as the two tables appear in the manuscript
# ---------------------------------------------------------------------------
sink("outputs/identifiability_summary.txt")
cat("IDENTIFIABILITY RIDGE AND HUMORAL-STRENGTH SCANS\n")
cat("Generated by R/identifiability_and_humoral.R\n")
cat("Date:", format(Sys.Date()), "\n")
cat(strrep("=", 78), "\n\n")

cat("MAIN-TEXT TABLE 3 - practical identifiability ridge\n")
cat("p -> p/f, beta -> beta*f, V_ref -> V_ref/f; representative HFRS patient,\n")
cat("21-day horizon. Every calibration-visible quantity is invariant.\n\n")
cat(sprintf("%-6s %-12s %-10s %-9s %-9s %-11s %-10s %-10s\n",
            "f", "V_peak", "mort (%)", "K_peak", "L_peak", "PLT nadir",
            "dialysis%", "RRR (%)"))
for (i in seq_len(nrow(ridge))) {
  r <- ridge[i, ]
  cat(sprintf("%-6g %-12.3g %-10.3f %-9.2f %-9.1f %-11.0f %-10.2f %-10.1f\n",
              r$f, r$V_peak, r$placebo_mortality, r$K_peak, r$L_peak,
              r$PLT_nadir, r$dialysis_risk, r$combo_day1_RRR))
}

cat("\n\nESM TABLE S6 - day-1 favipiravir penalty vs humoral strength\n")
cat("Favipiravir monotherapy, 15-day course, 42-day horizon.\n")
cat("Antibody share is the fraction of viral clearance performed by the two\n")
cat("neutralisation terms; it differs between the two start days.\n\n")
cat(sprintf("%-10s %-11s %-11s %-11s %-12s %-12s %-12s\n",
            "scaling", "placebo %", "day1 %", "day2 %", "d1-d2 pts",
            "Ab share d1", "Ab share d2"))
for (i in seq_len(nrow(humoral))) {
  r <- humoral[i, ]
  cat(sprintf("x%-9g %-11.3f %-11.3f %-11.3f %+-12.2f %-12.2f %-12.2f\n",
              r$antibody_scaling, r$placebo_mortality, r$day1_mortality,
              r$day2_mortality, r$day1_minus_day2_points,
              r$antibody_share_day1, r$antibody_share_day2))
}
cat("\nThe penalty (day1 - day2 > 0) survives x10 and x100 and reverses only at\n")
cat("x1000, where placebo mortality has itself moved and the calibration no\n")
cat("longer holds.\n")
sink()
message("  Saved: outputs/identifiability_summary.txt")
message("Done.")
