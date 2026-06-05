#' Fast recalibration harness
#'
#' Simulates the representative HFRS patient (base parameters, no vpop covariates)
#' for placebo + day-1 active arms and prints the key calibration metrics.
#' Optionally applies a named list of parameter overrides so candidate
#' calibrations can be evaluated in seconds before a full pipeline run.
#'
#' Usage:
#'   Rscript R/recalibrate_harness.R                 # current parameters
#'   (edit OVERRIDES below to test candidate values)

suppressWarnings(suppressMessages({
  source("model/hantavirus_qsp.R")
  source("model/parameters.R")
  source("R/pk_models.R")
  source("R/pd_models.R")
  source("R/virtual_population.R")
  source("R/simulate_trial.R")
}))

# ---- Candidate parameter overrides (edit to test calibrations) -------------
OVERRIDES <- list(
  # k_PLT_prod   = 25000,
  # k_PLT_cons   = 0.0018,
  # K50_dialysis = 80,
  # L50_ecmo     = 160
)

run_rep <- function(overrides = list()) {
  pars <- get_parameters()
  for (nm in names(overrides)) pars[[nm]] <- overrides[[nm]]

  arms <- c("placebo", "ribavirin", "favipiravir", "combination")
  rows <- list()
  for (arm in arms) {
    tstart <- if (arm == "placebo") 0 else 1
    sim <- simulate_patient(pars, arm, t_start = tstart, t_end = 21, dt = 0.5)
    if (is.null(sim)) { next }
    pk <- extract_peak_biomarkers(sim)
    ep <- extract_endpoints_from_simulation(sim, pars, syndrome = "HFRS")
    rows[[arm]] <- data.frame(
      arm          = arm,
      V_peak       = pk$V_peak,
      V_AUC        = compute_viral_AUC(sim),
      PLT_nadir    = pk$PLT_nadir,
      K_peak       = pk$K_peak,
      L_peak       = pk$L_peak,
      C_pro_peak   = pk$C_pro_peak,
      dialysis     = ep$dialysis_prob,
      ecmo         = ep$ecmo_prob,
      mortality    = ep$mortality_prob,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

res <- run_rep(OVERRIDES)
options(width = 200)
cat("\n=== Representative HFRS patient (base parameters) ===\n")
print(format(res, digits = 4, scientific = TRUE), row.names = FALSE)
cat(sprintf("\nPlacebo: mortality=%.1f%%  PLT_nadir=%.0f/uL  K_peak=%.1f  L_peak=%.1f  dialysis=%.1f%%  ecmo=%.1f%%\n",
            100*res$mortality[res$arm=="placebo"], res$PLT_nadir[res$arm=="placebo"],
            res$K_peak[res$arm=="placebo"], res$L_peak[res$arm=="placebo"],
            100*res$dialysis[res$arm=="placebo"], 100*res$ecmo[res$arm=="placebo"]))
cat(sprintf("Ribavirin d1: mortality=%.1f%%   Combination d1: mortality=%.1f%%\n",
            100*res$mortality[res$arm=="ribavirin"], 100*res$mortality[res$arm=="combination"]))
