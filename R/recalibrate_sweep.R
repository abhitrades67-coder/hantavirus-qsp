#' Calibration sweep: platelet consumption + endpoint thresholds
#'
#' Representative HFRS patient. Fixes k_PLT_prod = 25000 (baseline bug fix),
#' then sweeps k_PLT_cons to find a clinically realistic placebo nadir, and
#' reports K_peak/L_peak/mortality so endpoint thresholds can be set.

suppressWarnings(suppressMessages({
  source("model/hantavirus_qsp.R")
  source("model/parameters.R")
  source("R/pk_models.R")
  source("R/pd_models.R")
  source("R/simulate_trial.R")
}))

sim_metrics <- function(overrides, arm, tstart) {
  pars <- get_parameters()
  for (nm in names(overrides)) pars[[nm]] <- overrides[[nm]]
  sim <- simulate_patient(pars, arm, t_start = tstart, t_end = 21, dt = 0.5)
  pk <- extract_peak_biomarkers(sim)
  ep <- extract_endpoints_from_simulation(sim, pars, syndrome = "HFRS")
  list(PLT_nadir = pk$PLT_nadir, K_peak = pk$K_peak, L_peak = pk$L_peak,
       mortality = ep$mortality_prob)
}

cat("\n=== Platelet consumption sweep (k_PLT_prod = 25000, placebo) ===\n")
cat(sprintf("%-12s %12s %10s %10s %10s\n", "k_PLT_cons", "PLT_nadir", "K_peak", "L_peak", "mort%"))
for (kc in c(0.01, 0.005, 0.003, 0.002, 0.0015, 0.001)) {
  m <- sim_metrics(list(k_PLT_prod = 25000, k_PLT_cons = kc), "placebo", 0)
  cat(sprintf("%-12.4f %12.0f %10.1f %10.1f %10.2f\n",
              kc, m$PLT_nadir, m$K_peak, m$L_peak, 100*m$mortality))
}
