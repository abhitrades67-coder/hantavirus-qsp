# Quick test: does the new 23-compartment model run at all?
source("model/parameters.R")
source("model/hantavirus_qsp.R")
source("R/virtual_population.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/simulate_trial.R")

pop <- generate_virtual_population(N = 5, seed = 42)
pars <- get_parameters()

for (arm in c("placebo", "ribavirin")) {
  cat(sprintf("\n=== %s ===\n", toupper(arm)))
  row <- pop[1, ]
  pars_i <- apply_vpop_to_model(row, pars)
  
  sim <- tryCatch(
    simulate_patient(pars_i, arm, t_start = 3, t_end = 21),
    error = function(e) { cat(sprintf("ERROR: %s\n", e$message)); NULL }
  )
  
  if (!is.null(sim)) {
    cat(sprintf("  Columns: %s\n", paste(names(sim), collapse=", ")))
    cat(sprintf("  V_peak: %.2e\n", max(sim$V)))
    cat(sprintf("  CD8_N peak: %.1f, CD8_E peak: %.1f\n", max(sim$CD8_N), max(sim$CD8_E)))
    cat(sprintf("  CD8_E day 14: %.1f, day 21: %.1f\n", 
        sim$CD8_E[which.min(abs(sim$time - 14))],
        sim$CD8_E[which.min(abs(sim$time - 21))]))
  }
}
