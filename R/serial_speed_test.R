# Serial test with 5 patients
source("model/hantavirus_qsp.R")
source("model/parameters.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/virtual_population.R")
source("R/simulate_trial.R")

pars <- get_parameters()
pop <- generate_virtual_population(N = 50, seed = 42)

n <- 5
t0 <- Sys.time()
for (i in 1:n) {
  row <- pop[i, ]
  pars_i <- apply_vpop_to_model(row, pars)
  t_treat <- row$time_to_treatment_days
  cat(sprintf("\nPatient %d: t_treat=%.2f, eGFR=%.1f, immune=%.2f\n",
              i, t_treat, row$eGFR_mL_min, row$immune_strength))
  
  t_start <- Sys.time()
  sim <- simulate_patient(pars_i, "placebo", t_start = 0, t_end = 21, dt = 0.5)
  elapsed <- difftime(Sys.time(), t_start, units = "secs")
  
  if (is.null(sim)) {
    cat(sprintf("  FAILED! Elapsed: %.1fs\n", elapsed))
  } else {
    ep <- extract_endpoints_from_simulation(sim, pars_i, syndrome = "HFRS")
    cat(sprintf("  V_peak=%.0f I_peak=%.0f mort=%.4f time=%.1fs\n",
                max(sim$V), max(sim$I), ep$mortality_prob, elapsed))
  }
}

total <- difftime(Sys.time(), t0, units = "secs")
cat(sprintf("\nTotal: %.1fs for %d patients (%.1f s/patient)\n", total, n, total/n))
