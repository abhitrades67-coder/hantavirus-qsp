# Quick pipeline test with beta=7e-8
source("model/hantavirus_qsp.R")
source("model/parameters.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/virtual_population.R")
source("R/simulate_trial.R")
source("R/analysis.R")

pars <- get_parameters()
cat(sprintf("beta=%.1e, R0=%.2f\n", pars$beta, pars$beta*pars$T_0*pars$p/pars$c))

# Small virtual population
pop <- generate_virtual_population(N = 100, seed = 42)

# Run trial with 50 patients per arm
cat("\nRunning trial...\n")
trial <- run_virtual_trial(
  pop = pop, pars = pars,
  arms = c("placebo", "ribavirin", "favipiravir", "combination"),
  n_patients = 50, t_end = 21, dt = 0.5, seed = 123
)

cat(sprintf("\nCompleted: %d patient-simulations\n", nrow(trial$endpoints)))

# Summarize mortality by arm
for (arm in c("placebo", "ribavirin", "favipiravir", "combination")) {
  idx <- trial$endpoints$arm == arm
  mort <- trial$endpoints$mortality_prob[idx]
  cat(sprintf("%s: n=%d, mortality mean=%.4f, median=%.4f, range=[%.4f, %.4f]\n",
              arm, sum(idx), mean(mort), median(mort), min(mort), max(mort)))
}

# Write trial results
dir.create("outputs", showWarnings = FALSE)
write.csv(trial$endpoints, "outputs/quick_trial_endpoints.csv", row.names = FALSE)
cat("\nSaved: outputs/quick_trial_endpoints.csv\n")
