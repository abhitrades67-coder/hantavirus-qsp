#' Mortality recalibration after CD8_N baseline fix
#'
#' With CD8_E starting at 0 (baseline moved to CD8_N), the organ injury
#' pathways need to be strengthened to restore ~9.4% placebo mortality.
#' This script sweeps key parameters on single-patient simulations.

source("model/hantavirus_qsp.R")
source("model/parameters.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/simulate_trial.R")
source("R/virtual_population.R")

# --- Helper: run single patient and extract endpoints ---
run_single <- function(pars, t_end = 21) {
  sim <- simulate_patient(pars, "placebo", t_start = 0, t_end = t_end)
  if (is.null(sim)) return(NULL)

  pk  <- extract_peak_biomarkers(sim)
  ep  <- extract_endpoints_from_simulation(sim, pars, syndrome = "HFRS")

  # Get peak K, L, C_pro from the simulation
  K_peak  <- max(sim$K, na.rm = TRUE)
  L_peak  <- max(sim$L, na.rm = TRUE)
  C_peak  <- max(sim$C_pro, na.rm = TRUE)
  V_peak  <- max(sim$V, na.rm = TRUE)
  P_peak  <- max(sim$P, na.rm = TRUE)
  CD8_E_peak <- max(sim$CD8_E, na.rm = TRUE)

  c(
    K_peak         = K_peak,
    L_peak         = L_peak,
    C_pro_peak     = C_peak,
    V_peak         = V_peak,
    P_peak         = P_peak,
    CD8_E_peak     = CD8_E_peak,
    mortality_prob = ep$mortality_prob,
    dialysis_prob  = ep$dialysis_prob,
    ecmo_prob      = ep$ecmo_prob
  )
}

# --- Baseline diagnostic ---
cat("\n=== Baseline diagnostic (current parameters) ===\n")
pars_base <- get_parameters()
res_base <- run_single(pars_base, t_end = 60)
print(round(res_base, 4))

# --- Parameter sweep ---
cat("\n=== Parameter sweep ===\n")

# Key parameters to sweep:
# 1. k_CC: cytokine auto-amplification (was 0.8)
# 2. K_CC: auto-amplification half-max (was 200)
# 3. k_KP: renal injury per permeability (was 0.15)
# 4. k_LP: lung injury per permeability (was 0.2)
# 5. k_PV: viral permeability drive (was 0.5)
# 6. k_PC: cytokine permeability drive (was 0.1)
# 7. K50: renal half-max for mortality (was 30)
# 8. C50: cytokine half-max for mortality (was 300)
# 9. w_K_HFRS: renal mortality weight (was 0.12)
# 10. w_C_HFRS: cytokine mortality weight (was 0.05)
# 11. k_P_CD8: CD8 bystander permeability (was 0.005)
# 12. k_CI: cytokine prod by infected cells (was 2.0)

# Use a Latin hypercube approach — test combinations
set.seed(42)
n_trials <- 100

# Ranges for each parameter (log-uniform for wide exploration)
ranges <- list(
  k_CC        = c(0.5, 5.0),      # cytokine auto-amplification
  K_CC        = c(50, 500),       # auto-amplification half-max
  k_KP        = c(0.05, 2.0),     # renal injury per permeability
  k_LP        = c(0.1, 3.0),      # lung injury per permeability
  k_PV        = c(0.2, 5.0),      # viral permeability drive
  k_PC        = c(0.05, 2.0),     # cytokine permeability drive
  K50         = c(3, 50),         # renal half-max for mortality
  C50         = c(30, 500),       # cytokine half-max for mortality
  w_K_HFRS    = c(0.05, 0.30),    # renal mortality weight
  w_C_HFRS    = c(0.02, 0.15),    # cytokine mortality weight
  k_P_CD8     = c(0.001, 0.05),   # CD8 bystander damage
  k_CI        = c(0.5, 10.0)      # cytokine production by I
)

results <- data.frame()

for (i in 1:n_trials) {
  pars <- pars_base

  # Sample uniformly in log-space for each parameter
  for (pname in names(ranges)) {
    lo <- ranges[[pname]][1]
    hi <- ranges[[pname]][2]
    # Log-uniform sampling
    val <- exp(runif(1, log(lo), log(hi)))
    pars[[pname]] <- val
  }

  res <- run_single(pars, t_end = 60)
  if (is.null(res)) next

  row <- data.frame(
    trial = i,
    k_CC = pars$k_CC,
    K_CC = pars$K_CC,
    k_KP = pars$k_KP,
    k_LP = pars$k_LP,
    k_PV = pars$k_PV,
    k_PC = pars$k_PC,
    K50  = pars$K50,
    C50  = pars$C50,
    w_K_HFRS = pars$w_K_HFRS,
    w_C_HFRS = pars$w_C_HFRS,
    k_P_CD8 = pars$k_P_CD8,
    k_CI  = pars$k_CI,
    t(res)
  )

  results <- rbind(results, row)

  if (i %% 20 == 0) {
    cat(sprintf("  Completed %d/%d trials...\n", i, n_trials))
  }
}

cat(sprintf("\nCompleted %d trials\n", nrow(results)))

# --- Find trials with mortality closest to 9.4% ---
target <- 0.094
results$mortality_error <- abs(results$mortality_prob - target)
results <- results[order(results$mortality_error), ]

cat("\n=== Top 10 parameter sets (closest to 9.4% mortality) ===\n")
top10 <- head(results, 10)
print(top10[, c("trial", "mortality_prob", "K_peak", "L_peak", "C_pro_peak",
                 "V_peak", "P_peak", "CD8_E_peak")])

cat("\n=== Corresponding parameter values ===\n")
print(top10[, c("trial", "mortality_prob", "k_CC", "K_CC", "k_KP", "k_LP",
                 "k_PV", "k_PC", "K50", "C50", "w_K_HFRS", "w_C_HFRS",
                 "k_P_CD8", "k_CI")])

# Save results
write.csv(results, "outputs/mortality_calibration_sweep.csv", row.names = FALSE)
cat("\nSaved: outputs/mortality_calibration_sweep.csv\n")

# --- Recommend parameters ---
best <- top10[1, ]
cat("\n=== Recommended parameters ===\n")
cat(sprintf("k_CC:     %.3f  (was %.3f)\n", best$k_CC, pars_base$k_CC))
cat(sprintf("K_CC:     %.1f  (was %.1f)\n", best$K_CC, pars_base$K_CC))
cat(sprintf("k_KP:     %.3f  (was %.3f)\n", best$k_KP, pars_base$k_KP))
cat(sprintf("k_LP:     %.3f  (was %.3f)\n", best$k_LP, pars_base$k_LP))
cat(sprintf("k_PV:     %.3f  (was %.3f)\n", best$k_PV, pars_base$k_PV))
cat(sprintf("k_PC:     %.3f  (was %.3f)\n", best$k_PC, pars_base$k_PC))
cat(sprintf("K50:      %.1f  (was %.1f)\n", best$K50, pars_base$K50))
cat(sprintf("C50:      %.1f  (was %.1f)\n", best$C50, pars_base$C50))
cat(sprintf("w_K_HFRS: %.3f  (was %.3f)\n", best$w_K_HFRS, pars_base$w_K_HFRS))
cat(sprintf("w_C_HFRS: %.3f  (was %.3f)\n", best$w_C_HFRS, pars_base$w_C_HFRS))
cat(sprintf("k_P_CD8:  %.4f (was %.4f)\n", best$k_P_CD8, pars_base$k_P_CD8))
cat(sprintf("k_CI:     %.3f  (was %.3f)\n", best$k_CI, pars_base$k_CI))
cat(sprintf("\nAchieved mortality: %.4f (target: %.4f)\n", best$mortality_prob, target))
