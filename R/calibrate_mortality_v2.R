#' Mortality recalibration v2 — fix viral dynamics severity
#'
#' Root cause: V_peak ~4,285 copies/mL is far too low. With I_ref=1e6,
#' the cytokine term k_CI * I/I_ref is negligible (~0.001/day). The
#' cytokine cascade never activates, so organ injury is trivial.
#'
#' Strategy: Increase beta (infection rate) and p (viral production)
#' to get V_peak into the 1e4-1e6 range. Then tune organ injury params.

source("model/hantavirus_qsp.R")
source("model/parameters.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/simulate_trial.R")
source("R/virtual_population.R")

run_single <- function(pars, t_end = 60) {
  sim <- simulate_patient(pars, "placebo", t_start = 0, t_end = t_end)
  if (is.null(sim)) return(NULL)
  pk  <- extract_peak_biomarkers(sim)
  ep  <- extract_endpoints_from_simulation(sim, pars, syndrome = "HFRS")
  c(
    K_peak         = max(sim$K, na.rm = TRUE),
    L_peak         = max(sim$L, na.rm = TRUE),
    C_pro_peak     = max(sim$C_pro, na.rm = TRUE),
    V_peak         = max(sim$V, na.rm = TRUE),
    I_peak         = max(sim$I, na.rm = TRUE),
    P_peak         = max(sim$P, na.rm = TRUE),
    F_I_peak       = max(sim$F_I, na.rm = TRUE),
    NK_peak        = max(sim$NK, na.rm = TRUE),
    F_II_peak      = max(sim$F_II, na.rm = TRUE),
    CD8_E_peak     = max(sim$CD8_E, na.rm = TRUE),
    mortality_prob = ep$mortality_prob
  )
}

# --- Phase 1: Sweep beta and p to find viral severity ---
cat("=== Phase 1: Viral severity sweep ===\n")

beta_vals <- c(4.6e-8, 7e-8, 1e-7, 1.5e-7, 2e-7, 3e-7, 5e-7)
p_vals    <- c(100, 200, 500, 1000, 2000)

phase1 <- expand.grid(beta = beta_vals, p = p_vals)
phase1$V_peak <- NA
phase1$I_peak <- NA
phase1$C_pro_peak <- NA
phase1$K_peak <- NA
phase1$L_peak <- NA
phase1$P_peak <- NA
phase1$F_I_peak <- NA
phase1$mortality <- NA
phase1$R0 <- NA

for (i in seq_len(nrow(phase1))) {
  pars <- get_parameters()
  pars$beta <- phase1$beta[i]
  pars$p    <- phase1$p[i]
  phase1$R0[i] <- pars$beta * pars$T_0 * pars$p / pars$c

  res <- run_single(pars, t_end = 60)
  if (!is.null(res)) {
    phase1$V_peak[i]     <- res["V_peak"]
    phase1$I_peak[i]     <- res["I_peak"]
    phase1$C_pro_peak[i] <- res["C_pro_peak"]
    phase1$K_peak[i]     <- res["K_peak"]
    phase1$L_peak[i]     <- res["L_peak"]
    phase1$P_peak[i]     <- res["P_peak"]
    phase1$F_I_peak[i]   <- res["F_I_peak"]
    phase1$mortality[i]  <- res["mortality_prob"]
  }
  cat(sprintf("  beta=%.1e p=%d R0=%.2f V_peak=%.0f I_peak=%.0f C_pro=%.4f K=%.3f mort=%.4f\n",
              phase1$beta[i], phase1$p[i], phase1$R0[i],
              phase1$V_peak[i], phase1$I_peak[i],
              phase1$C_pro_peak[i], phase1$K_peak[i], phase1$mortality[i]))
}

cat("\n=== Phase 1 summary ===\n")
print(phase1[order(-phase1$V_peak), ])

# --- Phase 2: With best beta/p, sweep organ injury params ---
cat("\n=== Phase 2: Organ injury calibration ===\n")

# Pick beta/p that give reasonable V_peak (1e5-1e6 range)
best_idx <- which(phase1$V_peak > 50000 & phase1$V_peak < 5e6)
if (length(best_idx) == 0) {
  best_idx <- which.max(phase1$V_peak)
}
best_row <- phase1[best_idx[1], ]
cat(sprintf("Using beta=%.2e p=%d (V_peak=%.0f, R0=%.2f)\n",
            best_row$beta, best_row$p, best_row$V_peak, best_row$R0))

# Now sweep organ injury parameters
set.seed(123)
n_trials <- 200

ranges2 <- list(
  k_KP     = c(0.1, 5.0),     # renal injury per permeability
  k_LP     = c(0.2, 10.0),    # lung injury per permeability
  k_PV     = c(0.5, 10.0),    # viral permeability drive
  k_PC     = c(0.1, 5.0),     # cytokine permeability drive
  k_CC     = c(0.5, 5.0),     # cytokine auto-amplification
  K_CC     = c(20, 500),      # auto-amplification half-max
  I_ref    = c(1e3, 1e5),     # reference I (lower = more sensitive)
  C_ref    = c(100, 1000),    # reference C_pro
  K50      = c(1, 30),        # renal half-max for mortality
  C50      = c(20, 200),      # cytokine half-max for mortality
  w_K_HFRS = c(0.05, 0.25),   # renal mortality weight
  w_C_HFRS = c(0.02, 0.15),   # cytokine mortality weight
  k_P_CD8  = c(0.002, 0.05),  # CD8 bystander damage
  k_CI     = c(1.0, 20.0),    # cytokine production by I
  d_C      = c(1.0, 10.0),    # cytokine decay (slower = higher AUC)
  d_K      = c(0.1, 1.0),     # renal recovery (slower = higher AUC)
  d_L      = c(0.1, 1.0),     # lung recovery
  I_gate   = c(100, 50000)    # auto-amplification gate
)

results2 <- data.frame()
target <- 0.094

for (i in 1:n_trials) {
  pars <- get_parameters()
  pars$beta <- best_row$beta
  pars$p    <- best_row$p
  pars$symptom_onset_day <- 5  # keep for viral peak timing

  for (pname in names(ranges2)) {
    lo <- ranges2[[pname]][1]
    hi <- ranges2[[pname]][2]
    val <- exp(runif(1, log(lo), log(hi)))
    pars[[pname]] <- val
  }

  res <- run_single(pars, t_end = 60)
  if (is.null(res)) next

  row <- data.frame(trial = i)
  for (pname in names(ranges2)) {
    row[[pname]] <- pars[[pname]]
  }
  for (nm in names(res)) {
    row[[nm]] <- res[nm]
  }
  results2 <- rbind(results2, row)

  if (i %% 20 == 0) {
    morts <- results2$mortality_prob
    cat(sprintf("  %d/%d: mort range [%.4f, %.4f], median=%.4f, n>5%%=%d\n",
                i, n_trials, min(morts), max(morts),
                median(morts), sum(morts > 0.05)))
  }
}

results2$mortality_error <- abs(results2$mortality_prob - target)
results2 <- results2[order(results2$mortality_error), ]

cat(sprintf("\nCompleted %d trials\n", nrow(results2)))
cat(sprintf("Mortality range: [%.4f, %.4f]\n",
            min(results2$mortality_prob), max(results2$mortality_prob)))

cat("\n=== Top 10 parameter sets ===\n")
top10 <- head(results2, 10)
print(top10[, c("trial", "mortality_prob", "V_peak", "I_peak", "C_pro_peak",
                 "K_peak", "L_peak", "P_peak", "CD8_E_peak")])

cat("\n=== Parameter values for top 10 ===\n")
print(top10[, c("trial", "mortality_prob", "k_KP", "k_LP", "k_PV", "k_PC",
                 "k_CC", "K_CC", "I_ref", "C_ref", "K50", "C50",
                 "w_K_HFRS", "w_C_HFRS", "k_P_CD8", "k_CI", "d_C", "d_K", "d_L", "I_gate")])

write.csv(results2, "outputs/mortality_calibration_v2.csv", row.names = FALSE)

# --- Recommend ---
best <- top10[1, ]
pars_base <- get_parameters()
cat("\n=== Recommended changes ===\n")
cat(sprintf("beta:       %.2e (was %.2e)\n", best_row$beta, pars_base$beta))
cat(sprintf("p:          %.0f (was %.0f)\n", best_row$p, pars_base$p))
for (pname in names(ranges2)) {
  cat(sprintf("%-12s %.4f (was %.4f)\n", paste0(pname, ":"), best[[pname]], pars_base[[pname]]))
}
cat(sprintf("\nAchieved mortality: %.4f (target: %.4f)\n", best$mortality_prob, target))
