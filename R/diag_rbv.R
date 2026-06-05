source('R/simulate_trial.R')
source('R/pd_models.R')
source('R/pk_models.R')
source('model/hantavirus_qsp.R')
source('model/parameters.R')
source('R/virtual_population.R')
pars <- get_parameters()

cat("=== I_gate =", pars$I_gate, "===\n\n")

for (arm in c("placebo", "ribavirin", "favipiravir", "combination")) {
  t_start <- if (arm == "placebo") 0 else 1
  s <- simulate_patient(pars, arm, t_start = t_start, t_end = 21)
  
  V_max <- max(s$V)
  I_max <- max(s$I)
  C_max <- max(s$C_pro)
  K_max <- max(s$K)
  L_max <- max(s$L)
  
  ep <- extract_endpoints_from_simulation(s, pars, "HFRS")
  
  cat(sprintf("%-15s V_max=%12.0f  I_max=%10.1f  C_pro_max=%8.1f  K=%7.2f  L=%7.2f  mort=%.4f\n",
              arm, V_max, I_max, C_max, K_max, L_max, ep$mortality_prob))
  cat(sprintf("                gate exceeded: %s  (I_gate=%d)\n",
              I_max > pars$I_gate, pars$I_gate))
  
  # Find when auto-amplification first triggers
  gate_crossed <- which(s$I > pars$I_gate)[1]
  if (!is.na(gate_crossed)) {
    cat(sprintf("                auto-amp starts at day %.1f (I=%.0f, C_pro=%.1f, C_RBV=%.1f)\n",
                s$time[gate_crossed], s$I[gate_crossed],
                s$C_pro[gate_crossed], ifelse("C_RBV" %in% names(s), s$C_RBV[gate_crossed], 0)))
  } else {
    cat("                auto-amp never triggers (I never exceeds gate)\n")
  }
  cat("\n")
}

# Detailed trace: ribavirin at day 1 vs day 4
cat("=== Ribavirin day-by-day detail ===\n\n")
for (day in c(1, 2, 3, 4, 5, 6, 7)) {
  s <- simulate_patient(pars, "ribavirin", t_start = day, t_end = 21)
  ep <- extract_endpoints_from_simulation(s, pars, "HFRS")
  gate_crossed <- which(s$I > pars$I_gate)[1]
  
  cat(sprintf("  RBV day %d: C_pro_max=%7.1f  K=%6.2f  L=%6.2f  mort=%.4f",
              day, max(s$C_pro), max(s$K), max(s$L), ep$mortality_prob))
  if (!is.na(gate_crossed)) {
    cat(sprintf("  | gate crossed at day %.1f, C_RBV_at_gate=%.1f\n",
                s$time[gate_crossed], s$C_RBV[gate_crossed]))
  } else {
    cat("  | gate never crossed\n")
  }
}

cat("\n=== Ribavirin PK profile (Day 1 start) ===\n")
s <- simulate_patient(pars, "ribavirin", t_start = 1, t_end = 21)
idx <- which(s$time >= 1 & s$time <= 8)
cat("  Day   C_RBV    E_RBV     E_immuno\n")
for (i in idx[seq(1, length(idx), by = 12)]) {
  C_RBV <- s$C_RBV[i]
  E_RBV <- (0.95 * C_RBV^1.5) / (1^1.5 + C_RBV^1.5)
  E_immuno <- 0.8 * C_RBV / (3 + C_RBV)
  cat(sprintf("  %.1f   %7.2f   %6.3f    %6.3f\n",
              s$time[i], C_RBV, E_RBV, E_immuno))
}
