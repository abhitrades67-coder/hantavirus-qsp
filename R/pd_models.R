#' Pharmacodynamic Models — Antiviral Effect and Combination
#'
#' Emax models for Ribavirin and Favipiravir, Bliss independence,
#' synergy interaction, and clinical endpoint functions.
#'
#' @name pd_models
NULL

#' Ribavirin antiviral effect (Emax model)
#'
#' @param C_RBV Ribavirin plasma concentration (ug/mL)
#' @param Emax Maximum effect (0-1)
#' @param EC50 Half-maximal concentration (ug/mL)
#' @param gamma Hill coefficient
#' @return Effect fraction (0 to Emax)
#' @export
ribavirin_effect <- function(C_RBV, Emax = 0.95, EC50 = 5, gamma = 1.5) {
  (Emax * C_RBV^gamma) / (EC50^gamma + C_RBV^gamma)
}

#' Favipiravir antiviral effect via RTP metabolite (Emax model)
#'
#' @param C_RTP Favipiravir-RTP intracellular concentration (ug/mL)
#' @param Emax Maximum effect (0-1)
#' @param EC50 Half-maximal concentration (ug/mL)
#' @param gamma Hill coefficient
#' @return Effect fraction (0 to Emax)
#' @export
favipiravir_effect <- function(C_RTP, Emax = 0.98, EC50 = 0.5, gamma = 1.2) {
  (Emax * C_RTP^gamma) / (EC50^gamma + C_RTP^gamma)
}

#' Bliss independence combination effect
#'
#' E_Bliss = 1 - (1 - E1) * (1 - E2)
#'
#' @param E1 Effect of drug 1
#' @param E2 Effect of drug 2
#' @return Combined effect under Bliss independence
#' @export
bliss_combination <- function(E1, E2) {
  1 - (1 - E1) * (1 - E2)
}

#' Synergy-adjusted combination effect
#'
#' psi_combo = psi * E1 * E2 (additional synergistic effect)
#' Total effect = E_Bliss + psi_combo (capped at 1)
#'
#' @param E1 Effect of drug 1
#' @param E2 Effect of drug 2
#' @param psi Synergy interaction coefficient
#' @return Combined effect with synergy
#' @export
synergy_combination <- function(E1, E2, psi = 0.1) {
  E_bliss <- bliss_combination(E1, E2)
  psi_term <- psi * E1 * E2
  min(E_bliss + psi_term, 1.0)
}

#' Compute clinical endpoint probabilities
#'
#' @param K Renal injury index (0-1 scale)
#' @param L Lung injury index (0-1 scale)
#' @param C_pro Pro-inflammatory cytokine burden
#' @param C_pro_max Maximum cytokine burden observed
#' @param K_auc AUC of renal injury (day*AU, default 0)
#' @param L_auc AUC of lung injury (day*AU, default 0)
#' @param C_auc AUC of pro-inflammatory cytokines (day*AU, default 0)
#' @param syndrome Syndrome phenotype: "HFRS" or "HCPS"
#' @param pars Parameter list with threshold and weight parameters
#' @return Named list: dialysis_prob, ecmo_prob, mortality_prob
#' @export
compute_clinical_endpoints <- function(K, L, C_pro, C_pro_max, syndrome, pars,
                                        K_auc = 0, L_auc = 0, C_auc = 0) {
  w_auc <- if ("w_auc" %in% names(pars)) pars$w_auc else 0.30

  if (syndrome == "HFRS") {
    renal_risk    <- pars$w_K_HFRS * K / (K + pars$K50)
    lung_risk     <- pars$w_L_HFRS * L / (L + pars$L50)
    cytokine_risk <- pars$w_C_HFRS * C_pro / (C_pro + pars$C50)
    renal_auc    <- pars$w_K_HFRS * K_auc / (K_auc + pars$K_auc50)
    lung_auc     <- pars$w_L_HFRS * L_auc / (L_auc + pars$L_auc50)
    cytokine_auc <- pars$w_C_HFRS * C_auc / (C_auc + pars$C_auc50)
  } else {
    # HCPS: lung-predominant
    renal_risk    <- pars$w_K_HCPS * K / (K + pars$K50)
    lung_risk     <- pars$w_L_HCPS * L / (L + pars$L50)
    cytokine_risk <- pars$w_C_HCPS * C_pro / (C_pro + pars$C50)
    renal_auc    <- pars$w_K_HCPS * K_auc / (K_auc + pars$K_auc50)
    lung_auc     <- pars$w_L_HCPS * L_auc / (L_auc + pars$L_auc50)
    cytokine_auc <- pars$w_C_HCPS * C_auc / (C_auc + pars$C_auc50)
  }

  peak_component  <- min(renal_risk + lung_risk + cytokine_risk, 1)
  auc_component   <- min(renal_auc + lung_auc + cytokine_auc, 1)
  mortality_prob  <- min((1 - w_auc) * peak_component + w_auc * auc_component, 1)

  # Dialysis and ECMO as secondary outputs (for reporting). These use
  # half-saturation constants DECOUPLED from the mortality K50/L50 so the
  # renal-replacement and cardiopulmonary-support rates can be calibrated to
  # clinically plausible values independently of the mortality model.
  K50_d <- if ("K50_dialysis" %in% names(pars)) pars$K50_dialysis else pars$K50
  L50_e <- if ("L50_ecmo" %in% names(pars)) pars$L50_ecmo else pars$L50
  dialysis_prob <- K / (K + K50_d)
  ecmo_prob     <- L / (L + L50_e)

  list(
    dialysis_prob  = dialysis_prob,
    ecmo_prob      = ecmo_prob,
    mortality_prob = mortality_prob
  )
}

#' Extract clinical endpoints from simulation output
#'
#' @param sim_output Data frame from deSolve simulation
#' @param pars Parameter list
#' @param syndrome Syndrome phenotype: "HFRS" or "HCPS"
#' @return Named list of endpoint probabilities over the acute simulation
#' @export
extract_endpoints_from_simulation <- function(sim_output, pars, syndrome = "HFRS") {
  K         <- max(sim_output$K, na.rm = TRUE)
  L         <- max(sim_output$L, na.rm = TRUE)
  C_pro     <- max(sim_output$C_pro, na.rm = TRUE)
  C_pro_max <- max(sim_output$C_pro, na.rm = TRUE)

  # AUC of organ injury — trapezoidal integration
  t  <- sim_output$time
  K_auc  <- auc_trapezoidal(t, sim_output$K)
  L_auc  <- auc_trapezoidal(t, sim_output$L)
  C_auc  <- auc_trapezoidal(t, sim_output$C_pro)

  compute_clinical_endpoints(K, L, C_pro, C_pro_max, syndrome, pars,
                             K_auc = K_auc, L_auc = L_auc, C_auc = C_auc)
}

# Internal trapezoidal AUC helper
auc_trapezoidal <- function(t, y) {
  n <- length(t)
  if (n < 2) return(0)
  auc <- 0
  for (i in seq_len(n - 1)) {
    auc <- auc + 0.5 * (y[i + 1] + y[i]) * (t[i + 1] - t[i])
  }
  auc
}

#' Peak biomarker extraction
#'
#' @param sim_output Data frame from deSolve simulation
#' @return Named list of peak values
#' @export
extract_peak_biomarkers <- function(sim_output) {
  list(
    V_peak      = max(sim_output$V, na.rm = TRUE),
    C_pro_peak  = max(sim_output$C_pro, na.rm = TRUE),
    P_peak      = max(sim_output$P, na.rm = TRUE),
    PLT_nadir   = min(sim_output$PLT, na.rm = TRUE),
    K_peak      = max(sim_output$K, na.rm = TRUE),
    L_peak      = max(sim_output$L, na.rm = TRUE),
    CD8_E_peak  = max(sim_output$CD8_E, na.rm = TRUE),
    CD4_peak    = max(sim_output$CD4, na.rm = TRUE),
    IgM_peak    = max(sim_output$IgM, na.rm = TRUE),
    IgG_peak    = max(sim_output$IgG, na.rm = TRUE),
    Hgb_drop_max = max(sim_output$Hgb_drop, na.rm = TRUE),
    time_V_peak  = sim_output$time[which.max(sim_output$V)],
    time_C_pro_peak = sim_output$time[which.max(sim_output$C_pro)]
  )
}

#' Area under the viral load curve
#'
#' @param sim_output Data frame from deSolve simulation
#' @return AUC of viral load (copies * day / mL)
#' @export
compute_viral_AUC <- function(sim_output) {
  # Trapezoidal rule
  t <- sim_output$time
  V <- sim_output$V
  n <- length(t)
  auc <- 0
  for (i in 2:n) {
    auc <- auc + 0.5 * (V[i] + V[i - 1]) * (t[i] - t[i - 1])
  }
  auc
}

#' Time to viral clearance (below detectable threshold)
#'
#' @param sim_output Data frame from deSolve simulation
#' @param threshold Detection threshold (default 100 copies/mL)
#' @return Time in days, or NA if never cleared
#' @export
time_to_clearance <- function(sim_output, threshold = 100) {
  cleared <- which(sim_output$V < threshold)
  if (length(cleared) == 0) return(NA_real_)
  sim_output$time[min(cleared)]
}
