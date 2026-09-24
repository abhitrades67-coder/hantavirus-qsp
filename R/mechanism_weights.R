#' Which mechanisms carry the model output at the calibrated parameter values.
#'
#' Section 3.2 and supplementary table S5 report how much each mechanism in the
#' model actually contributes to the reported outcome. Those numbers used to be
#' computed by hand; this script derives them from the model so that the table
#' regenerates with everything else.
#'
#' Two kinds of quantity are produced.
#'
#'   1. Rate-term shares. Each removal or production term is re-evaluated along
#'      the stored placebo trajectory and integrated over the 21 days, then
#'      expressed as a percentage of the total for its process. This is a
#'      decomposition of the flux, not a sensitivity: it says what the terms did
#'      on this trajectory, not what would happen if one were removed.
#'   2. Mortality shares. The mortality probability is a weighted sum of a peak
#'      and an integral component, each of which adds a renal, a pulmonary and a
#'      cytokine contribution, so it decomposes exactly. The decomposition is
#'      valid only while neither min() in compute_clinical_endpoints() binds,
#'      which is checked below.
#'
#' The ablation rows of table S5 come from outputs/ablation_summary.csv and are
#' not recomputed here.
#'
#' Run from the project root: Rscript R/mechanism_weights.R

source("model/parameters.R")
source("model/hantavirus_qsp.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/simulate_trial.R")

pars <- get_parameters()

# dt = 0.01 rather than the 0.1 used elsewhere: these are integrals of rate
# terms over a trajectory with a sharp viral peak, so the quadrature error
# matters more here than it does for a peak or an endpoint.
sim <- simulate_patient(pars, arm = "placebo", t_start = 0, t_end = 21, dt = 0.01)
if (is.null(sim)) stop("placebo simulation failed")

tt <- sim$time
trap <- function(y) sum(0.5 * (y[-1] + y[-length(y)]) * diff(tt))
pct  <- function(x, total) 100 * x / total

# ---- viral clearance ------------------------------------------------------
# dV loses virus to non-specific clearance and to IgM and IgG neutralisation.
clr_c   <- trap(pars$c * sim$V)
clr_IgM <- trap(pars$k_neut_IgM * sim$IgM * sim$V)
clr_IgG <- trap(pars$k_neut_IgG * sim$IgG * sim$V)
clr_tot <- clr_c + clr_IgM + clr_IgG

# ---- infected-cell clearance ----------------------------------------------
# No drug is present on the placebo arm, so E_RBV_immuno is zero and the
# NK boost factor is 1. CD8 killing uses the same smooth weighting as the ODE.
I <- sim$I
CD8_clearance <- (I / (200 + I)) * (pars$delta_CD8 * sim$CD8_E * I) +
                 (200 / (200 + I)) * (5.0 * sim$CD8_E) * I
ic_nat <- trap(pars$delta_natural * I)
ic_NK  <- trap(pars$delta_NK * sim$NK * I)
ic_CD8 <- trap(CD8_clearance)
ic_tot <- ic_nat + ic_NK + ic_CD8

# ---- type I interferon suppression of viral production --------------------
# The antiviral factor is 1 - E_max_IFN * F_I / (EC50_F + F_I); report the
# largest suppression reached at any point on the trajectory.
ifn_sup <- max(pars$E_max_IFN * sim$F_I / (pars$EC50_F + sim$F_I))
NSs_max <- max(sim$NSs)

# ---- cytokine auto-amplification ------------------------------------------
# term_auto is gated on I > I_gate and saturates at k_CC; report the largest
# fraction of that capacity used.
gate <- sim$I > pars$I_gate
auto_frac <- if (any(gate)) {
  max((sim$C_pro[gate]^2) / (pars$K_CC^2 + sim$C_pro[gate]^2))
} else 0

# ---- mortality decomposition ----------------------------------------------
K <- max(sim$K); L <- max(sim$L); C_pro <- max(sim$C_pro)
K_auc <- trap(sim$K); L_auc <- trap(sim$L); C_auc <- trap(sim$C_pro)
w_auc <- if ("w_auc" %in% names(pars)) pars$w_auc else 0.30

peak_parts <- c(renal = pars$w_K_HFRS * K / (K + pars$K50),
                lung  = pars$w_L_HFRS * L / (L + pars$L50),
                cyto  = pars$w_C_HFRS * C_pro / (C_pro + pars$C50))
auc_parts  <- c(renal = pars$w_K_HFRS * K_auc / (K_auc + pars$K_auc50),
                lung  = pars$w_L_HFRS * L_auc / (L_auc + pars$L_auc50),
                cyto  = pars$w_C_HFRS * C_auc / (C_auc + pars$C_auc50))
if (sum(peak_parts) > 1 || sum(auc_parts) > 1) {
  stop("a min() in compute_clinical_endpoints() is binding; the mortality ",
       "probability no longer decomposes additively and the shares below ",
       "would be wrong")
}
mort_parts <- (1 - w_auc) * peak_parts + w_auc * auc_parts
mort_total <- sum(mort_parts)

# Cross-check against the function the rest of the pipeline calls.
ep <- extract_endpoints_from_simulation(sim, pars, syndrome = "HFRS")
if (abs(ep$mortality_prob - mort_total) > 1e-9) {
  stop(sprintf("decomposition sums to %.10f but the model reports %.10f",
               mort_total, ep$mortality_prob))
}

out <- data.frame(
  quantity = c(rep("Mortality probability", 3),
               rep("Infected-cell clearance", 3),
               rep("Viral clearance", 3),
               "Type I interferon", "Cytokine auto-amplification"),
  component = c("Renal injury (peak and integral)", "Pulmonary injury",
                "Pro-inflammatory cytokines",
                "CD8+ killing", "Natural turnover", "NK-mediated",
                "Non-specific (c)", "IgG neutralisation", "IgM neutralisation",
                "Maximum suppression of viral production",
                "Fraction of capacity used"),
  percent = c(pct(mort_parts[["renal"]], mort_total),
              pct(mort_parts[["lung"]],  mort_total),
              pct(mort_parts[["cyto"]],  mort_total),
              pct(ic_CD8, ic_tot), pct(ic_nat, ic_tot), pct(ic_NK, ic_tot),
              pct(clr_c, clr_tot), pct(clr_IgG, clr_tot), pct(clr_IgM, clr_tot),
              100 * ifn_sup, 100 * auto_frac),
  stringsAsFactors = FALSE
)

dir.create("outputs", showWarnings = FALSE)
utils::write.csv(out, "outputs/mechanism_weights.csv", row.names = FALSE)

cat(sprintf("\nRepresentative placebo patient, 21 days, dt = 0.01\n"))
cat(sprintf("  mortality probability      %.5f\n", mort_total))
cat(sprintf("  peak NSs                   %.3g (suppression constant %.3g)\n",
            NSs_max, pars$K_NSs_sup))
cat(sprintf("  antibody share of clearance %.3f%%\n",
            pct(clr_IgM + clr_IgG, clr_tot)))
cat("\n")
for (i in seq_len(nrow(out))) {
  cat(sprintf("  %-28s %-40s %10.5g%%\n",
              out$quantity[i], out$component[i], out$percent[i]))
}

# ---- descriptive values for the representative placebo patient ------------
# Section 3.1 quotes these. They are peaks and nadirs rather than integrals,
# so they are taken at dt = 0.1, the step the rest of the pipeline reports at,
# and not at the finer step used for the integrals above. The two differ in
# the last digit of the platelet nadir and in where the viral peak is placed.
sim_std <- simulate_patient(pars, arm = "placebo", t_start = 0, t_end = 21,
                            dt = 0.1)
if (is.null(sim_std)) stop("placebo simulation at dt = 0.1 failed")
ep_std <- extract_endpoints_from_simulation(sim_std, pars, syndrome = "HFRS")
pk <- extract_peak_biomarkers(sim_std)
# Target-cell limitation: at the viral peak, production p*I balances clearance
# c*V, so V should equal p*I/c there. Report the gap between the two.
j <- which.max(sim_std$V)
quasi <- pars$p * sim_std$I[j] / pars$c
gap <- 100 * abs(sim_std$V[j] - quasi) / quasi

# Supplementary table S3 also compares the platelet nadir timing and the
# ribavirin haemoglobin decline, so record those too. The haemoglobin drop
# needs a treated arm: the placebo patient never receives the drug.
sim_rbv <- simulate_patient(pars, arm = "ribavirin", t_start = 1, t_end = 21,
                            dt = 0.1)
if (is.null(sim_rbv)) stop("day-1 ribavirin simulation failed")

desc <- data.frame(
  quantity = c("V_peak", "t_V_peak", "I_peak", "V_peak_vs_pI_over_c_pct",
               "P_peak", "t_P_peak", "K_peak", "t_K_peak", "L_peak",
               "PLT_nadir", "t_PLT_nadir", "mortality_pct", "dialysis_pct", "ecmo_pct",
               "Hgb_drop_ribavirin_day1", "C_RBV_peak_ribavirin_day1"),
  value = c(pk$V_peak, pk$time_V_peak, max(sim_std$I), gap,
            pk$P_peak, sim_std$time[which.max(sim_std$P)],
            pk$K_peak, sim_std$time[which.max(sim_std$K)], pk$L_peak,
            pk$PLT_nadir, sim_std$time[which.min(sim_std$PLT)],
            100 * ep_std$mortality_prob,
            100 * ep_std$dialysis_prob, 100 * ep_std$ecmo_prob,
            max(sim_rbv$Hgb_drop), max(sim_rbv$C_RBV)),
  stringsAsFactors = FALSE)
utils::write.csv(desc, "outputs/representative_placebo.csv", row.names = FALSE)
cat("\nRepresentative placebo patient\n")
for (i in seq_len(nrow(desc))) {
  cat(sprintf("  %-26s %.6g\n", desc$quantity[i], desc$value[i]))
}
cat("\nSaved: outputs/representative_placebo.csv\n")
cat("\nSaved: outputs/mechanism_weights.csv\n")
