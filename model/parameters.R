#' Hantavirus QSP Model Parameters
#'
#' All parameters with default values, units, descriptions, sources, and
#' confidence ratings. Organised by module.
#'
#' Time base unit: **days** (PK rate constants converted from /h to /day
#' where needed).
#'
#' @return Named list of all model parameters
#' @export
get_parameters <- function() {
  pars <- list()

  # ===========================================================================
  # 1. VIRAL DYNAMICS
  # ===========================================================================
  pars$T_0       <- 1e6          # cells            # Initial susceptible target cells
  pars$I_0       <- 0            # cells            # Initial infected cells
  # V_0 increased from 100 to 1000 (2026-05-21): enables infection establishment
  # even at lower beta values. With V_0=100, beta must be ≥1e-7 to
  # overcome viral clearance before I builds up. V_0=1000 lowers the critical
  # I threshold proportionally and is clinically plausible (aerosolized rodent
  # excreta can deliver 10³-10⁴ virions; 1000 copies is ~0.003% of peak V).
  pars$V_0       <- 1000         # copies/mL        # Initial viral inoculum (increased for establishment)
  # beta recalibrated (2026-05-21, revised 2026-06): set to 1e-7 to ensure
  # >99% of virtual patients cross the infection bistable threshold (R0=3.33).
  # At R0 ~2.0-2.5, ~34% of patients with V0 variability burn out.
  # NOTE ON NOTATION: beta * T_0 * p / c = 3.33 is NOT the basic reproduction
  # number - it has units of 1/day. The dimensionless disease-free R0 for this
  # target-cell subsystem is beta*T_0*p/(c*delta) = 111 with delta =
  # delta_natural. The quantity below is an empirical viral-establishment
  # INDEX used for pre-screening; measured against this model its critical
  # value is about 1.5 (infection establishes at 2.0, fails at 1.33), not 1,
  # because early CD8 priming and IFN suppression, not delta_natural, set the
  # effective clearance rate. Do not report it as R0.
  # establishment index = beta * T_0 * p / c = 1e-7 * 1e6 * 100 / 3 = 3.33 /day.
  # Single-patient placebo mortality = 9.64% (cf. Huggins JID 1991: 9.4%).
  pars$beta      <- 1.0e-7       # mL/copy/day      # Infection rate constant (calibrated for 9.6% mortality)
  pars$p         <- 100          # copies/cell/day  # Viral production rate
  pars$c         <- 3            # /day             # Viral clearance rate
  pars$k_T_reg   <- 0.30         # /day             # Target cell regen rate (logistic, fast recovery)
  # 5-day infection-to-symptom burn-in: with the current beta=1e-7 and V_0=1000,
  # the viral peak occurs a few days after symptom onset, so the treatment
  # window (days 1-7) brackets peak viremia. (An earlier note here referenced a
  # superseded beta=4.6e-8 / "~day 7.5" calibration.)
  pars$symptom_onset_day <- 5    # days             # Infection-to-symptom-onset burn-in

  # ===========================================================================
  # 2. INNATE IMMUNE / CYTOKINE MODULE (expanded 6-variable)
  # ===========================================================================
  # NSs protein (hantavirus IFN antagonist)
  pars$k_NSs      <- 0.5         # /day             NSs production rate
  pars$d_NSs      <- 1.5         # /day             NSs decay (t1/2 ~ 11h)
  pars$K_NSs_sup  <- 10          # AU               NSs suppression half-max

  # Type I IFN (IFN-alpha/beta, antiviral)
  pars$k_FI       <- 5.0         # AU/day           IFN-I max production rate
  pars$K_FI       <- 1e5         # cells            IFN production saturation
  pars$d_F_I      <- 3.5         # /day             IFN-I decay (t1/2 ~ 4.8h)
  pars$E_max_IFN  <- 0.8         # dimensionless    Max IFN antiviral effect
  pars$EC50_F     <- 50          # AU               IFN half-max antiviral

  # NK cells
  pars$k_NK_I     <- 0.3         # AU/day           NK activation by infected cells
  pars$K_NK_I     <- 1e5         # cells            NK activation saturation (I)
  pars$k_NK_F     <- 0.2         # AU/day           NK activation by IFN-I
  pars$K_NK_F     <- 50          # AU               NK activation saturation (F_I)
  pars$d_NK       <- 0.7         # /day             NK decay (t1/2 ~ 24h)

  # Type II IFN (IFN-gamma, pro-inflammatory)
  pars$k_FII_NK   <- 0.2         # /day             IFN-gamma production by NK
  pars$d_F_II     <- 2.8         # /day             IFN-II decay (t1/2 ~ 6h)

  # Pro-inflammatory cytokines
  pars$k_CI       <- 2.0         # AU/day           Cytokine prod by infected cells
  pars$k_CF2      <- 1.5         # AU/day           Cytokine amplification by IFN-II
  pars$k_CC       <- 0.8         # /day             Cytokine auto-amplification max (calibrated)
  pars$K_CC       <- 200         # AU               Auto-amplification half-max (calibrated)
  pars$d_C        <- 5.5         # /day             Cytokine decay (t1/2 ~ 3h)

  # Anti-inflammatory cytokines (IL-10-like)
  pars$k_anti     <- 0.5         # AU/day           IL-10 production rate
  pars$K_anti     <- 100         # AU               IL-10 production saturation
  pars$d_anti     <- 0.5         # /day             IL-10 decay (t1/2 ~ 33h)
  pars$k_suppress <- 1.0         # /day             Max IL-10 suppression rate
  pars$K_suppress <- 30          # AU               IL-10 suppression potency

  # Reference values for normalization
  pars$I_ref      <- 1e6         # cells            Reference infected cells
  pars$F_I_ref    <- 100         # AU               Reference Type I IFN
  pars$F_II_ref   <- 50          # AU               Reference Type II IFN
  pars$C_ref      <- 500         # AU               Reference cytokines
  pars$V_ref      <- 1e5         # copies/mL        Reference viral load

  # Infected cell clearance (non-cytopathic)
  pars$delta_natural <- 0.03     # /day             Natural turnover (t1/2 ~ 23d)
  pars$delta_NK      <- 5e-6     # /(AU*day)        NK-mediated clearance

  # Auto-amplification gate threshold
  pars$I_gate     <- 10000       # cells            Gate for cytokine auto-amplification

  # ===========================================================================
  # 2b. ADAPTIVE IMMUNE MODULE (CD4, CD8, IgM, IgG)
  #
  # Literature-grounded parameters. Adaptive immunity develops over days 5-14
  # and provides delayed viral control. The system is intentionally modest
  # during the acute phase (first 21 days), consistent with hantavirus being
  # primarily controlled by innate mechanisms early on.
  # ===========================================================================
  # CD4+ helper T cells
  pars$k_CD4_I    <- 0.15          # AU/day           CD4 activation by infected cells
  pars$K_CD4_I    <- 500           # cells            CD4 activation saturation (I) (lowered to match CD8 priming sensitivity)
  pars$k_CD4_F    <- 0.10          # AU/day           CD4 activation by IFN-I
  pars$K_CD4_F    <- 50            # AU               CD4 activation saturation (F_I)
  pars$d_CD4      <- 0.10          # /day             CD4 decay (t1/2 ~ 7d, effector)

  # CD8+ cytotoxic T cells: DELAYED activation via maturation chain
  # CD8_N (primed) → CD8_E (effector) with ~7 day delay
  # (PMID: 19072554, 18814258, 21525363)
  pars$k_CD8_I    <- 5.0           # AU/day           CD8 priming by infected cells (upscaled for delay chain)
  pars$K_CD8_I    <- 1000          # cells            CD8 activation saturation (I) (balanced for low-antigen priming)
  pars$k_CD8_F    <- 3.0           # AU/day           CD8 priming by IFN-I (upscaled for delay chain)
  pars$K_CD8_F    <- 20            # AU               CD8 activation saturation (F_I) (balanced for low-antigen priming)
  pars$k_mat      <- 0.14          # /day             CD8_N → CD8_E maturation (delay ~7d)
  pars$d_CD8_N    <- 0.05          # /day             Primed CD8 decay (t1/2 ~14d)
  pars$d_CD8_E    <- 0.10          # /day             Effector CD8 decay (t1/2 ~7d)
  pars$alpha_CD4_help <- 2.0       # dimensionless    CD4 help max boost for CD8
  pars$K_CD4_help <- 5             # AU               CD4 help half-max for CD8

  # Steepness of the effector-to-memory decay transition at the 1 AU threshold
  # (applies to CD8_E and IgG). Numerical regularisation of what were if/else
  # steps: those steps are attracting discontinuities that pinned the solver at
  # the crossing. See the comment at dCD8_E in model/hantavirus_qsp.R.
  pars$n_mem_switch <- 50          # dimensionless    Memory-switch Hill exponent

  # CD8+ mediated infected cell clearance
  pars$delta_CD8  <- 1e-3          # /(AU*day)        CD8-mediated infected cell clearance (calibrated for acute + clearance)

  # CD8 bystander endothelial damage (IL-15/NKG2D axis)
  pars$k_P_CD8    <- 0.005         # /day             Permeability drive by CD8 bystander

  # IgM antibodies (early B cell response)
  pars$k_IgM      <- 0.3           # AU/day           IgM production rate
  pars$d_IgM      <- 0.03          # /day             IgM decay (t1/2 ~ 23d)
  pars$K_CD4_B    <- 5             # AU               CD4 help half-max for B cells

  # IgG antibodies (class-switched, neutralizing)
  pars$k_switch   <- 0.05          # /day             IgM→IgG class switch rate
  pars$d_IgG      <- 0.033         # /day             IgG decay (t1/2 ~ 21d)
  # Antigen-driven IgG: direct class-switch from infected cell + CD4 help
  # Complements IgM-dependent pathway — captures germinal center kinetics
  # where antigen persistence drives ongoing class-switch recombination
  pars$k_IgG_direct <- 0.02        # /day             Antigen-driven IgG production
  pars$K_IgG       <- 500          # cells            Antigen half-max for IgG

  # Antibody-mediated viral neutralization
  pars$k_neut_IgM <- 0.005         # /(AU*day)        IgM neutralization rate (increased)
  pars$k_neut_IgG <- 0.05          # /(AU*day)        IgG neutralization rate (10x increase for clearance)

  # ===========================================================================
  # 3. ENDOTHELIAL PERMEABILITY & PLATELETS
  # ===========================================================================
  pars$k_PV      <- 0.5          # /day             # Permeability drive by virus (norm to V_ref)
  pars$k_PC      <- 0.1          # /day             # Permeability drive by cytokines
  pars$d_P       <- 0.5          # /day             # Permeability recovery rate
  # Platelet production balances natural loss at the clinical baseline:
  # disease-free steady state = k_PLT_prod / k_PLT_loss = 25000 / 0.1 = 250,000/uL.
  # (Prior value 3500 collapsed the disease-free baseline to 35,000/uL and drove
  # the disease nadir to ~1,100/uL — a calibration bug fixed 2026-06-03.)
  pars$k_PLT_prod  <- 25000      # platelets/uL/day # Platelet production rate (balances baseline)
  pars$k_PLT_loss  <- 0.1        # /day             # Platelet natural loss
  # Consumption reduced 0.01 -> 0.002 so the placebo nadir lands in the
  # clinically realistic severe-thrombocytopenia range (~40,000-50,000/uL)
  # rather than the implausible ~1,100/uL produced by the old balance.
  pars$k_PLT_cons  <- 0.002      # /day             # Platelet consumption per permeability
  pars$PLT_0     <- 250000       # platelets/uL     # Baseline platelet count
  pars$P_thresh  <- 0.5          # dimensionless    # Permeability threshold for fluid leak

  # ===========================================================================
  # 4. ORGAN INJURY MODULES
  # ===========================================================================
  pars$k_KP      <- 0.15         # /day             # Renal injury per permeability unit
  pars$k_KC      <- 0.05         # /day             # Renal injury per cytokine unit (norm)
  pars$k_KH      <- 0.2          # /day             # Renal injury from thrombocytopenia
  pars$d_K       <- 0.8          # /day             # Renal recovery rate (increased)
  pars$k_LP      <- 0.2          # /day             # Lung injury per permeability unit
  pars$k_LC      <- 0.075        # /day             # Lung injury per cytokine unit (norm)
  pars$k_Lfluid  <- 0.2          # /day             # Lung injury from fluid leak
  pars$d_L       <- 0.6          # /day             # Lung recovery rate (increased)

  # ===========================================================================
  # 5. RIBAVIRIN PK
  # ===========================================================================
  # CRITICAL RECALIBRATION (2026-05-12):
  # Ribavirin is essentially a prodrug with massive RBC partitioning:
  #   - RBC:plasma ratio = 60:1
  #   - 87% of intracellular ribavirin is phosphorylated (active RTP)
  #   - Chronic dosing half-life = 151 h (~6.3 days) due to RBC saturation
  #   - Apparent Vd = 800-5,000 L (not 45 L) due to RBC sequestration
  #
  # The 1-compartment model uses Vd = 45 L (plasma volume) but the
  # CLEARANCE is recalibrated to match the EFFECTIVE half-life observed
  # during chronic dosing, which captures the slow RBC uptake/release
  # without adding an explicit RBC compartment.
  #
  # k_el = ln(2) / 6.3 days = 0.11 /day
  # CL = k_el * Vd = 0.11 * 45 = 5.0 L/day (NOT 288 L/day)
  #
  # EC50_RBV = 1-2 ug/mL is retained as an "effective" value that accounts
  # for the 60:1 RBC:plasma partitioning - plasma concentrations underestimate
  # intracellular RTP concentrations by ~60x.
  #
  # NOT the plasma volume of distribution (800-5000 L). This defines a
  # model-internal EFFECT-SITE exposure scale: ribavirin is transported into
  # cells and phosphorylated, so intracellular exposure greatly exceeds plasma.
  # EC50_RBV is calibrated against this same scale, so the antiviral effect is
  # invariant to it; only the concentration axis of the PK figure is affected.
  pars$Vd_RBV    <- 45           # L                # Effect-site volume (internal scale)
  pars$CL_RBV_std <- 5.0        # L/day            # Effective clearance at eGFR=90 (t1/2 = 6.3 days)
  pars$CL_RBV    <- 5.0         # L/day            # Default clearance (eGFR=90)
  pars$eGFR_ref  <- 90           # mL/min           # Reference eGFR for CL scaling
  pars$k_hgb     <- 0.0009       # (g/dL)/(ug/mL)/day  # Hgb drop rate constant (bounded toxicity model)
  pars$d_hgb     <- 0.05         # /day             # Hgb recovery rate
  pars$Hgb_drop_max <- 4.0       # g/dL             # Maximum modeled RBV hemoglobin decline

  # ===========================================================================
  # 6. FAVIPIRAVIR PK
  # ===========================================================================
  pars$ka_FAV    <- 1.5 * 24     # /day             # Absorption rate (1.5/h * 24)
  pars$Vd_FAV    <- 30           # L                # Volume of distribution
  # Favipiravir plasma clearance is linear (first-order). A second,
  # concentration-dependent term (CL_FAV_nl = 96 L/day, Km_FAV = 50 ug/mL) was
  # removed: it was documented as saturable metabolism but made total clearance
  # RISE with concentration, whereas favipiravir inhibits aldehyde oxidase and
  # so its clearance falls with exposure. See the comment at dC_FAV in
  # model/hantavirus_qsp.R.
  pars$CL_FAV    <- 8 * 24       # L/day            # Linear clearance (8 L/h * 24)
  pars$k_form    <- 0.1 * 24     # /day             # RTP formation rate (0.1/h * 24)
  pars$k_elim_RTP <- 0.05 * 24  # /day             # RTP elimination rate (0.05/h * 24)

  # ===========================================================================
  # 7. PHARMACODYNAMICS (Emax models)
  # ===========================================================================
  # Ribavirin PD — recalibrated (2026-05-12) to match Huggins trial:
  #   - Huggins et al. (J Infect Dis 1991, PMID:1683355): 7-fold mortality RRR
  #     (placebo ~7-10% vs ribavirin ~1-2%, P=0.01)
  #   - Literature EC50 for HTNV: ~40-50 μM (~10-12 μg/mL)
  #   - With RBC partitioning (60:1), effective plasma EC50 is higher than
  #     originally assumed. Calibrated to reproduce clinical RRR.
  #   - Emax reduced to 0.70 to prevent complete viral suppression,
  #     consistent with partial clinical efficacy (not sterilizing cure)
  pars$Emax_RBV  <- 0.70         # dimensionless    # Max ribavirin effect (reduced for partial efficacy)
  pars$EC50_RBV  <- 8            # ug/mL            # Ribavirin half-max (calibrated for Huggins 7-fold RRR)
  pars$gamma_RBV <- 1.5          # dimensionless    # Hill coefficient RBV
  pars$Emax_FAV  <- 0.98         # dimensionless    # Max favipiravir effect
  pars$EC50_FAV  <- 0.5          # ug/mL            # Favipiravir RTP half-max
  pars$gamma_FAV <- 1.2          # dimensionless    # Hill coefficient FAV
  pars$psi       <- 0.1          # dimensionless    # Synergy interaction coefficient
  pars$rbv_immuno_scale <- 1.0   # dimensionless    # Scale for RBV immunomodulatory effects
  pars$rbv_endothelial_scale <- 1.0 # dimensionless # Scale for RBV endothelial protection
  pars$rbv_recovery_scale <- 1.0 # dimensionless    # Scale for RBV organ recovery enhancement

  # ===========================================================================
  # 8. CLINICAL ENDPOINT THRESHOLDS — syndrome-specific half-saturation model
  # ===========================================================================
  # Half-saturation constants for mortality risk components
  pars$K50       <- 30           # K value for 50% renal mortality risk contribution
  pars$L50       <- 100          # L value for 50% lung mortality risk contribution
  pars$C50       <- 300          # C_pro value for 50% cytokine mortality risk contribution

  # Dialysis / ECMO half-saturation constants — DECOUPLED from the mortality
  # K50/L50 (added 2026-06-03). Previously dialysis/ECMO reused K50=30 and
  # L50=100, which forced placebo rates of ~64% for both. Separate, higher
  # thresholds give clinically plausible and differentiated placebo rates
  # (dialysis ~40%, ECMO ~28-30%) without perturbing the mortality calibration.
  pars$K50_dialysis <- 70        # K value for 50% dialysis (renal replacement) risk
  pars$L50_ecmo     <- 400       # L value for 50% ECMO (cardiopulmonary support) risk

  # HFRS mortality weights (renal-dominant; target ~10% placebo mortality)
  pars$w_K_HFRS  <- 0.12         # Renal weight for HFRS
  pars$w_L_HFRS  <- 0.02         # Lung weight for HFRS (minor)
  pars$w_C_HFRS  <- 0.05         # Cytokine weight for HFRS

  # HCPS mortality weights removed: the lung-weighted mapping was intended to
  # reproduce ~40% placebo mortality and delivered 26%, so it was never a
  # calibrated arm of this model. This analysis is HFRS-only.

  # AUC-based mortality components — captures cumulative organ injury burden
  # (not just peak severity). This makes treatment-accelerated recovery
  # translate to mortality benefit. AUC units are day·AU. Weights are
  # normalized so that peak + AUC contributions sum appropriately.
  pars$w_auc      <- 0.30         # Weight for AUC component (peak gets 1-w_auc)
  pars$K_auc50    <- 200          # day·AU          K AUC for 50% renal risk
  pars$L_auc50    <- 500          # day·AU          L AUC for 50% lung risk
  pars$C_auc50    <- 1500         # day·AU          C_pro AUC for 50% cytokine risk

  pars
}

#' Build a parameter table for export
#'
#' @param pars Named parameter list (from [get_parameters()])
#' @return Data frame with columns: parameter, value, units, description,
#'   source, confidence
#' @export
build_parameter_table <- function(pars = get_parameters()) {
  tbl <- data.frame(
    parameter   = character(),
    value       = numeric(),
    units       = character(),
    description = character(),
    source      = character(),
    confidence  = character(),
    stringsAsFactors = FALSE
  )

  # --- Viral dynamics ---
  tbl <- rbind(tbl, data.frame(parameter="T_0", value=1e6, units="cells",
    description="Initial susceptible target cells",
    source="Assumed (typical epithelial pool)", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="I_0", value=0, units="cells",
    description="Initial infected cells", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="V_0", value=1000, units="copies/mL",
    description="Initial viral inoculum (increased for establishment)", source="Estimated from clinical viraemia", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="beta", value=1.0e-7, units="mL/copy/day",
    description="Infection rate constant (calibrated to ~9.6% placebo mortality; gives a viral-establishment index beta*T_0*p/c = 3.33 /day, which is NOT a reproduction number - see the note at pars$beta)", source="Calibrated to Huggins JID 1991 placebo mortality", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="p", value=100, units="copies/cell/day",
    description="Viral production rate", source="Fitted to peak viraemia", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="c", value=3, units="/day",
    description="Viral clearance rate", source="Literature (RNA virus half-life ~5.5h)", confidence="high"))
  tbl <- rbind(tbl, data.frame(parameter="k_T_reg", value=0.30, units="/day",
    description="Target-cell regeneration rate (logistic, carrying capacity T_0)", source="Assumed (epithelial turnover)", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="symptom_onset_day", value=5, units="days",
    description="Infection-to-symptom-onset burn-in; therapeutic simulations start at this point", source="Assumed (hantavirus incubation 7-42 d)", confidence="low"))

  # --- NSs protein ---
  tbl <- rbind(tbl, data.frame(parameter="k_NSs", value=0.5, units="/day",
    description="NSs production rate", source="Assumed (viral protein kinetics)", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="d_NSs", value=1.5, units="/day",
    description="NSs decay rate (t1/2 ~ 11h)", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="K_NSs_sup", value=10, units="AU",
    description="NSs suppression half-max", source="Assumed", confidence="low"))

  # --- Type I IFN ---
  tbl <- rbind(tbl, data.frame(parameter="k_FI", value=5.0, units="AU/day",
    description="IFN-I max production rate", source="Assumed (innate response)", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="K_FI", value=1e5, units="cells",
    description="IFN production saturation", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="d_F_I", value=3.5, units="/day",
    description="IFN-I decay (t1/2 ~ 4.8h)", source="Literature (IFN-alpha PK)", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="E_max_IFN", value=0.8, units="dimensionless",
    description="Max IFN antiviral effect", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="EC50_F", value=50, units="AU",
    description="IFN half-max antiviral", source="Assumed", confidence="low"))

  # --- NK cells ---
  tbl <- rbind(tbl, data.frame(parameter="k_NK_I", value=0.3, units="AU/day",
    description="NK activation by infected cells", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="K_NK_I", value=1e5, units="cells",
    description="NK activation saturation (I)", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="k_NK_F", value=0.2, units="AU/day",
    description="NK activation by IFN-I", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="K_NK_F", value=50, units="AU",
    description="NK activation saturation (F_I)", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="d_NK", value=0.7, units="/day",
    description="NK decay (t1/2 ~ 24h)", source="Literature (NK cell kinetics)", confidence="medium"))

  # --- Type II IFN ---
  tbl <- rbind(tbl, data.frame(parameter="k_FII_NK", value=0.2, units="/day",
    description="IFN-gamma production by NK", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="d_F_II", value=2.8, units="/day",
    description="IFN-II decay (t1/2 ~ 6h)", source="Literature (IFN-gamma PK)", confidence="medium"))

  # --- Cytokines ---
  tbl <- rbind(tbl, data.frame(parameter="k_CI", value=2.0, units="AU/day",
    description="Cytokine prod by infected cells", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="k_CF2", value=1.5, units="AU/day",
    description="Cytokine amplification by IFN-II", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="k_CC", value=0.8, units="/day",
    description="Cytokine auto-amplification max (calibrated)", source="Calibrated to mortality targets", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="K_CC", value=200, units="AU",
    description="Auto-amplification half-max (calibrated)", source="Calibrated to mortality targets", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="d_C", value=5.5, units="/day",
    description="Cytokine decay (t1/2 ~ 3h)", source="Literature (IL-6/TNF half-life)", confidence="medium"))

  # --- Anti-inflammatory ---
  tbl <- rbind(tbl, data.frame(parameter="k_anti", value=0.5, units="AU/day",
    description="IL-10 production rate", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="K_anti", value=100, units="AU",
    description="IL-10 production saturation", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="d_anti", value=0.5, units="/day",
    description="IL-10 decay (t1/2 ~ 33h)", source="Literature (IL-10 half-life)", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="k_suppress", value=1.0, units="/day",
    description="Max IL-10 suppression rate", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="K_suppress", value=30, units="AU",
    description="IL-10 suppression potency", source="Assumed", confidence="low"))

  # --- Reference values ---
  tbl <- rbind(tbl, data.frame(parameter="I_ref", value=1e6, units="cells",
    description="Reference infected cells", source="Defined scale", confidence="high"))
  tbl <- rbind(tbl, data.frame(parameter="F_I_ref", value=100, units="AU",
    description="Reference Type I IFN", source="Defined scale", confidence="high"))
  tbl <- rbind(tbl, data.frame(parameter="F_II_ref", value=50, units="AU",
    description="Reference Type II IFN", source="Defined scale", confidence="high"))
  tbl <- rbind(tbl, data.frame(parameter="C_ref", value=500, units="AU",
    description="Reference cytokines", source="Defined scale", confidence="high"))
  tbl <- rbind(tbl, data.frame(parameter="V_ref", value=1e5, units="copies/mL",
    description="Reference viral load", source="Defined scale", confidence="high"))

  # --- Infected cell clearance ---
  tbl <- rbind(tbl, data.frame(parameter="delta_natural", value=0.03, units="/day",
    description="Natural infected cell turnover (t1/2 ~ 23d)",
    source="Literature (non-cytopathic hantavirus)", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="delta_NK", value=5e-6, units="1/(AU*day)",
    description="NK-mediated infected cell clearance", source="Assumed", confidence="low"))

  # --- Auto-amplification gate ---
  tbl <- rbind(tbl, data.frame(parameter="I_gate", value=10000, units="cells",
    description="Auto-amplification gate threshold (raised to match clinical RBV efficacy)", source="Calibrated", confidence="medium"))

  # --- Adaptive immunity: CD4+ helper T cells ---
  tbl <- rbind(tbl, data.frame(parameter="k_CD4_I", value=1.5, units="AU/day",
    description="CD4 activation by infected cells", source="Assumed (T cell kinetics, PMID:19072554)", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="K_CD4_I", value=500, units="cells",
    description="CD4 activation saturation (I), lowered to match CD8 priming sensitivity", source="Calibrated (CD4/CD8 consistency)", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="k_CD4_F", value=1.0, units="AU/day",
    description="CD4 activation by IFN-I", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="K_CD4_F", value=50, units="AU",
    description="CD4 activation saturation (F_I)", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="d_CD4", value=0.03, units="/day",
    description="CD4 decay (t1/2 ~ 7d, effector)", source="General viral immunology", confidence="medium"))

  # --- Adaptive immunity: CD8+ cytotoxic T cells (delayed maturation chain) ---
  tbl <- rbind(tbl, data.frame(parameter="k_CD8_I", value=5.0, units="AU/day",
    description="CD8 priming by infected cells (upscaled for delay chain)", source="Calibrated (T cell kinetics, PMID:19072554)", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="K_CD8_I", value=1000, units="cells",
    description="CD8 activation saturation (I), lowered for low-antigen priming", source="Calibrated (prevents viral rebound)", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="k_CD8_F", value=3.0, units="AU/day",
    description="CD8 priming by IFN-I (upscaled for delay chain)", source="Calibrated", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="K_CD8_F", value=20, units="AU",
    description="CD8 activation saturation (F_I), lowered for low-antigen priming", source="Calibrated (prevents viral rebound)", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="k_mat", value=0.14, units="/day",
    description="CD8_N to CD8_E maturation rate (delay ~7d)", source="PMID:21525363 (T cell expansion kinetics)", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="d_CD8_N", value=0.05, units="/day",
    description="Primed CD8 decay (t1/2 ~14d)", source="General viral immunology", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="d_CD8_E", value=0.10, units="/day",
    description="Effector CD8 decay (t1/2 ~7d)", source="General viral immunology", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="alpha_CD4_help", value=2.0, units="dimensionless",
    description="CD4 help max boost for CD8", source="Assumed (T cell help)", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="K_CD4_help", value=5, units="AU",
    description="CD4 help half-max for CD8", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="n_mem_switch", value=50, units="dimensionless",
    description="Hill exponent of the effector-to-memory decay transition at the 1 AU threshold (CD8_E, IgG); numerical regularisation of a step function, not a fitted biological quantity",
    source="Numerical (agrees with the step function it replaces to within 1% of the decay rate outside +/-9% of the threshold)", confidence="structural"))
  tbl <- rbind(tbl, data.frame(parameter="delta_CD8", value=1e-3, units="1/(AU*day)",
    description="CD8-mediated infected cell clearance", source="Calibrated (CD8 killing)", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="k_P_CD8", value=0.005, units="/day",
    description="Permeability drive by CD8_E bystander (IL-15/NKG2D)", source="PMID:36590594", confidence="low"))

  # --- Adaptive immunity: IgM antibodies ---
  tbl <- rbind(tbl, data.frame(parameter="k_IgM", value=0.3, units="AU/day",
    description="IgM production rate", source="Assumed (B cell kinetics, PMID:2902106)", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="d_IgM", value=0.03, units="/day",
    description="IgM decay (t1/2 ~ 23d)", source="PUUV IgM kinetics (PMID:2902106)", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="K_CD4_B", value=5, units="AU",
    description="CD4 help half-max for B cells", source="Assumed", confidence="low"))

  # --- Adaptive immunity: IgG antibodies ---
  tbl <- rbind(tbl, data.frame(parameter="k_switch", value=0.05, units="/day",
    description="IgM to IgG class switch rate", source="Assumed (class switch)", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="k_IgG_direct", value=0.02, units="/day",
    description="Antigen-driven (I) IgG production rate", source="Assumed (germinal centre kinetics)", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="K_IgG", value=500, units="cells",
    description="Antigen half-max for direct IgG production", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="d_IgG", value=0.033, units="/day",
    description="IgG decay (t1/2 ~ 21d)", source="General antibody kinetics", confidence="medium"))

  # --- Adaptive immunity: antibody neutralization ---
  tbl <- rbind(tbl, data.frame(parameter="k_neut_IgM", value=0.005, units="1/(AU*day)",
    description="IgM neutralisation rate", source="Calibrated (lower affinity)", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="k_neut_IgG", value=0.05, units="1/(AU*day)",
    description="IgG neutralisation rate (higher affinity)", source="Calibrated (PMID:30463919)", confidence="low"))

  # --- Permeability / platelets ---
  tbl <- rbind(tbl, data.frame(parameter="k_PV", value=0.5, units="/day",
    description="Permeability drive by virus (norm to V_ref)", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="k_PC", value=0.1, units="/day",
    description="Permeability drive by cytokines (norm)", source="Literature (cytokine storm)", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="d_P", value=0.5, units="/day",
    description="Permeability recovery rate", source="Assumed (endothelial repair)", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="k_PLT_prod", value=25000, units="platelets/uL/day",
    description="Platelet production rate (balances baseline: prod/loss = PLT_0)", source="Physiology (normal turnover ~7-10d)", confidence="high"))
  tbl <- rbind(tbl, data.frame(parameter="k_PLT_loss", value=0.1, units="/day",
    description="Platelet natural loss", source="Physiology", confidence="high"))
  tbl <- rbind(tbl, data.frame(parameter="k_PLT_cons", value=0.002, units="/day",
    description="Platelet consumption per permeability (calibrated to ~40-50k/uL nadir)", source="Calibrated (consumption coagulopathy)", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="PLT_0", value=250000, units="platelets/uL",
    description="Baseline platelet count", source="Clinical reference range", confidence="high"))
  tbl <- rbind(tbl, data.frame(parameter="P_thresh", value=0.5, units="dimensionless",
    description="Permeability threshold for fluid leak", source="Assumed", confidence="low"))

  # --- Organ injury ---
  tbl <- rbind(tbl, data.frame(parameter="k_KP", value=0.15, units="/day",
    description="Renal injury per permeability unit", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="k_KC", value=0.05, units="/day",
    description="Renal injury per cytokine unit (norm)", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="k_KH", value=0.2, units="/day",
    description="Renal injury from thrombocytopenia", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="d_K", value=0.5, units="/day",
    description="Renal recovery rate", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="k_LP", value=0.2, units="/day",
    description="Lung injury per permeability unit", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="k_LC", value=0.075, units="/day",
    description="Lung injury per cytokine unit (norm)", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="k_Lfluid", value=0.2, units="/day",
    description="Lung injury from fluid leak", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="d_L", value=0.4, units="/day",
    description="Lung recovery rate", source="Assumed", confidence="low"))

  # --- Ribavirin PK ---
  tbl <- rbind(tbl, data.frame(parameter="Vd_RBV", value=45, units="L",
    description="Ribavirin effect-site volume, defining a model-internal exposure scale. NOT the plasma volume of distribution, which is 800-5000 L; ribavirin is transported into cells and phosphorylated, so the effect-site concentration is far above plasma. EC50_RBV is calibrated on this same scale, so the antiviral effect is unaffected by it", source="Model-internal exposure scale", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="CL_RBV_std", value=5.0, units="L/day",
    description="Ribavirin effective clearance at eGFR=90 (t1/2=6.2d, RBC-adjusted)", source="Recalibrated from chronic dosing PK (t1/2 ~150 h)", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="CL_RBV", value=5.0, units="L/day",
    description="Ribavirin clearance actually used; rescaled per patient by eGFR", source="Set to CL_RBV_std at the reference eGFR", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="eGFR_ref", value=90, units="mL/min",
    description="Reference eGFR at which CL_RBV_std applies (linear renal scaling)", source="Clinical reference", confidence="high"))
  tbl <- rbind(tbl, data.frame(parameter="k_hgb", value=0.0009, units="(g/dL)/(ug/mL)/day",
    description="Hgb drop rate per RBV concentration in bounded toxicity submodel", source="Calibrated to clinical RBV-associated decline (~2-4 g/dL)", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="d_hgb", value=0.05, units="/day",
    description="Hgb recovery rate", source="Physiology (RBC lifespan ~120d)", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="Hgb_drop_max", value=4.0, units="g/dL",
    description="Upper bound on modelled ribavirin-associated haemoglobin decline", source="Clinical constraint for 7-10 day IV RBV courses", confidence="medium"))

  # --- Favipiravir PK ---
  tbl <- rbind(tbl, data.frame(parameter="ka_FAV", value=36, units="/day",
    description="Favipiravir absorption rate (1.5/h)", source="Madelain et al. 2016 (Clin Pharmacokinet)", confidence="high"))
  tbl <- rbind(tbl, data.frame(parameter="Vd_FAV", value=30, units="L",
    description="Favipiravir volume of distribution", source="FDA review", confidence="high"))
  tbl <- rbind(tbl, data.frame(parameter="CL_FAV", value=192, units="L/day",
    description="Favipiravir linear clearance (8 L/h)", source="Literature", confidence="high"))
  tbl <- rbind(tbl, data.frame(parameter="k_form", value=2.4, units="/day",
    description="RTP formation rate (0.1/h)", source="Assumed (intracellular activation)", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="k_elim_RTP", value=1.2, units="/day",
    description="RTP elimination rate (0.05/h)", source="Assumed (intracellular half-life)", confidence="low"))

  # --- PD ---
  tbl <- rbind(tbl, data.frame(parameter="Emax_RBV", value=0.70, units="dimensionless",
    description="Max ribavirin effect (reduced for partial efficacy)", source="Calibrated for Huggins 7-fold RRR (JID 1991)", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="EC50_RBV", value=8, units="ug/mL",
    description="Ribavirin half-max effect concentration on the Vd_RBV effect-site scale (calibrated to the day-1 window of the Huggins trial)", source="Calibrated to Huggins JID 1991", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="gamma_RBV", value=1.5, units="dimensionless",
    description="Hill coefficient for RBV", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="Emax_FAV", value=0.98, units="dimensionless",
    description="Max favipiravir antiviral effect", source="In vitro (broad-spectrum)", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="EC50_FAV", value=0.5, units="ug/mL",
    description="Favipiravir RTP half-max concentration", source="In vitro IC50 for RTP", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="gamma_FAV", value=1.2, units="dimensionless",
    description="Hill coefficient for FAV", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="psi", value=0.1, units="dimensionless",
    description="Synergy interaction coefficient", source="Assumed (Bliss + synergy)", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="rbv_immuno_scale", value=1.0, units="dimensionless",
    description="Scale factor for ribavirin immunomodulatory effects", source="Ablation-control parameter", confidence="defined"))
  tbl <- rbind(tbl, data.frame(parameter="rbv_endothelial_scale", value=1.0, units="dimensionless",
    description="Scale factor for ribavirin direct endothelial protection", source="Ablation-control parameter", confidence="defined"))
  tbl <- rbind(tbl, data.frame(parameter="rbv_recovery_scale", value=1.0, units="dimensionless",
    description="Scale factor for ribavirin organ-injury recovery enhancement", source="Ablation-control parameter", confidence="defined"))

  # --- Syndrome-specific mortality half-saturation model ---
  tbl <- rbind(tbl, data.frame(parameter="K50", value=30, units="AU",
    description="K value for 50% renal mortality risk", source="Calibrated", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="L50", value=100, units="AU",
    description="L value for 50% lung mortality risk", source="Calibrated", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="C50", value=300, units="AU",
    description="C_pro value for 50% cytokine mortality risk", source="Calibrated", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="K50_dialysis", value=70, units="AU",
    description="K value for 50% dialysis risk (decoupled from mortality K50)", source="Calibrated (~40% placebo dialysis)", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="L50_ecmo", value=400, units="AU",
    description="L value for 50% ECMO risk (decoupled from mortality L50)", source="Calibrated (~28-30% placebo ECMO)", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="w_K_HFRS", value=0.12, units="dimensionless",
    description="Renal mortality weight for HFRS", source="Calibrated for ~10% HFRS placebo", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="w_L_HFRS", value=0.02, units="dimensionless",
    description="Lung mortality weight for HFRS", source="Calibrated", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="w_C_HFRS", value=0.05, units="dimensionless",
    description="Cytokine mortality weight for HFRS", source="Calibrated", confidence="medium"))
  tbl <- rbind(tbl, data.frame(parameter="w_auc", value=0.30, units="dimensionless",
    description="Weight for AUC-based mortality component", source="Assumed (cumulative organ injury)", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="K_auc50", value=200, units="day*AU",
    description="K AUC for 50% renal risk", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="L_auc50", value=500, units="day*AU",
    description="L AUC for 50% lung risk", source="Assumed", confidence="low"))
  tbl <- rbind(tbl, data.frame(parameter="C_auc50", value=1500, units="day*AU",
    description="C_pro AUC for 50% cytokine risk", source="Assumed", confidence="low"))

  # get_parameters() is the SINGLE SOURCE OF TRUTH for parameter values. The
  # value= literals in the rows above are legacy placeholders only — they are
  # NOT authoritative. We discard them here and re-populate every value directly
  # from get_parameters(), so the exported table can never drift from the model.
  tbl$value <- NA_real_
  value_rows <- tbl$parameter %in% names(pars)
  tbl$value[value_rows] <- vapply(
    tbl$parameter[value_rows],
    function(nm) as.numeric(pars[[nm]]),
    numeric(1)
  )

  missing_rows <- setdiff(names(pars), tbl$parameter)
  if (length(missing_rows) > 0) {
    tbl <- rbind(tbl, data.frame(
      parameter = missing_rows,
      value = vapply(missing_rows, function(nm) as.numeric(pars[[nm]]), numeric(1)),
      units = "see code",
      description = "Parameter defined in get_parameters(); metadata pending",
      source = "Defined in model/parameters.R",
      confidence = "low",
      stringsAsFactors = FALSE
    ))
  }

  tbl
}
