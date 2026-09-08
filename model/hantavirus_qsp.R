#' Hantavirus QSP Full ODE Model
#'
#' Coupled viral dynamics, innate and adaptive immune, endothelial
#' permeability, organ injury, and antiviral PK/PD model.
#'
#' State vector (23 compartments):
#'   T   - Susceptible target cells
#'   I   - Infected cells
#'   V   - Viral RNA (copies/mL)
#'   NSs - Hantavirus NSs protein (IFN antagonist)
#'   F_I - Type I IFN (IFN-alpha/beta, antiviral)
#'   NK  - NK cell activation
#'   F_II - Type II IFN (IFN-gamma, pro-inflammatory)
#'   C_pro - Pro-inflammatory cytokine burden
#'   C_anti - Anti-inflammatory cytokine burden (IL-10-like)
#'   P   - Vascular permeability index
#'   PLT - Platelet count (/uL)
#'   K   - Renal injury index
#'   L   - Lung/cardiopulmonary injury index
#'   CD8_N - CD8+ naive/primed T cells (antigen-primed, not yet expanded)
#'   CD8_E - CD8+ effector T cells (expanded, cytotoxic)
#'   CD4 - CD4+ helper T cells (adaptive)
#'   IgM - IgM antibodies (early adaptive)
#'   IgG - Neutralizing IgG antibodies (late adaptive)
#'   C_RBV      - Ribavirin plasma concentration (ug/mL)
#'   C_FAV_gut  - Favipiravir gut compartment (mg)
#'   C_FAV      - Favipiravir plasma concentration (ug/mL)
#'   C_FAVI_RTP - Favipiravir active metabolite (ug/mL)
#'   Hgb_drop   - Hemoglobin drop (g/dL)
#'
#' @param t    Time (days)
#' @param y    State vector
#' @param pars Named parameter list
#' @return List with dydt
#' @export
hantavirus_qsp_ode <- function(t, y, pars) {
  with(as.list(c(y, pars)), {

    # --- Antiviral effect calculations ---
    # Ribavirin effect (handle C_RBV = 0 gracefully)
    E_RBV <- if (C_RBV > 0) {
      (Emax_RBV * C_RBV^gamma_RBV) / (EC50_RBV^gamma_RBV + C_RBV^gamma_RBV)
    } else {
      0
    }

    # Favipiravir effect via RTP metabolite
    E_FAV <- if (C_FAVI_RTP > 0) {
      (Emax_FAV * C_FAVI_RTP^gamma_FAV) / (EC50_FAV^gamma_FAV + C_FAVI_RTP^gamma_FAV)
    } else {
      0
    }

    # Synergy interaction term
    psi_combo <- psi * E_RBV * E_FAV

    # Ribavirin immunomodulatory effect (reduces cytokine storm independently
    # of antiviral activity — consistent with RBV effects on T-cell response)
    rbv_immuno_factor <- if ("rbv_immuno_scale" %in% names(pars)) rbv_immuno_scale else 1.0
    E_RBV_immuno <- if (C_RBV > 0) {
      rbv_immuno_factor * 0.9 * C_RBV / (2 + C_RBV)
    } else {
      0
    }

    # --- Viral Dynamics ---
    # Target cell regeneration: logistic growth with carrying capacity T_0
    # and intrinsic growth rate k_T_reg. The logistic term k_T_reg * T * (1 - T/T_0)
    # goes to zero as T → 0, making T = 0 an absorbing state that cannot recover.
    # To prevent this numerical artifact, we add a small seed term k_T_reg * eps
    # that represents minimal epithelial progenitor activity (stem-cell-driven
    # turnover at ~0.0003 cells/day when T ≈ 0). This is mathematically equivalent
    # to logistic growth from a tiny seed population and becomes negligible once
    # T grows above eps. When T >> eps (nearly always), the infection term
    # is the standard beta * V * T.
    eps <- 1e-3
    dT  <- k_T_reg * (T + eps) * (1 - T / T_0) - beta * V * T

    # Infected cell clearance: non-cytopathic with NK-mediated clearance
    # Ribavirin enhances NK cell-mediated clearance (immunomodulatory effect)
    NK_clearance_boost <- 1 + 0.8 * E_RBV_immuno

    # CD8+ T cell clearance (uses effector compartment CD8_E): two regimes
    # 1. High I (bulk infection): mass-action killing (delta_CD8 * CD8_E * I)
    # 2. Low I (residual infection): CTLs scan tissue and kill individual
    #    infected cells at a per-cell rate proportional to CD8_E density.
    #    This regime is what clears the last infected cells when mass-action
    #    killing becomes ineffective.
    #    CAUTION: both branches are proportional to CD8_E * I, so this is NOT
    #    a switch between two functional forms - it is a single rate constant
    #    that rises from delta_CD8 (high I) to 5.0 (low I) as I falls.
    #    The per-cell clearance rate is 5.0 * CD8_E per day, so at CD8_E = 50
    #    it is 250/day, not 5.0/day as an earlier version of this comment said.
    CD8_mass_action <- delta_CD8 * CD8_E * I
    CD8_first_order_rate <- 5.0 * CD8_E  # /day per infected cell (very strong for sterilizing clearance)
    # Smooth transition: mass-action dominates when I > 200 cells,
    # first-order dominates when I < 200 cells
    w_mass <- I / (200 + I)  # weight for mass-action
    w_first <- 200 / (200 + I)  # weight for first-order
    CD8_clearance <- w_mass * CD8_mass_action + w_first * CD8_first_order_rate * I

    dI  <- beta * V * T - delta_natural * I -
           delta_NK * NK * I * NK_clearance_boost -
           CD8_clearance

    # IFN antiviral effect on viral production (E_max model)
    IFN_antiviral <- 1 - (E_max_IFN * F_I) / (EC50_F + F_I)

    dV  <- p * I * IFN_antiviral *
             (1 - E_RBV) * (1 - E_FAV) * (1 - psi_combo) -
           c * V - k_neut_IgM * IgM * V - k_neut_IgG * IgG * V
    # Antibody neutralization: stoichiometric binding (IgM + IgG).
    # No catalytic enhancement — the prior (1 + IgG/50) term was removed
    # because it caused over-clearance at high IgG levels seen in placebo.

    # --- NSs protein (hantavirus IFN antagonist) ---
    dNSs <- k_NSs * I - d_NSs * NSs

    # --- Type I IFN (antiviral, suppressed by NSs) ---
    # Ribavirin enhances IFN production (immunomodulatory effect)
    # consistent with RBV shifting T-cell response toward Th1
    IFN_boost <- 1 + 1.5 * E_RBV_immuno

    # IFN production saturated by infected cells, suppressed by NSs
    F_I_prod <- if (I > 0) {
      (k_FI * I / (K_FI + I)) * (1 - NSs / (K_NSs_sup + NSs)) * IFN_boost
    } else {
      0
    }
    dF_I <- F_I_prod - d_F_I * F_I

    # --- NK cell activation (by infected cells and Type I IFN) ---
    # Ribavirin enhances NK cell activation (immunomodulatory effect)
    NK_boost <- 1 + 1.2 * E_RBV_immuno

    NK_act_I <- if (I > 0) k_NK_I * I / (K_NK_I + I) * NK_boost else 0
    NK_act_F <- if (F_I > 0) k_NK_F * F_I / (K_NK_F + F_I) * NK_boost else 0
    dNK  <- NK_act_I + NK_act_F - d_NK * NK

    # --- Type II IFN (pro-inflammatory, from NK cells) ---
    dF_II <- k_FII_NK * NK - d_F_II * F_II

    # --- Pro-inflammatory cytokines (auto-amplification + IL-10 suppression) ---
    # Normalized production terms; RBV immunomodulatory effect reduces all three
    term_I   <- if (I > 0) k_CI * (I / I_ref) * (1 - E_RBV_immuno) else 0
    term_FII <- if (F_II > 0) k_CF2 * (F_II / F_II_ref) * (1 - E_RBV_immuno) else 0
    # Hill-2 auto-amplification (cytokine storm) -- GATED: only active when
    # infected cells exceed threshold, preventing runaway amplification after
    # viral suppression clears the infection trigger
    term_auto <- if (I > I_gate) k_CC * (C_pro^2 / (K_CC^2 + C_pro^2)) *
                   (1 - E_RBV_immuno) else 0
    # IL-10 anti-inflammatory suppression
    term_suppress <- if (C_anti > 0) {
      k_suppress * (C_anti / (K_suppress + C_anti)) * C_pro
    } else {
      0
    }
    dC_pro <- term_I + term_FII + term_auto - term_suppress - d_C * C_pro

    # --- Anti-inflammatory cytokines (delayed resolution, IL-10-like) ---
    C_anti_prod <- if (C_pro > 0) {
      k_anti * (C_pro / (K_anti + C_pro))
    } else {
      0
    }
    dC_anti <- C_anti_prod - d_anti * C_anti

    # --- Adaptive Immunity ---
    # CD4+ helper T cells: activated by infected cells + Type I IFN
    # Delayed onset (~day 5-7), peak ~day 14 (PMID: 19072554, 17641066)
    CD4_act_I <- if (I > 0) k_CD4_I * I / (K_CD4_I + I) else 0
    CD4_act_F <- if (F_I > 0) k_CD4_F * F_I / (K_CD4_F + F_I) else 0
    # NOTE: a CD4 memory floor was computed here and never used in dCD4. The
    # dead assignment has been removed; CD4 decays as a pure effector pool.
    dCD4 <- CD4_act_I + CD4_act_F - d_CD4 * CD4

    # CD8+ cytotoxic T cells: DELAYED activation via maturation chain
    # Antigen signal → CD8_N (primed pool) → CD8_E (effector pool)
    # This captures the 7-day biological lag: DC antigen presentation →
    # lymph node priming → clonal expansion → tissue migration
    # (PMID: 19072554, 18814258, 21525363)
    #
    # CD8_N: antigen-primed cells accumulating in lymph nodes
    # CD8_E: expanded effector cells that migrate to tissue and kill
    #
    # The maturation chain acts as a gamma-distributed delay (shape=2),
    # producing CD8_E peak at day 14-21 — consistent with literature
    # reporting PUUV-specific CD8+ cells encompassing up to 50% of
    # blood CD8+ T cells during weeks 2-3 (PMID: 21525363).
    CD4_help <- 1 + alpha_CD4_help * CD4 / (K_CD4_help + CD4)
    CD8_prime_I <- if (I > 0) k_CD8_I * I / (K_CD8_I + I) * CD4_help else 0
    CD8_prime_F <- if (F_I > 0) k_CD8_F * F_I / (K_CD8_F + F_I) * CD4_help else 0

    # Primed pool: accumulates antigen signal, transitions to effector
    dCD8_N <- CD8_prime_I + CD8_prime_F - k_mat * CD8_N - d_CD8_N * CD8_N

    # Effector pool: receives matured cells from CD8_N
    # Memory: reduced decay when CD8_E drops below 1 AU (t1/2 ~ months)
    CD8_E_decay_eff <- if (CD8_E < 1) {
      d_CD8_E * 0.1  # Memory half-life ~70 days vs 7 days for effectors
    } else {
      d_CD8_E
    }
    dCD8_E <- k_mat * CD8_N - CD8_E_decay_eff * CD8_E

    # IgM antibodies: early B cell response driven by antigen + CD4 help
    # Detectable day 1-3, peak ~day 30 (PMID: 2902106)
    CD4_help_B <- CD4 / (K_CD4_B + CD4)
    IgM_prod <- if (I > 0) k_IgM * (I / I_ref) * CD4_help_B else 0
    dIgM <- IgM_prod - d_IgM * IgM

    # IgG antibodies: class-switched from IgM (CD4-dependent) AND
    # directly from antigen + CD4 help (germinal center kinetics).
    # Onset ~day 10, sustained (PMID: 22787210, 2902106)
    # Memory B cells persist and rapidly reactivate (PMID: 18814258)
    IgG_switch <- if (IgM > 0) k_switch * IgM * CD4_help_B else 0
    IgG_direct <- if (I > 0) k_IgG_direct * I / (K_IgG + I) * CD4_help_B else 0
    # Memory B cell floor: reduced decay at low IgG levels
    IgG_decay_eff <- if (IgG < 1) {
      d_IgG * 0.1  # Memory half-life ~210 days vs 21 days for plasma cells
    } else {
      d_IgG
    }
    dIgG <- IgG_switch + IgG_direct - IgG_decay_eff * IgG

    # --- Endothelial Permeability & Platelets ---
    # Coupling uses normalized C_pro / C_ref; viral load normalized to V_ref
    # to prevent runaway permeability from high absolute viral loads
    # CD8_E bystander damage via IL-15/NKG2D axis (PMID: 36590594)
    C_norm <- C_pro / C_ref
    dP   <- k_PV * (V / V_ref) + k_PC * C_norm + k_P_CD8 * CD8_E - d_P * P
    dPLT <- k_PLT_prod - k_PLT_loss * PLT - k_PLT_cons * P * PLT

    # --- Organ Injury ---
    Hgb_threshold_factor <- if (PLT < PLT_0) (1 - PLT / PLT_0) else 0
    P_fluid_factor <- if (P > P_thresh) (P - P_thresh) else 0

    # Ribavirin has direct anti-inflammatory effects on endothelium
    # (consistent with in vitro data showing RBV reduces endothelial permeability)
    # AND enhances recovery from established organ injury (per Huggins trial data
    # showing reduced progression to oliguric phase)
    rbv_endothelial_factor <- if ("rbv_endothelial_scale" %in% names(pars)) rbv_endothelial_scale else 1.0
    rbv_recovery_factor <- if ("rbv_recovery_scale" %in% names(pars)) rbv_recovery_scale else 1.0
    E_RBV_endothelial <- if (C_RBV > 0) {
      rbv_endothelial_factor * 0.7 * C_RBV / (2 + C_RBV)
    } else {
      0
    }

    # Recovery enhancement: RBV speeds up organ injury resolution
    d_K_recovered <- d_K * (1 + rbv_recovery_factor * 0.5 * E_RBV_endothelial)
    d_L_recovered <- d_L * (1 + rbv_recovery_factor * 0.5 * E_RBV_endothelial)

    dK <- k_KP * P * (1 - E_RBV_endothelial) +
           k_KC * C_norm + k_KH * Hgb_threshold_factor - d_K_recovered * K
    dL <- k_LP * P * (1 - E_RBV_endothelial) +
           k_LC * C_norm + k_Lfluid * P_fluid_factor - d_L_recovered * L

    # --- Ribavirin PK (one-compartment IV, time in days) ---
    dC_RBV <- -(CL_RBV / Vd_RBV) * C_RBV

    # --- Favipiravir PK with active metabolite ---
    dC_FAV_gut <- -ka_FAV * C_FAV_gut
    dC_FAV <- (ka_FAV * C_FAV_gut / Vd_FAV) -
      (CL_FAV + CL_FAV_nl * C_FAV / (Km_FAV + C_FAV)) / Vd_FAV * C_FAV
    dC_FAVI_RTP <- k_form * C_FAV - k_elim_RTP * C_FAVI_RTP

    # --- Hemolytic Anemia (Ribavirin toxicity) ---
    Hgb_cap <- if ("Hgb_drop_max" %in% names(pars)) Hgb_drop_max else 4.0
    Hgb_drive <- k_hgb * C_RBV * max(0, 1 - Hgb_drop / Hgb_cap)
    dHgb_drop <- Hgb_drive - d_hgb * Hgb_drop

    list(c(dT, dI, dV, dNSs, dF_I, dNK, dF_II, dC_pro, dC_anti, dP, dPLT, dK, dL,
           dCD8_N, dCD8_E, dCD4, dIgM, dIgG,
           dC_RBV, dC_FAV_gut, dC_FAV, dC_FAVI_RTP, dHgb_drop))
  })
}

# NOTE: compute_clinical_endpoints() is defined once, in R/pd_models.R.
# It was previously duplicated here; the duplicate was removed so behaviour can
# no longer depend on the order in which the two files are sourced.
