#' Virtual Population Generation
#'
#' Generates a cohort of virtual patients with covariate distributions
#' reflecting the clinical epidemiology of Hantavirus infection.
#'
#' @name virtual_population
NULL

#' Truncated normal sampler
#'
#' Genuine truncation by inverse transform: a uniform draw is taken on the
#' probability scale BETWEEN the bounds and mapped back through qnorm(). The
#' previous implementation clamped rnorm() draws with pmax()/pmin(), which does
#' not truncate: it piles the tail probability onto the bounds themselves
#' (in the archived population this produced 39 patients aged exactly 18.000,
#' 9 aged exactly 80.000 and 15 with V0 exactly 10.0).
#'
#' Named `rtnorm_bounded` rather than `rtruncnorm` because the CRAN truncnorm
#' package exports `rtruncnorm(n, a, b, mean, sd)` with a DIFFERENT argument
#' order; if that package were ever attached it would silently shadow this
#' helper and misparameterise every call.
#'
#' @param n Number of samples
#' @param mean Mean of normal distribution
#' @param sd Standard deviation
#' @param lower Lower bound
#' @param upper Upper bound
#' @return Vector of n samples
#' @export
rtnorm_bounded <- function(n, mean, sd, lower = -Inf, upper = Inf) {
  p_lo <- pnorm(lower, mean, sd)
  p_hi <- pnorm(upper, mean, sd)
  qnorm(runif(n, p_lo, p_hi), mean, sd)
}

#' Truncated log-normal sampler
#'
#' Genuine truncation by inverse transform on the log scale (see
#' [rtnorm_bounded()] for why clamping was replaced and why the helper was
#' renamed away from `rtrunclnorm`).
#'
#' @param n Number of samples
#' @param meanlog Mean of log-scale (mu)
#' @param sdlog SD of log-scale (sigma)
#' @param lower Lower bound
#' @param upper Upper bound
#' @return Vector of n samples
#' @export
rtlnorm_bounded <- function(n, meanlog, sdlog, lower = -Inf, upper = Inf) {
  log_lo <- if (lower <= 0) -Inf else log(lower)
  log_hi <- if (is.infinite(upper)) Inf else log(upper)
  p_lo <- pnorm(log_lo, meanlog, sdlog)
  p_hi <- pnorm(log_hi, meanlog, sdlog)
  exp(qnorm(runif(n, p_lo, p_hi), meanlog, sdlog))
}

#' Generate virtual population
#'
#' Creates N virtual patients. All nine covariates are drawn INDEPENDENTLY;
#' no correlation structure is imposed.
#'
#' Known simplification: age and eGFR are physiologically correlated (eGFR
#' declines with age), and drawing them independently ignores that.
#'
#' Covariates:
#'   - Age: Normal(45, 15), truncated [18, 80]
#'   - Body weight: Normal(75, 12), truncated [40, 120] kg
#'   - eGFR: Normal(90, 25), truncated [15, 150] mL/min
#'   - Time to treatment: Uniform(1, 7) days from symptom onset
#'   - V0 (initial viral load): LogNormal(log(100), 1), truncated [10, 10000]
#'   - Immune strength multiplier: LogNormal(0, 0.3)
#'   - Endothelial sensitivity: LogNormal(0, 0.3)
#'   - Adaptive immune strength: LogNormal(0, 0.3)
#'   - Baseline PLT: Normal(250000, 50000), truncated [100000, 400000]
#'   - Syndrome phenotype: HFRS (this analysis is HFRS-only)
#'
#' @param N Number of virtual patients (default 1000)
#' @param seed Random seed for reproducibility
#' @return Data frame with one row per patient
#' @export
generate_virtual_population <- function(N = 1000, seed = 42) {
  set.seed(seed)

  pop <- data.frame(
    patient_id       = seq_len(N),
    age_years        = rtnorm_bounded(N, mean = 45, sd = 15, lower = 18, upper = 80),
    body_weight_kg   = rtnorm_bounded(N, mean = 75, sd = 12, lower = 40, upper = 120),
    eGFR_mL_min      = rtnorm_bounded(N, mean = 90, sd = 25, lower = 15, upper = 150),
    time_to_treatment_days = runif(N, min = 1, max = 7),
    V0_copies_mL     = rtlnorm_bounded(N, meanlog = log(100), sdlog = 1,
                                    lower = 10, upper = 10000),
    immune_strength  = rtlnorm_bounded(N, meanlog = 0, sdlog = 0.3,
                                    lower = 0.1, upper = 5),
    endothelial_sensitivity = rtlnorm_bounded(N, meanlog = 0, sdlog = 0.3,
                                           lower = 0.1, upper = 5),
    adaptive_strength = rtlnorm_bounded(N, meanlog = 0, sdlog = 0.3,
                                     lower = 0.1, upper = 5),
    PLT0_baseline    = rtnorm_bounded(N, mean = 250000, sd = 50000,
                                   lower = 100000, upper = 400000),
    # HFRS only. The lung-weighted HCPS mapping was never calibrated to its
    # intended ~40% placebo mortality and is not reported, so it has been
    # removed rather than shipped uncalibrated. This was the last draw in the
    # frame, so removing its runif() leaves every other covariate unchanged.
    syndrome         = rep("HFRS", N),
    stringsAsFactors = FALSE
  )

  pop
}

#' Apply virtual population covariates to model parameters
#'
#' Maps patient-level covariates to parameter modifications.
#'
#' Covariate-parameter mappings:
#'   - eGFR -> CL_RBV (renal clearance scaling)
#'   - body_weight_kg -> Vd_RBV, Vd_FAV (allometric scaling)
#'   - immune_strength -> k_FI, k_CI, k_CF2 (multiplied)
#'   - endothelial_sensitivity -> k_PV, k_PC (multiplied)
#'   - adaptive_strength -> k_CD4_I, k_CD8_I, k_IgM, k_switch (multiplied)
#'   - PLT0_baseline -> PLT_0 (initial condition)
#'   - V0 -> V_0 (initial condition)
#'   - syndrome -> modifies relative weights w_K and w_L
#'
#' @param pop_row Single row from virtual population data frame
#' @param pars Base parameter list from [get_parameters()]
#' @return Modified parameter list
#' @export
apply_vpop_to_model <- function(pop_row, pars) {
  pars_mod <- pars

  # Renal clearance scaling for Ribavirin
  pars_mod$CL_RBV <- pars$CL_RBV_std * (pop_row$eGFR_mL_min / pars$eGFR_ref)
  pars_mod$body_weight_kg <- pop_row$body_weight_kg
  pars_mod$adaptive_strength <- pop_row$adaptive_strength

  # Allometric volume scaling (reference 75 kg)
  wt_factor <- pop_row$body_weight_kg / 75
  pars_mod$Vd_RBV <- pars$Vd_RBV * wt_factor
  pars_mod$Vd_FAV <- pars$Vd_FAV * wt_factor

  # Immune strength multiplier (scaled production rates)
  pars_mod$k_FI  <- pars$k_FI * pop_row$immune_strength
  pars_mod$k_CI  <- pars$k_CI * pop_row$immune_strength
  pars_mod$k_CF2 <- pars$k_CF2 * pop_row$immune_strength

  # Endothelial sensitivity
  pars_mod$k_PV <- pars$k_PV * pop_row$endothelial_sensitivity
  pars_mod$k_PC <- pars$k_PC * pop_row$endothelial_sensitivity

  # Adaptive immune strength multiplier (T cell and B cell activation rates)
  pars_mod$k_CD4_I <- pars$k_CD4_I * pop_row$adaptive_strength
  pars_mod$k_CD8_I <- pars$k_CD8_I * pop_row$adaptive_strength
  pars_mod$k_IgM   <- pars$k_IgM * pop_row$adaptive_strength
  pars_mod$k_switch <- pars$k_switch * pop_row$adaptive_strength

  # Initial conditions
  pars_mod$PLT_0 <- pop_row$PLT0_baseline
  pars_mod$V_0   <- pop_row$V0_copies_mL

  # Syndrome-specific mortality is now handled in compute_clinical_endpoints()

  pars_mod
}

#' Summary statistics for virtual population
#'
#' @param pop Virtual population data frame
#' @return Data frame with mean, sd, min, max, median per covariate
#' @export
summarise_virtual_population <- function(pop) {
  covariates <- c("age_years", "body_weight_kg", "eGFR_mL_min",
                   "time_to_treatment_days", "V0_copies_mL",
                   "immune_strength", "endothelial_sensitivity",
                   "adaptive_strength", "PLT0_baseline")

  do.call(rbind, lapply(covariates, function(col) {
    data.frame(
      covariate = col,
      mean      = mean(pop[[col]], na.rm = TRUE),
      sd        = sd(pop[[col]], na.rm = TRUE),
      min       = min(pop[[col]], na.rm = TRUE),
      max       = max(pop[[col]], na.rm = TRUE),
      median    = median(pop[[col]], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}
