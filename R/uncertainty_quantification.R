#!/usr/bin/env Rscript
# ============================================================================
# Uncertainty Quantification — Hantavirus QSP Model (clean rewrite)
# ============================================================================

suppressPackageStartupMessages({
  library(deSolve)
  library(dplyr)
  library(ggplot2)
})

source("model/parameters.R")
source("model/hantavirus_qsp.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/simulate_trial.R")

set.seed(2026)
N_SAMPLES <- 50
t_end  <- 21
dt     <- 0.1

# --- Define parameter confidence classes --------------------------------------
pars_base <- get_parameters()

# Read parameter table for confidence ratings
param_table <- read.csv("outputs/parameter_table.csv", stringsAsFactors = FALSE)

# Every confidence class in the table must have an explicit treatment. Silently
# matching only "high"/"medium"/"low" meant the three parameters rated "defined"
# fell through and were NEVER perturbed, with no warning.
KNOWN_CONFIDENCE <- c("high", "medium", "low", "defined", "structural")
unknown_rows <- !(param_table$confidence %in% KNOWN_CONFIDENCE)
if (any(unknown_rows)) {
  stop(sprintf(
    paste0("Unrecognised confidence class(es) in outputs/parameter_table.csv: %s\n",
           "  offending parameters: %s\n",
           "  Add an explicit treatment for each class before running."),
    paste(unique(param_table$confidence[unknown_rows]), collapse = ", "),
    paste(param_table$parameter[unknown_rows], collapse = ", ")))
}

high_conf   <- param_table$parameter[param_table$confidence == "high"]
medium_conf <- param_table$parameter[param_table$confidence == "medium"]
low_conf    <- param_table$parameter[param_table$confidence == "low"]

# Confidence class "defined": ablation-control switches whose nominal value is
# 1.0. They are not measured quantities but on/off scale factors for the
# ribavirin non-antiviral terms, so perturbing them would blur the ablation
# contrast rather than express uncertainty. They are therefore excluded
# deliberately and visibly (below) instead of being dropped by accident.
ablation_controls <- param_table$parameter[param_table$confidence == "defined"]

# Confidence class "structural": numerical regularisation constants that have no
# measured counterpart (currently n_mem_switch, the steepness of the smoothed
# effector-to-memory decay transition). Perturbing them would vary the numerical
# treatment of a threshold rather than a biological quantity, so they are
# excluded deliberately and visibly, like the ablation controls.
numerical_constants <- param_table$parameter[param_table$confidence == "structural"]

# Exclude structural constants and the ablation-control switches from perturbation
exclude <- c("I_0", "T_0", "I_ref", "F_I_ref", "F_II_ref", "C_ref", "V_ref",
             "eGFR_ref", "P_0", "C_pro_0", "K_0", "L_0", "PLT_0",
             ablation_controls, numerical_constants)
high_conf   <- setdiff(high_conf, exclude)
medium_conf <- setdiff(medium_conf, exclude)
low_conf    <- setdiff(low_conf, exclude)

# Sample in a fixed (alphabetical) name order, NOT the row order of
# parameter_table.csv. Positional assignment of the RNG stream meant that
# regenerating or reordering that table silently changed every sampled
# parameter set even under the same seed.
high_conf   <- sort(high_conf)
medium_conf <- sort(medium_conf)
low_conf    <- sort(low_conf)

n_perturbed <- length(unique(c(high_conf, medium_conf, low_conf)))
n_excluded  <- length(intersect(unique(param_table$parameter), exclude))

cat("=== Uncertainty Quantification ===\n")
cat(sprintf("  High confidence:   %d parameters (CV=10%%)\n", length(high_conf)))
cat(sprintf("  Medium confidence: %d parameters (CV=30%%)\n", length(medium_conf)))
cat(sprintf("  Low confidence:    %d parameters (CV=50%%)\n", length(low_conf)))
cat(sprintf("  Perturbed: %d parameters; deliberately excluded: %d (%d structural constants + %d ablation controls: %s)\n",
            n_perturbed, n_excluded,
            n_excluded - length(intersect(unique(param_table$parameter), ablation_controls)),
            length(intersect(unique(param_table$parameter), ablation_controls)),
            paste(ablation_controls, collapse = ", ")))

# --- Sample parameters --------------------------------------------------------
# mu = log(base_val), NOT log(base_val) - sigma^2/2. The mean-preserving form
# gives E[X] = base_val but median(X) = base_val * exp(-sigma^2/2), i.e. it
# shifts EVERY median down (-1.2% at CV=10%, -4.3% at CV=30%, -10.6% at CV=50%)
# simultaneously and all in the same direction. The intent is that the
# calibrated value is the central (median) value of the sampling distribution.
sample_param <- function(base_val, cv) {
  if (base_val <= 0) return(base_val)
  sigma <- sqrt(log(1 + cv^2))
  mu <- log(base_val)
  lo <- base_val * exp(-3 * sigma)
  hi <- base_val * exp( 3 * sigma)
  # Genuine truncation at +/-3 SD by rejection sampling. Clamping with
  # max()/min() piles the tail probability onto the bounds instead of
  # removing it. Rejection probability here is ~0.3%, so 100 attempts is
  # ample; the clamp survives only as an unreachable fallback.
  for (attempt in seq_len(100)) {
    val <- rlnorm(1, mu, sigma)
    if (val >= lo && val <= hi) return(val)
  }
  warning(sprintf("sample_param: no in-bounds draw in 100 attempts (base=%g, cv=%g); falling back to the bound",
                  base_val, cv))
  min(max(rlnorm(1, mu, sigma), lo), hi)
}

param_sets <- vector("list", N_SAMPLES)
for (i in seq_len(N_SAMPLES)) {
  pars <- pars_base
  for (p in high_conf) {
    if (p %in% names(pars)) pars[[p]] <- sample_param(pars[[p]], 0.10)
  }
  for (p in medium_conf) {
    if (p %in% names(pars)) pars[[p]] <- sample_param(pars[[p]], 0.30)
  }
  for (p in low_conf) {
    if (p %in% names(pars)) pars[[p]] <- sample_param(pars[[p]], 0.50)
  }
  # Restore platelet homeostasis. model/parameters.R (~line 152) requires
  # k_PLT_prod / k_PLT_loss == PLT_0 for the disease-free steady state, but the
  # two rates are perturbed independently while PLT_0 is held fixed, so the
  # disease-free baseline drifts. The renal-injury term k_KH * (1 - PLT/PLT_0)
  # fires ONLY when PLT < PLT_0, so that drift is one-directional: it can only
  # inflate renal injury and mortality, never deflate them. Recomputing PLT_0
  # keeps the identity while retaining the uncertainty in both rates.
  if (!is.null(pars$k_PLT_prod) && !is.null(pars$k_PLT_loss) &&
      is.finite(pars$k_PLT_loss) && pars$k_PLT_loss > 0) {
    pars$PLT_0 <- pars$k_PLT_prod / pars$k_PLT_loss
  }
  param_sets[[i]] <- pars
}

# Save parameter sets
param_sets_df <- do.call(rbind, lapply(param_sets, as.data.frame))
write.csv(param_sets_df, "outputs/uq_parameter_sets.csv", row.names = FALSE)
cat(sprintf("  Saved: outputs/uq_parameter_sets.csv\n"))

# --- Run simulations (PARALLEL) ----------------------------------------------
cat("  Running simulations...\n")

n_cores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterExport(cl,
  varlist = c("param_sets", "N_SAMPLES", "t_end", "dt"),
  envir = environment()
)

endpoint_list <- parallel::parLapply(cl, seq_len(N_SAMPLES), function(i) {
  source("model/hantavirus_qsp.R", local = TRUE)
  source("model/parameters.R", local = TRUE)
  source("R/pk_models.R", local = TRUE)
  source("R/pd_models.R", local = TRUE)
  source("R/simulate_trial.R", local = TRUE)

  results <- list()
  idx <- 1
  for (arm in c("placebo", "ribavirin")) {
    t_start <- if (arm == "placebo") 0 else 3
    sim_out <- tryCatch(
      simulate_patient(pars = param_sets[[i]], arm = arm, t_start = t_start,
                       t_end = t_end, dt = dt),
      error = function(e) NULL
    )
    if (is.null(sim_out) || is.null(sim_out$V)) {
      results[[idx]] <- data.frame(
        sample = i, arm = arm, V_peak = NA_real_, V_AUC = NA_real_,
        K_peak = NA_real_, L_peak = NA_real_, C_pro_max = NA_real_,
        mortality_prob = NA_real_
      )
      idx <- idx + 1
      next
    }
    pk  <- extract_peak_biomarkers(sim_out)
    auc <- compute_viral_AUC(sim_out)
    ep  <- extract_endpoints_from_simulation(sim_out, param_sets[[i]],
                                             syndrome = "HFRS")
    results[[idx]] <- data.frame(
      sample = i, arm = arm, V_peak = pk$V_peak, V_AUC = auc,
      K_peak = pk$K_peak, L_peak = pk$L_peak, C_pro_max = pk$C_pro_peak,
      mortality_prob = ep$mortality_prob
    )
    idx <- idx + 1
  }
  do.call(rbind, results)
})

uq_df <- do.call(rbind, Filter(function(x) is.data.frame(x), endpoint_list))

# --- Classify runs: established vs non-established (burn-out) -----------------
# A burn-out run is NOT a failure: simulate_patient() returns a valid disease-free
# trajectory with V identically 0, which is finite and so was previously counted
# as "successful". The endpoint distribution is therefore a MIXTURE of a point
# mass at zero and the established distribution; summarising it as if unimodal is
# what drove p5 = 0 for every endpoint and CVs of 69-162%.
# Threshold matches the repo convention in R/preexposure_analysis.R.
V_ESTABLISH_THRESHOLD <- 1e4   # copies/mL
# Classify by the PLACEBO arm of the same parameter sample, not per run.
# Classifying each run on its own V_peak would drop exactly those treated runs
# whose viral load was pushed below the threshold BY THE TREATMENT - i.e. the
# best responders - and would leave the arms summarising different subsets of
# the design. That is the same collider-conditioning error the GSA script now
# reports explicitly. Establishment is a property of the parameter set.
placebo_est <- uq_df$sample[uq_df$arm == "placebo" &
                              is.finite(uq_df$V_peak) &
                              uq_df$V_peak > V_ESTABLISH_THRESHOLD]
uq_df$established <- uq_df$sample %in% placebo_est

n_solved <- sum(!is.na(uq_df$V_peak))
n_est    <- sum(uq_df$established)
cat(sprintf("  Total simulations: %d (solver returned a result: %d)\n",
            nrow(uq_df), n_solved))
cat(sprintf("  Infection established: %d/%d (%.1f%%); burn-out (V_peak <= %.0e): %d\n",
            n_est, n_solved, 100 * n_est / max(n_solved, 1),
            V_ESTABLISH_THRESHOLD, n_solved - n_est))

# --- Compute uncertainty intervals --------------------------------------------
cat("  Computing uncertainty intervals...\n")

endpoints <- c("V_peak", "V_AUC", "K_peak", "L_peak", "C_pro_max", "mortality_prob")

# Statistics are computed on the ESTABLISHED subset only, so they describe one
# distribution rather than a bimodal mixture. The zeros are not dropped
# silently: n_total, n_established and establishment_fraction are carried in
# every row of the output file, and the subset column labels the conditioning.
uq_results <- data.frame()
for (a in c("placebo", "ribavirin")) {
  arm_data  <- uq_df %>% filter(arm == a)
  n_arm_tot <- sum(!is.na(arm_data$V_peak))
  arm_est   <- arm_data %>% filter(established)
  n_arm_est <- nrow(arm_est)
  est_frac  <- n_arm_est / max(n_arm_tot, 1)
  for (ep_name in endpoints) {
    vals <- arm_est[[ep_name]]
    vals <- vals[is.finite(vals) & !is.na(vals)]
    if (length(vals) < 3) next

    uq_results <- rbind(uq_results, data.frame(
      arm      = a,
      endpoint = ep_name,
      subset   = "established",
      median   = median(vals),
      p5       = quantile(vals, 0.05, names = FALSE),
      p95      = quantile(vals, 0.95, names = FALSE),
      cv       = sd(vals) / mean(vals) * 100,
      n        = length(vals),
      n_total  = n_arm_tot,
      n_established = n_arm_est,
      establishment_fraction = est_frac
    ))
  }
}

write.csv(uq_results, "outputs/uncertainty_intervals.csv", row.names = FALSE)

cat("\n=== Uncertainty Intervals (CONDITIONAL ON INFECTION ESTABLISHMENT) ===\n")
for (a in c("placebo", "ribavirin")) {
  cat(sprintf("\n  --- %s ---\n", toupper(a)))
  arm_res <- uq_results %>% filter(arm == a)
  if (nrow(arm_res) > 0) {
    cat(sprintf("    established: %d/%d runs (%.1f%%); the remaining %d burn-out runs\n",
                arm_res$n_established[1], arm_res$n_total[1],
                100 * arm_res$establishment_fraction[1],
                arm_res$n_total[1] - arm_res$n_established[1]))
    cat("    are a point mass at V_peak = 0 and are excluded from the statistics below.\n")
  }
  for (j in seq_len(nrow(arm_res))) {
    r <- arm_res[j, ]
    cat(sprintf("    %-15s: median=%.3g  [%.3g, %.3g]  CV=%.1f%%\n",
                r$endpoint, r$median, r$p5, r$p95, r$cv))
  }
}

# --- Scatter plot -------------------------------------------------------------
p1 <- ggplot(uq_df, aes(x = factor(arm), y = V_peak, fill = arm)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  scale_y_log10() +
  scale_fill_manual(values = c("placebo" = "#888888", "ribavirin" = "#2166AC")) +
  labs(x = "Arm", y = "Peak Viral Load (copies/mL)",
       title = "Uncertainty Quantification: V_peak",
       subtitle = sprintf("N=%d Monte Carlo samples", N_SAMPLES)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

p2 <- ggplot(uq_df, aes(x = factor(arm), y = mortality_prob, fill = arm)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = c("placebo" = "#888888", "ribavirin" = "#2166AC")) +
  labs(x = "Arm", y = "Mortality Probability",
       title = "Uncertainty Quantification: Mortality",
       subtitle = sprintf("N=%d Monte Carlo samples", N_SAMPLES)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

library(gridExtra)
ggsave("outputs/uq_scatter.png", plot = arrangeGrob(p1, p2, nrow = 2),
       width = 8, height = 10, dpi = 200)
cat("\n  Saved: outputs/uq_scatter.png\n")

cat("\n=== Uncertainty quantification complete ===\n")
