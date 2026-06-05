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
high_conf   <- param_table$parameter[param_table$confidence == "high"]
medium_conf <- param_table$parameter[param_table$confidence == "medium"]
low_conf    <- param_table$parameter[param_table$confidence == "low"]

# Exclude structural constants from perturbation
exclude <- c("I_0", "T_0", "I_ref", "F_I_ref", "F_II_ref", "C_ref", "V_ref",
             "eGFR_ref", "P_0", "C_pro_0", "K_0", "L_0", "PLT_0")
high_conf   <- setdiff(high_conf, exclude)
medium_conf <- setdiff(medium_conf, exclude)
low_conf    <- setdiff(low_conf, exclude)

cat("=== Uncertainty Quantification ===\n")
cat(sprintf("  High confidence:   %d parameters (CV=10%%)\n", length(high_conf)))
cat(sprintf("  Medium confidence: %d parameters (CV=30%%)\n", length(medium_conf)))
cat(sprintf("  Low confidence:    %d parameters (CV=50%%)\n", length(low_conf)))

# --- Sample parameters --------------------------------------------------------
sample_param <- function(base_val, cv) {
  if (base_val <= 0) return(base_val)
  sigma <- sqrt(log(1 + cv^2))
  mu <- log(base_val) - sigma^2 / 2
  val <- rlnorm(1, mu, sigma)
  # Truncate at ±3 SD
  val <- max(val, base_val * exp(-3 * sigma))
  val <- min(val, base_val * exp(3 * sigma))
  val
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
cat(sprintf("  Total simulations: %d (successful: %d)\n",
            nrow(uq_df), sum(!is.na(uq_df$V_peak))))

# --- Compute uncertainty intervals --------------------------------------------
cat("  Computing uncertainty intervals...\n")

endpoints <- c("V_peak", "V_AUC", "K_peak", "L_peak", "C_pro_max", "mortality_prob")

uq_results <- data.frame()
for (a in c("placebo", "ribavirin")) {
  arm_data <- uq_df %>% filter(arm == a)
  for (ep_name in endpoints) {
    vals <- arm_data[[ep_name]]
    vals <- vals[is.finite(vals) & !is.na(vals)]
    if (length(vals) < 3) next

    uq_results <- rbind(uq_results, data.frame(
      arm      = a,
      endpoint = ep_name,
      median   = median(vals),
      p5       = quantile(vals, 0.05, names = FALSE),
      p95      = quantile(vals, 0.95, names = FALSE),
      cv       = sd(vals) / mean(vals) * 100
    ))
  }
}

write.csv(uq_results, "outputs/uncertainty_intervals.csv", row.names = FALSE)

cat("\n=== Uncertainty Intervals ===\n")
for (a in c("placebo", "ribavirin")) {
  cat(sprintf("\n  --- %s ---\n", toupper(a)))
  arm_res <- uq_results %>% filter(arm == a)
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
