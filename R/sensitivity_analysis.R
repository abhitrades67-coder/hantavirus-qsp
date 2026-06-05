#!/usr/bin/env Rscript
# ============================================================================
# Local Sensitivity Analysis — Hantavirus QSP Model
# ============================================================================
#
# Perturbs key parameters by ±10%, ±25%, ±50% and computes normalized
# sensitivity coefficients for viral, organ-injury, and mortality endpoints.
#
# Usage: Rscript R/sensitivity_analysis.R
#
# Outputs:
#   outputs/sensitivity_local.csv       — full results table
#   outputs/sensitivity_tornado.png     — tornado plot (V_peak & mortality)
#   outputs/sensitivity_summary.txt     — ranked parameter summary
# ============================================================================

suppressPackageStartupMessages({
  library(deSolve)
  library(dplyr)
  library(ggplot2)
})

# --- Source model components ------------------------------------------------
source("model/parameters.R")
source("model/hantavirus_qsp.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/simulate_trial.R")

# --- Configuration -----------------------------------------------------------
set.seed(42)

perturbation_pcts <- c(10, 25, 50)
t_end             <- 21
dt                <- 0.1
syndrome          <- "HFRS"

sensitivity_params <- c(
  "beta",          # infection rate
  "p",             # viral production
  "c",             # viral clearance
  "delta_natural", # infected cell natural clearance
  "delta_NK",      # NK-mediated clearance
  "delta_CD8",     # CD8-mediated clearance (adaptive)
  "k_CC",          # cytokine auto-amplification
  "K_CC",          # cytokine auto-amplification half-max
  "d_C",           # cytokine decay
  "k_CD8_I",       # CD8 activation by infected cells (adaptive)
  "k_neut_IgG",    # IgG neutralization rate (adaptive)
  "k_KP",          # renal injury per permeability
  "k_LP",          # lung injury per permeability
  "d_K",           # renal recovery
  "d_L",           # lung recovery
  "EC50_RBV",      # ribavirin potency
  "Emax_RBV",      # ribavirin max effect
  "CL_RBV"         # ribavirin clearance
)

arms <- c("placebo", "ribavirin")
rbavirin_t_start <- 3

# --- Baseline simulation -----------------------------------------------------
cat("=== Baseline simulations ===\n")
pars_base <- get_parameters()

run_sim <- function(pars, arm, t_start = rbavirin_t_start) {
  simulate_patient(
    pars    = pars,
    arm     = arm,
    t_start = t_start,
    t_end   = t_end,
    dt      = dt
  )
}

extract_endpoints <- function(sim_out, pars, syndrome = "HFRS") {
  if (is.null(sim_out)) return(NULL)

  pk  <- extract_peak_biomarkers(sim_out)
  auc <- compute_viral_AUC(sim_out)
  ep  <- extract_endpoints_from_simulation(sim_out, pars, syndrome = syndrome)

  list(
    V_peak       = pk$V_peak,
    V_AUC        = auc,
    K_peak       = pk$K_peak,
    L_peak       = pk$L_peak,
    C_pro_max    = pk$C_pro_peak,
    mortality_prob = ep$mortality_prob
  )
}

sim_base_placebo <- run_sim(pars_base, "placebo", t_start = 0)
sim_base_rbv     <- run_sim(pars_base, "ribavirin", t_start = rbavirin_t_start)

base_placebo <- extract_endpoints(sim_base_placebo, pars_base, syndrome)
base_rbv     <- extract_endpoints(sim_base_rbv, pars_base, syndrome)

cat(sprintf("  Placebo: V_peak=%.1f, K_peak=%.3f, L_peak=%.3f, mort=%.4f\n",
            base_placebo$V_peak, base_placebo$K_peak, base_placebo$L_peak,
            base_placebo$mortality_prob))
cat(sprintf("  Ribavirin: V_peak=%.1f, K_peak=%.3f, L_peak=%.3f, mort=%.4f\n",
            base_rbv$V_peak, base_rbv$K_peak, base_rbv$L_peak,
            base_rbv$mortality_prob))

# --- Sensitivity loop (PARALLEL) ---------------------------------------------
cat("\n=== Sensitivity analysis ===\n")

# Build task grid
task_grid <- expand.grid(
  param     = sensitivity_params,
  pct       = perturbation_pcts,
  direction = c(-1, 1),
  arm       = arms,
  stringsAsFactors = FALSE
)
cat(sprintf("  Total tasks: %d\n", nrow(task_grid)))

# Parallel setup
n_cores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterExport(cl,
  varlist = c("task_grid", "pars_base", "t_end", "dt", "syndrome",
               "rbavirin_t_start", "arms", "sensitivity_params",
               "perturbation_pcts"),
  envir = environment()
)

endpoint_list <- parallel::parLapply(cl, seq_len(nrow(task_grid)), function(i) {
  source("model/hantavirus_qsp.R", local = TRUE)
  source("model/parameters.R", local = TRUE)
  source("R/pk_models.R", local = TRUE)
  source("R/pd_models.R", local = TRUE)
  source("R/simulate_trial.R", local = TRUE)

  row <- task_grid[i, ]
  pert_val <- pars_base[[row$param]] * (1 + row$direction * row$pct / 100)
  if (pert_val < 0) pert_val <- 1e-12

  pars_pert <- pars_base
  pars_pert[[row$param]] <- pert_val

  t_start_arm <- if (row$arm == "placebo") 0 else rbavirin_t_start
  sim_out <- simulate_patient(pars = pars_pert, arm = row$arm,
                              t_start = t_start_arm, t_end = t_end, dt = dt)
  if (is.null(sim_out)) return(NULL)

  pk  <- extract_peak_biomarkers(sim_out)
  auc <- compute_viral_AUC(sim_out)
  ep  <- extract_endpoints_from_simulation(sim_out, pars_pert, syndrome = syndrome)

  data.frame(
    parameter        = row$param,
    perturbation_pct = row$direction * row$pct,
    arm              = row$arm,
    V_peak           = pk$V_peak,
    V_AUC            = auc,
    K_peak           = pk$K_peak,
    L_peak           = pk$L_peak,
    C_pro_max        = pk$C_pro_peak,
    mortality_prob   = ep$mortality_prob,
    stringsAsFactors = FALSE
  )
})

sens_df <- do.call(rbind, Filter(function(x) is.data.frame(x), endpoint_list))
cat(sprintf("  Completed: %d/%d simulations\n", nrow(sens_df), nrow(task_grid)))

# --- Compute sensitivity coefficients ----------------------------------------
# Compute baseline endpoints for normalization
sim_base_placebo <- simulate_patient(pars = pars_base, arm = "placebo",
                                     t_start = 0, t_end = t_end, dt = dt)
sim_base_rbv     <- simulate_patient(pars = pars_base, arm = "ribavirin",
                                     t_start = rbavirin_t_start, t_end = t_end, dt = dt)

get_endpoints <- function(sim_out, pars) {
  pk  <- extract_peak_biomarkers(sim_out)
  auc <- compute_viral_AUC(sim_out)
  ep  <- extract_endpoints_from_simulation(sim_out, pars, syndrome = syndrome)
  c(V_peak = pk$V_peak, V_AUC = auc, K_peak = pk$K_peak,
    L_peak = pk$L_peak, C_pro_max = pk$C_pro_peak,
    mortality_prob = ep$mortality_prob)
}

base_placebo <- get_endpoints(sim_base_placebo, pars_base)
base_rbv     <- get_endpoints(sim_base_rbv, pars_base)

sens_df$S_V_peak    <- NA_real_
sens_df$S_V_AUC     <- NA_real_
sens_df$S_K_peak    <- NA_real_
sens_df$S_L_peak    <- NA_real_
sens_df$S_C_pro_max <- NA_real_
sens_df$S_mortality <- NA_real_

for (r in seq_len(nrow(sens_df))) {
  base_vals <- if (sens_df$arm[r] == "placebo") base_placebo else base_rbv
  param <- sens_df$parameter[r]
  delta_ratio <- (pars_base[[param]] * (1 + sens_df$perturbation_pct[r] / 100) -
                  pars_base[[param]]) / pars_base[[param]]
  if (abs(delta_ratio) < .Machine$double.eps) next

  for (ep in c("V_peak", "V_AUC", "K_peak", "L_peak", "C_pro_max", "mortality_prob")) {
    base_v <- base_vals[ep]
    if (abs(base_v) < .Machine$double.eps) next
    pert_v <- sens_df[[ep]][r]
    s_col <- if (ep == "mortality_prob") "S_mortality" else paste0("S_", ep)
    sens_df[[s_col]][r] <- ((pert_v - base_v) / base_v) / delta_ratio
  }
}

# --- Save CSV ----------------------------------------------------------------
out_dir <- "outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

write.csv(sens_df, file.path(out_dir, "sensitivity_local.csv"), row.names = FALSE)
cat(sprintf("\nSaved: %s\n", file.path(out_dir, "sensitivity_local.csv")))

# --- Tornado plot ------------------------------------------------------------
cat("Generating tornado plots...\n")

tornado_data <- sens_df %>%
  filter(abs(perturbation_pct) == 50) %>%
  group_by(parameter, arm) %>%
  summarise(
    S_V_peak_avg     = mean(abs(S_V_peak), na.rm = TRUE),
    S_mortality_avg  = mean(abs(S_mortality), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ungroup()

# V_peak tornado
tornado_data_v <- tornado_data %>%
  arrange(desc(S_V_peak_avg)) %>%
  mutate(parameter = factor(parameter, levels = unique(parameter)))

p_v <- ggplot(tornado_data_v, aes(x = S_V_peak_avg, y = parameter, fill = arm)) +
  geom_col(position = "dodge", width = 0.7, alpha = 0.85) +
  scale_fill_manual(
    values = c("placebo" = "#888888", "ribavirin" = "#2166AC"),
    labels = c("Placebo", "Ribavirin")
  ) +
  labs(
    x = "|S_Vpeak|  (normalized sensitivity coefficient)",
    y = "Parameter",
    title = "Local Sensitivity: Peak Viral Load (V_peak)",
    subtitle = "Absolute normalized sensitivity at ±50% perturbation",
    fill = "Arm"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "sensitivity_tornado_vpeak.png"), plot = p_v,
       width = 9, height = 6, dpi = 200)
cat(sprintf("  Saved: %s\n", file.path(out_dir, "sensitivity_tornado_vpeak.png")))

# Mortality tornado
tornado_data_m <- tornado_data %>%
  arrange(desc(S_mortality_avg)) %>%
  mutate(parameter = factor(parameter, levels = unique(parameter)))

p_m <- ggplot(tornado_data_m, aes(x = S_mortality_avg, y = parameter, fill = arm)) +
  geom_col(position = "dodge", width = 0.7, alpha = 0.85) +
  scale_fill_manual(
    values = c("placebo" = "#888888", "ribavirin" = "#2166AC"),
    labels = c("Placebo", "Ribavirin")
  ) +
  labs(
    x = "|S_mortality|  (normalized sensitivity coefficient)",
    y = "Parameter",
    title = "Local Sensitivity: Mortality Probability",
    subtitle = "Absolute normalized sensitivity at ±50% perturbation",
    fill = "Arm"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "sensitivity_tornado_mortality.png"), plot = p_m,
       width = 9, height = 6, dpi = 200)
cat(sprintf("  Saved: %s\n", file.path(out_dir, "sensitivity_tornado_mortality.png")))

# Combined tornado
p_combined <- ggplot(
  bind_rows(
    tornado_data %>% select(parameter, arm, sensitivity = S_V_peak_avg) %>% mutate(endpoint = "V_peak"),
    tornado_data %>% select(parameter, arm, sensitivity = S_mortality_avg) %>% mutate(endpoint = "Mortality")
  ) %>%
    mutate(
      endpoint = factor(endpoint, levels = c("V_peak", "Mortality")),
      parameter = factor(parameter)
    ),
  aes(x = sensitivity, y = parameter, fill = arm)
) +
  geom_col(position = "dodge", width = 0.7, alpha = 0.85) +
  facet_wrap(~ endpoint, scales = "free_x", ncol = 2) +
  scale_fill_manual(
    values = c("placebo" = "#888888", "ribavirin" = "#2166AC"),
    labels = c("Placebo", "Ribavirin")
  ) +
  labs(
    x = "|Normalized sensitivity coefficient|",
    y = "Parameter",
    title = "Local Sensitivity Analysis — Hantavirus QSP Model",
    subtitle = "Absolute sensitivity at ±50% perturbation for V_peak and Mortality",
    fill = "Arm"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )

ggsave(file.path(out_dir, "sensitivity_tornado.png"), plot = p_combined,
       width = 12, height = 7, dpi = 200)
cat(sprintf("  Saved: %s\n", file.path(out_dir, "sensitivity_tornado.png")))

# --- Text summary ------------------------------------------------------------
cat("Writing summary...\n")

summary_v <- tornado_data %>%
  select(parameter, arm, S_V_peak_avg) %>%
  arrange(desc(S_V_peak_avg))

summary_m <- tornado_data %>%
  select(parameter, arm, S_mortality_avg) %>%
  arrange(desc(S_mortality_avg))

sink(file.path(out_dir, "sensitivity_summary.txt"))

cat("============================================================\n")
cat("  LOCAL SENSITIVITY ANALYSIS — HANTAVIRUS QSP MODEL\n")
cat("============================================================\n")
cat(sprintf("Date: %s\n", Sys.Date()))
cat(sprintf("Perturbations: ±%s%%\n", paste(perturbation_pcts, collapse = ", ")))
cat(sprintf("Arms: %s\n", paste(arms, collapse = ", ")))
cat(sprintf("Syndrome: %s\n", syndrome))
cat(sprintf("Simulation horizon: 0–%d days\n", t_end))
cat("\n")

cat("------------------------------------------------------------\n")
cat("  RANKING by |S| for V_peak (±50% perturbation)\n")
cat("------------------------------------------------------------\n")
cat(sprintf("%-20s  %-12s  %s\n", "Parameter", "Arm", "|S_V_peak|"))
cat(sprintf("%-20s  %-12s  %s\n", "---------", "---", "----------"))
for (i in seq_len(nrow(summary_v))) {
  cat(sprintf("%-20s  %-12s  %.4f\n",
              summary_v$parameter[i], summary_v$arm[i], summary_v$S_V_peak_avg[i]))
}

cat("\n")
cat("------------------------------------------------------------\n")
cat("  RANKING by |S| for mortality_prob (±50% perturbation)\n")
cat("------------------------------------------------------------\n")
cat(sprintf("%-20s  %-12s  %s\n", "Parameter", "Arm", "|S_mortality|"))
cat(sprintf("%-20s  %-12s  %s\n", "---------", "---", "-------------"))
for (i in seq_len(nrow(summary_m))) {
  cat(sprintf("%-20s  %-12s  %.4f\n",
              summary_m$parameter[i], summary_m$arm[i], summary_m$S_mortality_avg[i]))
}

cat("\n")
cat("------------------------------------------------------------\n")
cat("  INTERPRETATION\n")
cat("------------------------------------------------------------\n")
cat("  |S| > 1.0  : highly sensitive (output changes more than parameter)\n")
cat("  0.1 < |S| < 1.0 : moderately sensitive\n")
cat("  |S| < 0.1  : weakly sensitive\n")
cat("\n")

top5_v <- head(summary_v %>% filter(arm == "placebo"), 5)
top5_m <- head(summary_m %>% filter(arm == "placebo"), 5)

cat("  Top 5 drivers of V_peak (placebo):\n")
for (i in seq_len(nrow(top5_v))) {
  cat(sprintf("    %d. %s  (|S| = %.3f)\n", i, top5_v$parameter[i], top5_v$S_V_peak_avg[i]))
}

cat("\n  Top 5 drivers of mortality (placebo):\n")
for (i in seq_len(nrow(top5_m))) {
  cat(sprintf("    %d. %s  (|S| = %.3f)\n", i, top5_m$parameter[i], top5_m$S_mortality_avg[i]))
}

cat("\n============================================================\n")

sink()

cat(sprintf("  Saved: %s\n", file.path(out_dir, "sensitivity_summary.txt")))

cat("\n=== Sensitivity analysis complete ===\n")
