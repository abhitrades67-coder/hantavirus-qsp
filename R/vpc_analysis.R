#' Visual Predictive Check (VPC) Analysis
#'
#' Generates a visual predictive check for the Hantavirus QSP model,
#' showing prediction uncertainty bands for viral load trajectories.
#'
#' Requires: outputs/uq_parameter_sets.csv (from uncertainty_quantification.R)
#'
#' Usage: Rscript R/vpc_analysis.R
#'
#' @return Saves outputs/vpc_plot.png and outputs/vpc_data.csv

# ============================================================================
# 0. SETUP
# ============================================================================
library(deSolve)
library(dplyr)
library(ggplot2)

# NOTE: this script previously called setwd() on a hardcoded Google Drive path,
# which made it unrunnable on any machine but the author's -- including from the
# Zenodo archive. Run it from the project root, like every other script here.
if (!file.exists("DESCRIPTION") || !dir.exists("model")) {
  stop("Run this script from the Hantavirus_QSP project root.")
}

# Source all model components
source("model/parameters.R")
source("model/hantavirus_qsp.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/virtual_population.R")
source("R/simulate_trial.R")

cat("=== Visual Predictive Check (VPC) Analysis ===\n")

# ============================================================================
# 1. LOAD PARAMETER SETS
# ============================================================================
# Load pre-generated parameter sets from UQ analysis
uq_file <- "outputs/uq_parameter_sets.csv"
if (!file.exists(uq_file)) {
  cat("Parameter sets not found. Running uncertainty quantification first...\n")
  source("R/uncertainty_quantification.R")
}

param_sets_df <- read.csv(uq_file, stringsAsFactors = FALSE)
N_SAMPLES <- nrow(param_sets_df)
cat(sprintf("Loaded %d parameter sets from %s\n", N_SAMPLES, uq_file))

pars <- get_parameters()

# ============================================================================
# 2. RUN SIMULATIONS FOR VPC
# ============================================================================
run_single_sim_vpc <- function(pars_i, arm, t_start = 3, t_end = 21, dt = 0.5) {
  simulate_patient(pars_i, arm, t_start = t_start, t_end = t_end, dt = dt)
}

# Run simulations for both arms (PARALLEL)
cat("\nRunning VPC simulations...\n")

n_cores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterExport(cl,
  varlist = c("param_sets_df", "N_SAMPLES", "pars"),
  envir = environment()
)

sim_list <- parallel::parLapply(cl, seq_len(N_SAMPLES), function(s) {
  source("model/hantavirus_qsp.R", local = TRUE)
  source("model/parameters.R", local = TRUE)
  source("R/pk_models.R", local = TRUE)
  source("R/pd_models.R", local = TRUE)
  source("R/virtual_population.R", local = TRUE)
  source("R/simulate_trial.R", local = TRUE)

  pars_i <- pars
  for (pname in setdiff(names(param_sets_df), "sample_id")) {
    if (pname %in% names(pars_i)) pars_i[[pname]] <- param_sets_df[s, pname]
  }

  results <- list()
  for (arm in c("placebo", "ribavirin")) {
    t_start_arm <- if (arm == "placebo") 0 else 3
    out <- tryCatch({
      simulate_patient(pars_i, arm, t_start = t_start_arm, t_end = 21, dt = 0.5)
    }, error = function(e) NULL)

    if (!is.null(out)) {
      if (!any(is.nan(out$V)) && !any(is.infinite(out$V))) {
        out$arm <- arm
        out$sample_id <- s
        results[[arm]] <- out
      }
    }
  }
  do.call(rbind, results)
})

vpc_all <- do.call(rbind, Filter(function(x) is.data.frame(x), sim_list))

# Split by arm for compatibility with plotting code
vpc_data_list <- list(
  placebo   = vpc_all[vpc_all$arm == "placebo", ],
  ribavirin = vpc_all[vpc_all$arm == "ribavirin", ]
)

# Save VPC data
write.csv(vpc_all, "outputs/vpc_data.csv", row.names = FALSE)
cat("\nSaved: outputs/vpc_data.csv\n")

# ============================================================================
# 4. GENERATE VPC PLOT
# ============================================================================
cat("\nGenerating VPC plot...\n")

# Compute prediction bands for each arm
compute_vpc_bands <- function(arm_data) {
  # Group by time, compute quantiles
  arm_data %>%
    group_by(time) %>%
    summarise(
      p5      = quantile(V, 0.05, na.rm = TRUE),
      p25     = quantile(V, 0.25, na.rm = TRUE),
      median  = quantile(V, 0.50, na.rm = TRUE),
      p75     = quantile(V, 0.75, na.rm = TRUE),
      p95     = quantile(V, 0.95, na.rm = TRUE),
      mean_v  = mean(V, na.rm = TRUE),
      n_sim   = n(),
      .groups = "drop"
    )
}

bands_placebo <- compute_vpc_bands(vpc_data_list$placebo)
bands_placebo$arm <- "placebo"

bands_ribavirin <- compute_vpc_bands(vpc_data_list$ribavirin)
bands_ribavirin$arm <- "ribavirin"

bands_all <- bind_rows(bands_placebo, bands_ribavirin)

# Add log10-transformed values (handle zeros by setting to minimum positive)
min_positive_v <- min(vpc_all$V[vpc_all$V > 0], na.rm = TRUE)
vpc_all$V_safe  <- pmax(vpc_all$V, min_positive_v / 10)
vpc_all$log10V  <- log10(vpc_all$V_safe)

bands_all$median_safe  <- pmax(bands_all$median, min_positive_v / 10)
bands_all$p5_safe      <- pmax(bands_all$p5, min_positive_v / 10)
bands_all$p95_safe     <- pmax(bands_all$p95, min_positive_v / 10)
bands_all$median_log10 <- log10(bands_all$median_safe)
bands_all$p5_log10     <- log10(bands_all$p5_safe)
bands_all$p95_log10    <- log10(bands_all$p95_safe)

# Clinical reference ranges (approximate, from literature)
# HFRS placebo: peak ~10^7-10^8 copies/mL, declining over 7-14 days
clinical_ref <- data.frame(
  arm = c("placebo", "placebo", "ribavirin", "ribavirin"),
  label = c("Clinical peak range (HFRS)", "Clinical peak range (HFRS)",
            "Expected: accelerated clearance", "Expected: accelerated clearance"),
  day_min = c(3, 7, 5, 7),
  day_max = c(5, 14, 7, 14),
  log10V_min = c(7, 5, 5, 3),
  log10V_max = c(8, 7, 7, 4)
)

# Build the two-panel VPC plot
arm_labels <- c("placebo" = "Placebo", "ribavirin" = "Ribavirin (t_start = 3 d)")

# Treatment start annotation
treatment_line <- data.frame(
  arm = c("placebo", "ribavirin"),
  t_start = c(NA, 3)
)

# Individual trajectories: downsample for visual clarity (show every 5th)
individual_ids <- sort(unique(vpc_all$sample_id))
# Show a subset to avoid overplotting — every Nth trajectory
show_every <- max(1, floor(length(individual_ids) / 15))
show_ids <- individual_ids[seq(1, length(individual_ids), by = show_every)]

vpc_individual <- vpc_all[vpc_all$sample_id %in% show_ids, ]

# Create the plot
p <- ggplot() +
  # Prediction interval band (5-95%)
  geom_ribbon(
    data = bands_all,
    aes(x = time, ymin = p5_log10, ymax = p95_log10, fill = arm),
    alpha = 0.25
  ) +
  # Median line
  geom_line(
    data = bands_all,
    aes(x = time, y = median_log10, color = arm, linetype = "Median"),
    linewidth = 1.2
  ) +
  # Individual simulation trajectories (thin, transparent)
  geom_line(
    data = vpc_individual,
    aes(x = time, y = log10V, group = interaction(arm, sample_id)),
    color = "gray60", alpha = 0.3, linewidth = 0.3
  ) +
  # Treatment start line for ribavirin
  geom_vline(
    data = treatment_line[!is.na(treatment_line$t_start), ],
    aes(xintercept = t_start),
    color = "red", linetype = "dashed", linewidth = 0.8
  ) +
  # Clinical reference annotation
  annotate("rect",
    xmin = 3, xmax = 5, ymin = 7, ymax = 8,
    fill = "blue", alpha = 0.08
  ) +
  annotate("text",
    x = 4, y = 8.2,
    label = "Clinical peak\nrange (HFRS)",
    size = 3, color = "blue"
  ) +
  # Facet by arm
  facet_wrap(~arm, labeller = as_labeller(arm_labels), ncol = 1, scales = "free_y") +
  # Formatting
  scale_y_continuous(
    name = expression("Viral load" ~ (log["10"] ~ "copies/mL")),
    breaks = seq(0, 10, by = 1),
    labels = function(x) paste0("10^", x)
  ) +
  scale_x_continuous(
    name = "Time (days)",
    breaks = seq(0, 21, by = 3)
  ) +
  scale_fill_manual(
    values = c("placebo" = "#2166AC", "ribavirin" = "#B2182B"),
    guide = "none"
  ) +
  scale_color_manual(
    values = c("placebo" = "#2166AC", "ribavirin" = "#B2182B"),
    name = NULL
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c("Median" = 1),
    labels = c("Model median")
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  labs(
    title = "Visual Predictive Check — Viral Load Trajectories",
    subtitle = sprintf("N = %d Monte Carlo parameter sets | Shaded: 5th-95th percentile prediction band | Gray lines: individual simulations",
                       N_SAMPLES)
  )

# Save the plot
ggsave("outputs/vpc_plot.png", p, width = 8, height = 10, dpi = 300, bg = "white")
cat("Saved: outputs/vpc_plot.png\n")

# ============================================================================
# 5. VPC SUMMARY STATISTICS
# ============================================================================
cat("\n=== VPC Summary ===\n")
for (arm in c("placebo", "ribavirin")) {
  cat(sprintf("\n--- %s ---\n", arm))
  arm_bands <- bands_all[bands_all$arm == arm, ]

  # Find peak
  peak_idx <- which.max(arm_bands$median)
  cat(sprintf("  Peak viral load (median): %.2e copies/mL at day %.1f\n",
              arm_bands$median[peak_idx], arm_bands$time[peak_idx]))
  cat(sprintf("  Peak 95%% PI: [%.2e, %.2e]\n",
              arm_bands$p5[peak_idx], arm_bands$p95[peak_idx]))

  # Time to clearance (< 100 copies/mL)
  clearance_time <- min(arm_bands$time[arm_bands$p95 < 100], na.rm = TRUE)
  if (is.finite(clearance_time)) {
    cat(sprintf("  Time to 95%% clearance (<100 copies/mL): day %.1f\n", clearance_time))
  } else {
    cat("  Time to 95%% clearance: not reached within simulation window\n")
  }
}

cat("\n=== Done ===\n")
