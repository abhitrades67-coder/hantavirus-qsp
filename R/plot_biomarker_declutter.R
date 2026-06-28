#' Biomarker Trajectories — Decluttered
#'
#' Option B: Facet by treatment day. 3 biomarkers (rows) × 7 days (cols) = 21 panels.
#'   Each panel shows 4 median lines (placebo + 3 active arms).
#'   No individual patient spaghetti.
#'
#' Option C: Delta-from-placebo curves.
#'   3 panels (one per biomarker). Each shows (treatment − placebo) curves for
#'   all treatment days × active arms. Color = treatment day, line type = arm.
#'   Directly visualizes when treatment benefit vanishes.
#'
#' Data cached to outputs/biomarker_traj_data.csv for fast re-generation.

suppressPackageStartupMessages({
  library(deSolve)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

source("model/hantavirus_qsp.R")
source("model/parameters.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/simulate_trial.R")
source("R/analysis.R")

pars <- get_parameters()
cache_path <- "outputs/biomarker_traj_data.csv"

# --- Step 1: Simulate or load ---
if (file.exists(cache_path)) {
  cat("Loading cached biomarker trajectory data...\n")
  all_traj <- read.csv(cache_path, stringsAsFactors = FALSE)
} else {
  treatment_days <- 1:7
  arms <- c("placebo", "ribavirin", "favipiravir", "combination")
  t_end <- 21
  dt <- 0.5

  cat("Running biomarker trajectory simulations...\n")
  cat(sprintf("  %d arms x %d days = %d simulations\n",
              length(arms), length(treatment_days),
              length(arms) * length(treatment_days)))

  all_traj <- data.frame(
    time = numeric(), treatment_day = integer(),
    arm = character(),
    C_pro = numeric(), P = numeric(), PLT = numeric(),
    stringsAsFactors = FALSE
  )

  # Helper: simulate and interpolate to common grid
  run_clean_sim <- function(arm, t_start, t_end, dt) {
    sim_out <- simulate_patient(pars, arm, t_start = t_start, t_end = t_end, dt = dt)
    if (is.null(sim_out)) return(NULL)
    common_times <- seq(0, t_end, by = dt)
    interp <- function(x) approx(sim_out$time, x, xout = common_times, rule = 2)$y
    data.frame(
      time  = common_times,
      C_pro = interp(sim_out$C_pro),
      P     = interp(sim_out$P),
      PLT   = interp(sim_out$PLT)
    )
  }

  for (arm in arms) {
    for (tday in treatment_days) {
      cat(sprintf("  %s day %d...\n", arm, tday))
      sim <- run_clean_sim(arm, t_start = tday, t_end = t_end, dt = dt)
      if (!is.null(sim)) {
        sim$treatment_day <- tday
        sim$arm <- arm
        all_traj <- rbind(all_traj, sim)
      }
    }
  }

  # Scale PLT to K/uL for readability
  all_traj$PLT <- all_traj$PLT / 1000

  write.csv(all_traj, cache_path, row.names = FALSE)
  cat(sprintf("  Cached: %s\n", cache_path))
}

# --- Shared aesthetics ---
arm_colors <- c(
  placebo     = "#757575",
  ribavirin   = "#1f78b4",
  favipiravir = "#E69F00",
  combination = "#D55E00"
)

arm_labels <- c(
  placebo     = "Placebo",
  ribavirin   = "Ribavirin",
  favipiravir = "Favipiravir",
  combination = "Combination"
)

all_traj$arm <- factor(all_traj$arm,
  levels = names(arm_colors))

biomarker_labels <- c(
  C_pro = "Pro-inflammatory cytokines (AU)",
  P     = "Permeability index",
  PLT   = "Platelets (K/µL)"
)

# --- Pivot to long for faceting ---
traj_long <- all_traj |>
  pivot_longer(cols = c(C_pro, P, PLT),
               names_to = "biomarker", values_to = "value") |>
  mutate(
    biomarker = factor(biomarker, levels = c("C_pro", "P", "PLT"),
                       labels = biomarker_labels),
    day_label = factor(paste("Day", treatment_day),
                       levels = paste("Day", 1:7))
  )

# ============================================================================
# OPTION B: Facet by treatment day (3 biomarkers × 7 days = 21 panels)
# ============================================================================

cat("Generating Option B (facet by treatment day)...\n")

p_b <- ggplot(traj_long, aes(x = time, y = value, colour = arm)) +
  # Thin lines for each arm (no individual patient spaghetti)
  stat_summary(fun = median, geom = "line", linewidth = 1.0,
               aes(group = arm)) +
  facet_grid(biomarker ~ day_label, scales = "free_y") +
  scale_colour_manual(values = arm_colors, labels = arm_labels) +
  labs(
    x        = "Time (days)",
    y        = "Value (median)",
    colour   = "Treatment arm",
    title    = "Biomarker Trajectories by Treatment Start Day",
    subtitle = "Rows: biomarkers | Columns: treatment start day | Lines: arm medians"
  ) +
  theme_qsp() +
  theme(
    panel.spacing = unit(0.4, "lines"),
    strip.text.y  = element_text(size = 9),
    strip.text.x  = element_text(size = 8),
    legend.position = "bottom"
  )

ggsave_safe("outputs/biomarker_trajectories_faceted.png",
            width = 16, height = 9, dpi = 300)
cat("  Saved: outputs/biomarker_trajectories_faceted.png\n")

# ============================================================================
# OPTION C: Delta-from-placebo curves
# ============================================================================

cat("Generating Option C (delta from placebo)...\n")

# Compute placebo median at each time point for each treatment day
placebo_med <- traj_long |>
  filter(arm == "placebo") |>
  group_by(biomarker, day_label, time) |>
  summarise(placebo_val = median(value), .groups = "drop")

# Join placebo onto treatment data, compute delta
delta_data <- traj_long |>
  filter(arm != "placebo") |>
  left_join(placebo_med,
            by = c("biomarker", "day_label", "time")) |>
  mutate(delta = value - placebo_val) |>
  # Re-level arm for display
  mutate(arm = factor(arm, levels = c("ribavirin", "favipiravir", "combination")))

# Re-label day for color scale
delta_data$day_num <- factor(delta_data$treatment_day, levels = 1:7)

p_c <- ggplot(delta_data,
              aes(x = time, y = delta,
                  colour = day_num, linetype = arm)) +
  geom_hline(yintercept = 0, colour = "#757575", linewidth = 0.5, linetype = "dotted") +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  facet_wrap(~ biomarker, scales = "free_y", ncol = 1) +
  scale_colour_brewer(palette = "YlOrRd", name = "Treatment\nstart day",
                      direction = -1) +
  scale_linetype_manual(
    values = c(ribavirin = "solid", favipiravir = "dashed", combination = "dotdash"),
    labels = arm_labels
  ) +
  labs(
    x        = "Time (days)",
    y        = expression(Delta * " from placebo (treatment " - " placebo)"),
    colour   = "Start day",
    linetype = "Arm",
    title    = "Treatment Effect vs. Placebo: Biomarker Delta Curves",
    subtitle = paste0(
      "Δ = treatment − placebo | ",
      "Below zero = biomarker reduced | Above zero = biomarker elevated | ",
      "Curves converging to zero = loss of treatment effect"
    )
  ) +
  theme_qsp() +
  theme(
    legend.position  = "right",
    legend.key.width = unit(1.2, "cm"),
    strip.text       = element_text(size = 10, face = "bold")
  )

ggsave_safe("outputs/biomarker_trajectories_delta.png",
            width = 12, height = 10, dpi = 300)
cat("  Saved: outputs/biomarker_trajectories_delta.png\n")
cat("Done.\n")
