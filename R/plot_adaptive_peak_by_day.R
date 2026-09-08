#' Peak Adaptive Immunity vs Treatment Start Day
#'
#' For each treatment arm and day, extracts peak values of CD8_E, CD4, IgM, IgG.
#' Generates a 4-panel figure: one panel per compartment, X = treatment day, Y = peak.
#' Placebo shown as a horizontal reference line.
#'
#' This captures the "earlier = better" dose-response relationship cleanly,
#' with zero line overlap.

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

treatment_days <- 1:7
arms <- c("ribavirin", "favipiravir", "combination")
t_end <- 21
dt <- 0.5

# One patient per arm/day. `simulate_patient()` is called with the same
# `pars <- get_parameters()` every time and nothing on this path is stochastic
# (no virtual population is drawn here), so any number of replicates would be
# byte-identical and every SD/SE would be exactly 0. Running one deterministic
# representative patient is the honest description of what this figure shows.
n_patients <- 1

cat("Running adaptive immunity simulations for peak analysis...\n")
cat(sprintf("  %d patients x %d arms x %d days = %d simulations\n",
            n_patients, length(arms), length(treatment_days),
            n_patients * length(arms) * length(treatment_days)))

# Collect peak values
results <- data.frame(
  treatment_day = integer(),
  arm = character(),
  patient = integer(),
  CD8_E_peak = numeric(),
  CD4_peak = numeric(),
  IgM_peak = numeric(),
  IgG_peak = numeric(),
  stringsAsFactors = FALSE
)

# Helper: run simulation and extract peaks
run_and_extract <- function(arm, t_start, dt) {
  sim <- simulate_patient(pars, arm, t_start = t_start, t_end = t_end, dt = dt)
  if (is.null(sim)) return(NULL)
  data.frame(
    CD8_E_peak = max(sim$CD8_E),
    CD4_peak   = max(sim$CD4),
    IgM_peak   = max(sim$IgM),
    IgG_peak   = max(sim$IgG)
  )
}

# Run placebo once for reference
cat("  Placebo baseline... ")
placebo_peaks <- data.frame()
for (p in 1:n_patients) {
  pk <- run_and_extract("placebo", t_start = 0, dt = dt)
  if (!is.null(pk)) placebo_peaks <- rbind(placebo_peaks, pk)
}
placebo_summary <- data.frame(
  CD8_E_mean = mean(placebo_peaks$CD8_E_peak),
  CD4_mean   = mean(placebo_peaks$CD4_peak),
  IgM_mean   = mean(placebo_peaks$IgM_peak),
  IgG_mean   = mean(placebo_peaks$IgG_peak)
)
cat("done\n")

# Run treatment arms
for (tday in treatment_days) {
  cat(sprintf("  Day %d: ", tday))
  for (arm in arms) {
    for (p in 1:n_patients) {
      pk <- run_and_extract(arm, t_start = tday, dt = dt)
      if (!is.null(pk)) {
        pk$treatment_day <- tday
        pk$arm <- arm
        pk$patient <- p
        results <- rbind(results, pk)
      }
    }
    cat(sprintf("%s ", arm))
  }
  cat("\n")
}

# Collect the single deterministic peak per arm per day
summary_tbl <- results |>
  group_by(arm, treatment_day) |>
  summarise(
    CD8_E_mean = mean(CD8_E_peak),
    CD4_mean   = mean(CD4_peak),
    IgM_mean   = mean(IgM_peak),
    IgG_mean   = mean(IgG_peak),
    .groups = "drop"
  )

# Convert to long format
summary_long <- summary_tbl |>
  pivot_longer(
    cols = c(CD8_E_mean, CD4_mean, IgM_mean, IgG_mean),
    names_to = "compartment",
    values_to = "peak"
  ) |>
  mutate(
    compartment_raw = gsub("_mean$", "", compartment)
  )

summary_long$arm <- factor(summary_long$arm,
  levels = c("ribavirin", "favipiravir", "combination"))

# Placebo reference data (long format)
placebo_long <- data.frame(
  compartment_raw = c("CD8_E", "CD4", "IgM", "IgG"),
  placebo_mean = c(placebo_summary$CD8_E_mean, placebo_summary$CD4_mean,
                   placebo_summary$IgM_mean, placebo_summary$IgG_mean)
)

# Merge placebo reference
summary_long <- summary_long |>
  left_join(placebo_long, by = "compartment_raw")

# Compartment labels (clean facet names)
comp_labels <- c(
  "CD8_E" = "CD8+ effector T cells",
  "CD4"   = "CD4+ helper T cells",
  "IgM"   = "IgM antibodies",
  "IgG"   = "IgG antibodies"
)

summary_long$comp_label <- factor(comp_labels[summary_long$compartment_raw],
  levels = comp_labels)

# The placebo reference layers are drawn from `placebo_long`, so it must carry
# the faceting variable too. Without `comp_label`, ggplot2 replicates the whole
# reference data frame into every panel and the CD8_E reference line lands in
# the CD4 / IgM / IgG panels, breaking scales = "free_y".
placebo_long$comp_label <- factor(comp_labels[placebo_long$compartment_raw],
  levels = comp_labels)

arm_colors <- c(
  ribavirin    = "#1f78b4",
  favipiravir  = "#E69F00",
  combination  = "#D55E00"
)

arm_labels <- c(
  ribavirin   = "Ribavirin",
  favipiravir = "Favipiravir",
  combination = "Combination"
)

cat("Generating figure...\n")

p <- ggplot(summary_long, aes(x = treatment_day, y = peak, colour = arm)) +
  # Placebo reference line (horizontal, spans full X range)
  geom_hline(data = placebo_long,
             aes(yintercept = placebo_mean),
             linetype = "dashed", colour = "#757575", linewidth = 0.7) +

  # Treatment lines
  geom_line(linewidth = 1.0, position = position_dodge(width = 0.3)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.3)) +

  facet_wrap(~ comp_label, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = arm_colors, labels = arm_labels) +
  scale_x_continuous(
    breaks = treatment_days,
    labels = paste("Day", treatment_days)
  ) +
  labs(
    x = "Treatment start (days post symptom onset)",
    y = "Peak response (AU)",
    colour = "Treatment arm",
    title = "Peak Adaptive Immune Response by Treatment Start Day",
    subtitle = "Single deterministic representative patient per arm and start day; grey dashed = placebo reference"
  ) +
  theme_qsp() +
  theme(
    strip.text = element_text(size = 10),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "bottom"
  )

ggsave_safe("outputs/adaptive_immunity_peak_by_day.png",
            width = 10, height = 8, dpi = 300)

cat("Saved: outputs/adaptive_immunity_peak_by_day.png\n")
cat("Done.\n")
