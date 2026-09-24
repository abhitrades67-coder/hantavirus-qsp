#' Heatmap: Peak Adaptive Immune Response by Treatment Arm and Day
#'
#' 4 panels (one per immune compartment). Rows = treatment start day, columns = arm.
#' Color = mean peak response. Placebo level shown as a separate reference row.
#' Cell text shows the peak value.
#'
#' Data cached to outputs/adaptive_peak_data.csv to avoid re-simulation.

suppressPackageStartupMessages({
  library(deSolve)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

source("model/hantavirus_qsp.R")
source("model/parameters.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/simulate_trial.R")
source("R/analysis.R")

pars <- get_parameters()
cache_path <- "outputs/adaptive_peak_data.csv"

# --- Step 1: Load or generate data ---
if (file.exists(cache_path)) {
  cat("Loading cached peak data...\n")
  results <- read.csv(cache_path, stringsAsFactors = FALSE)
  placebo_peaks <- read.csv("outputs/adaptive_peak_placebo.csv", stringsAsFactors = FALSE)
} else {
  treatment_days <- 1:7
  arms <- c("ribavirin", "favipiravir", "combination")
  t_end <- 21
  dt <- 0.5
  # One representative patient at nominal parameter values. Earlier
  # versions ran ten identical copies of it and described them as ten
  # patients, which made every spread statistic structurally zero.
  n_patients <- 1

  cat("Running adaptive immunity simulations for heatmap...\n")
  cat(sprintf("  %d patients x %d arms x %d days = %d simulations\n",
              n_patients, length(arms), length(treatment_days),
              n_patients * length(arms) * length(treatment_days)))

  results <- data.frame(
    treatment_day = integer(), arm = character(), patient = integer(),
    CD8_E_peak = numeric(), CD4_peak = numeric(),
    IgM_peak = numeric(), IgG_peak = numeric(),
    stringsAsFactors = FALSE
  )

  run_and_extract <- function(arm, t_start, dt) {
    sim <- simulate_patient(pars, arm, t_start = t_start, t_end = t_end, dt = dt)
    if (is.null(sim)) return(NULL)
    data.frame(CD8_E_peak = max(sim$CD8_E), CD4_peak = max(sim$CD4),
               IgM_peak = max(sim$IgM), IgG_peak = max(sim$IgG))
  }

  # Placebo
  cat("  Placebo... ")
  placebo_peaks <- data.frame()
  for (p in 1:n_patients) {
    pk <- run_and_extract("placebo", t_start = 0, dt = dt)
    if (!is.null(pk)) placebo_peaks <- rbind(placebo_peaks, pk)
  }
  write.csv(placebo_peaks, "outputs/adaptive_peak_placebo.csv", row.names = FALSE)
  cat("done\n")

  # Treatment arms
  for (tday in treatment_days) {
    cat(sprintf("  Day %d: ", tday))
    for (arm in arms) {
      for (p in 1:n_patients) {
        pk <- run_and_extract(arm, t_start = tday, dt = dt)
        if (!is.null(pk)) {
          pk$treatment_day <- tday; pk$arm <- arm; pk$patient <- p
          results <- rbind(results, pk)
        }
      }
      cat(sprintf("%s ", arm))
    }
    cat("\n")
  }
  write.csv(results, cache_path, row.names = FALSE)
  cat(sprintf("  Cached: %s\n", cache_path))
}

# --- Step 2: Summarise means ---
summary_tbl <- results |>
  group_by(arm, treatment_day) |>
  summarise(
    CD8_E = mean(CD8_E_peak), CD4 = mean(CD4_peak),
    IgM   = mean(IgM_peak),   IgG = mean(IgG_peak),
    .groups = "drop"
  )

# Placebo means
placebo_means <- data.frame(
  CD8_E = mean(placebo_peaks$CD8_E_peak),
  CD4   = mean(placebo_peaks$CD4_peak),
  IgM   = mean(placebo_peaks$IgM_peak),
  IgG   = mean(placebo_peaks$IgG_peak)
)

# --- Step 3: Build heatmap data ---
# Pivot to long, then compute % of placebo for color
heat_data <- summary_tbl |>
  pivot_longer(cols = c(CD8_E, CD4, IgM, IgG),
               names_to = "compartment", values_to = "peak") |>
  mutate(
    placebo_val = case_when(
      compartment == "CD8_E" ~ placebo_means$CD8_E,
      compartment == "CD4"   ~ placebo_means$CD4,
      compartment == "IgM"   ~ placebo_means$IgM,
      compartment == "IgG"   ~ placebo_means$IgG
    ),
    pct_of_placebo = peak / placebo_val * 100
  )

# Factor ordering
heat_data$arm <- factor(heat_data$arm,
  levels = c("ribavirin", "favipiravir", "combination"),
  labels = c("Ribavirin", "Favipiravir", "Combination"))
heat_data$treatment_day <- factor(heat_data$treatment_day, levels = 1:7,
  labels = paste("Day", 1:7))
heat_data$compartment <- factor(heat_data$compartment,
  levels = c("CD8_E", "CD4", "IgM", "IgG"),
  labels = c("CD8+ effector\nT cells", "CD4+ helper\nT cells",
             "IgM\nantibodies", "IgG\nantibodies"))

# Dual color: blue-white-red diverge around 100% of placebo
# Above 100% = immune preservation (blue), Below 100% = immune suppression (red)
max_pct <- max(heat_data$pct_of_placebo, na.rm = TRUE)
min_pct <- min(heat_data$pct_of_placebo, na.rm = TRUE)
# Asymmetric limits about the midpoint (see the note in plot_organ_heatmap.R).
fill_lo <- min(min_pct, 95)
fill_hi <- max(max_pct, 105)

cat(sprintf("  Pct of placebo range: %.0f%% – %.0f%%\n", min_pct, max_pct))

# --- Step 4: Heatmap ---
p <- ggplot(heat_data,
            aes(x = arm, y = treatment_day, fill = pct_of_placebo)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = sprintf("%.0f%%", pct_of_placebo)),
            size = 3.5, fontface = "bold") +
  facet_wrap(~ compartment, ncol = 2) +
  scale_fill_gradient2(
    low      = "#d73027",   # red: immune suppression below placebo
    mid      = "#f7f7f7",   # white: at placebo level
    high     = "#4575b4",   # blue: immune preservation above placebo
    midpoint = 100,
    limits   = c(fill_lo, fill_hi),
    # NOTE: the colour direction here is the OPPOSITE of Figure 3 -- in this figure
    # blue means ABOVE placebo (response preserved). Say so, or a reader
    # comparing the two adjacent heatmaps will misread this one.
    name     = "% of placebo\npeak response\n(blue = above placebo)",
    guide    = guide_colourbar(barwidth = 1.2, barheight = 12)
  ) +
  labs(
    x        = "Treatment arm",
    y        = "Treatment start day",
    # No in-figure title/subtitle: the caption is supplied in the manuscript
    # (journal style). Legend + axis labels carry the interpretation.
    title    = NULL
  ) +
  theme_qsp() +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right"
  )

ggsave_safe("outputs/adaptive_immunity_heatmap.png",
            width = 9, height = 7, dpi = 600)

cat("Saved: outputs/adaptive_immunity_heatmap.png\n")
cat("Done.\n")
