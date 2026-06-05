#' Heatmap: Peak Organ Injury by Treatment Start Day and Arm
#'
#' Extracts peak P (permeability), K (renal), L (lung) injury values.
#' Generates a 3-panel heatmap: rows = treatment day, columns = arm.
#' Color = % of placebo peak (blue = preserved/protected, red = worse than placebo).
#'
#' Data cached to outputs/organ_peak_data.csv for fast re-generation.

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
cache_path     <- "outputs/organ_peak_data.csv"
placebo_cache  <- "outputs/organ_peak_placebo.csv"

# --- Step 1: Load or generate data ---
if (file.exists(cache_path)) {
  cat("Loading cached organ peak data...\n")
  results <- read.csv(cache_path, stringsAsFactors = FALSE)
  placebo_peaks <- read.csv(placebo_cache, stringsAsFactors = FALSE)
} else {
  treatment_days <- 1:7
  arms <- c("ribavirin", "favipiravir", "combination")
  t_end <- 21
  dt <- 0.5
  n_patients <- 10

  cat("Running organ injury simulations for heatmap...\n")
  cat(sprintf("  %d patients x %d arms x %d days = %d simulations\n",
              n_patients, length(arms), length(treatment_days),
              n_patients * length(arms) * length(treatment_days)))

  results <- data.frame(
    treatment_day = integer(), arm = character(), patient = integer(),
    P_peak = numeric(), K_peak = numeric(), L_peak = numeric(),
    stringsAsFactors = FALSE
  )

  run_and_extract <- function(arm, t_start, dt) {
    sim <- simulate_patient(pars, arm, t_start = t_start, t_end = t_end, dt = dt)
    if (is.null(sim)) return(NULL)
    data.frame(
      P_peak = max(sim$P),
      K_peak = max(sim$K),
      L_peak = max(sim$L)
    )
  }

  # Placebo
  cat("  Placebo... ")
  placebo_peaks <- data.frame()
  for (p in 1:n_patients) {
    pk <- run_and_extract("placebo", t_start = 0, dt = dt)
    if (!is.null(pk)) placebo_peaks <- rbind(placebo_peaks, pk)
  }
  write.csv(placebo_peaks, placebo_cache, row.names = FALSE)
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
    P_mean = mean(P_peak),
    K_mean = mean(K_peak),
    L_mean = mean(L_peak),
    .groups = "drop"
  )

# Placebo means
placebo_means <- colMeans(placebo_peaks[, c("P_peak", "K_peak", "L_peak")])
names(placebo_means) <- c("P", "K", "L")

# --- Step 3: Build heatmap data ---
heat_data <- summary_tbl |>
  pivot_longer(cols = c(P_mean, K_mean, L_mean),
               names_to = "compartment", values_to = "peak") |>
  mutate(
    compartment_raw = gsub("_mean$", "", compartment),
    placebo_val = case_when(
      compartment_raw == "P" ~ placebo_means["P"],
      compartment_raw == "K" ~ placebo_means["K"],
      compartment_raw == "L" ~ placebo_means["L"]
    ),
    pct_of_placebo = peak / placebo_val * 100
  )

# Factor ordering
heat_data$arm <- factor(heat_data$arm,
  levels = c("ribavirin", "favipiravir", "combination"),
  labels = c("Ribavirin", "Favipiravir", "Combination"))
heat_data$treatment_day <- factor(heat_data$treatment_day, levels = 1:7,
  labels = paste("Day", 1:7))
heat_data$compartment <- factor(heat_data$compartment_raw,
  levels = c("P", "K", "L"),
  labels = c("Endothelial\npermeability (P)",
             "Renal injury\nindex (K)",
             "Lung injury\nindex (L)"))

# For organ injury, LOWER % of placebo = BETTER (more protection)
# Flip color scale: blue = protection (below 100%), red = worse (above 100%)
max_pct <- max(heat_data$pct_of_placebo, na.rm = TRUE)
min_pct <- min(heat_data$pct_of_placebo, na.rm = TRUE)
lim <- max(abs(c(max_pct - 100, 100 - min_pct)), 5)

cat(sprintf("  Placebo means: P=%.4f  K=%.4f  L=%.4f\n",
            placebo_means["P"], placebo_means["K"], placebo_means["L"]))
cat(sprintf("  Pct of placebo range: %.0f%% – %.0f%%\n", min_pct, max_pct))

# --- Step 4: Generate heatmap ---
p <- ggplot(heat_data,
            aes(x = arm, y = treatment_day, fill = pct_of_placebo)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = sprintf("%.0f%%", pct_of_placebo)),
            size = 3.5, fontface = "bold") +
  facet_wrap(~ compartment, ncol = 3,
             labeller = label_wrap_gen(width = 14)) +
  scale_fill_gradient2(
    low      = "#4575b4",   # blue: protection (lower injury than placebo)
    mid      = "#f7f7f7",   # white: at placebo level
    high     = "#d73027",   # red: worse injury than placebo
    midpoint = 100,
    limits   = c(100 - lim, 100 + lim),
    name     = "% of placebo\npeak injury\n(lower = better)",
    guide    = guide_colourbar(barwidth = 1.2, barheight = 12)
  ) +
  labs(
    x        = "Treatment arm",
    y        = "Treatment start day",
    title    = "Peak Organ Injury: Heatmap by Treatment Arm and Day",
    subtitle = paste0(
      "Values = % of placebo peak injury | ",
      "Blue = organ-protective | Red = worse than placebo | White = at placebo level"
    )
  ) +
  theme_qsp() +
  theme(
    panel.grid   = element_blank(),
    axis.ticks   = element_blank(),
    strip.text   = element_text(size = 9.5, face = "bold"),
    legend.position = "right"
  )

ggsave_safe("outputs/organ_injury_heatmap.png",
            width = 12, height = 5.5, dpi = 300)

cat("Saved: outputs/organ_injury_heatmap.png\n")
cat("Done.\n")
