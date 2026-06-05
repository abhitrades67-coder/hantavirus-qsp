#!/usr/bin/env Rscript
# Fix the combined tornado plot and summary from sensitivity analysis
library(dplyr)
library(ggplot2)

sens_df <- read.csv("outputs/sensitivity_local.csv")

# Tornado data
tornado_data <- sens_df %>%
  filter(abs(perturbation_pct) == 50) %>%
  group_by(parameter, arm) %>%
  summarise(
    S_V_peak_avg     = mean(abs(S_V_peak), na.rm = TRUE),
    S_mortality_avg  = mean(abs(S_mortality), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ungroup()

# Combined tornado plot
p_combined <- ggplot(
  bind_rows(
    tornado_data %>% mutate(endpoint = "V_peak", sensitivity = S_V_peak_avg) %>% select(parameter, arm, endpoint, sensitivity),
    tornado_data %>% mutate(endpoint = "Mortality", sensitivity = S_mortality_avg) %>% select(parameter, arm, endpoint, sensitivity)
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

ggsave("outputs/sensitivity_tornado.png", plot = p_combined,
       width = 12, height = 7, dpi = 200)
cat("Saved: outputs/sensitivity_tornado.png\n")

# Text summary
summary_v <- tornado_data %>%
  select(parameter, arm, S_V_peak_avg) %>%
  arrange(desc(S_V_peak_avg))

summary_m <- tornado_data %>%
  select(parameter, arm, S_mortality_avg) %>%
  arrange(desc(S_mortality_avg))

sink("outputs/sensitivity_summary.txt")

cat("============================================================\n")
cat("  LOCAL SENSITIVITY ANALYSIS — HANTAVIRUS QSP MODEL\n")
cat("============================================================\n")
cat(sprintf("Date: %s\n", Sys.Date()))
cat(sprintf("Perturbations: ±10%%, ±25%%, ±50%%\n"))
cat(sprintf("Arms: placebo, ribavirin\n"))
cat("Syndrome: HFRS\n")
cat("Simulation horizon: 0-21 days\n\n")

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
cat("  |S| < 0.1  : weakly sensitive\n\n")

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
cat("Saved: outputs/sensitivity_summary.txt\n")
cat("\n=== Combined tornado and summary complete ===\n")
