#!/usr/bin/env Rscript
# Generate VPC plot from pre-computed vpc_data.csv
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
})

vpc_all <- read.csv("outputs/vpc_data.csv")
cat(sprintf("Loaded %d rows from vpc_data.csv\n", nrow(vpc_all)))

# Compute prediction bands per arm
bands_all <- vpc_all %>%
  group_by(arm, time) %>%
  summarise(
    p5      = quantile(V, 0.05, na.rm = TRUE),
    median  = quantile(V, 0.50, na.rm = TRUE),
    p95     = quantile(V, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

# Handle log10 zeros
min_v <- min(bands_all$median[bands_all$median > 0], na.rm = TRUE)
bands_all$median_safe <- pmax(bands_all$median, min_v / 10)
bands_all$p5_safe     <- pmax(bands_all$p5, min_v / 10)
bands_all$p95_safe    <- pmax(bands_all$p95, min_v / 10)
bands_all$median_log  <- log10(bands_all$median_safe)
bands_all$p5_log      <- log10(bands_all$p5_safe)
bands_all$p95_log     <- log10(bands_all$p95_safe)

# Individual trajectories (downsample)
vpc_all$logV <- log10(pmax(vpc_all$V, min_v / 10))
show_ids <- sort(unique(vpc_all$sample_id))
show_ids <- show_ids[seq(1, length(show_ids), by = max(1, floor(length(show_ids)/15)))]
vpc_sub <- vpc_all[vpc_all$sample_id %in% show_ids, ]

# Placebo VPC
p_placebo <- ggplot() +
  geom_ribbon(data = bands_all %>% filter(arm == "placebo"),
              aes(x = time, ymin = p5_log, ymax = p95_log),
              fill = "#2166AC", alpha = 0.25) +
  geom_line(data = vpc_sub %>% filter(arm == "placebo"),
            aes(x = time, y = logV, group = sample_id),
            color = "gray60", alpha = 0.3, linewidth = 0.3) +
  geom_line(data = bands_all %>% filter(arm == "placebo"),
            aes(x = time, y = median_log),
            color = "#2166AC", linewidth = 1.5) +
  annotate("rect", xmin = 3, xmax = 5, ymin = 7, ymax = 8,
           fill = "blue", alpha = 0.08) +
  annotate("text", x = 4, y = 8.2,
           label = "Clinical peak\nrange (HFRS)",
           size = 3, color = "blue") +
  scale_y_continuous(name = "Viral load (log10 copies/mL)",
                     breaks = seq(0, 10, by = 1),
                     labels = function(x) paste0("10^", x)) +
  scale_x_continuous(name = "Time (days)", breaks = seq(0, 21, by = 3)) +
  ggtitle("Placebo", subtitle = "N=50 parameter sets") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

# Ribavirin VPC
p_ribavirin <- ggplot() +
  geom_ribbon(data = bands_all %>% filter(arm == "ribavirin"),
              aes(x = time, ymin = p5_log, ymax = p95_log),
              fill = "#B2182B", alpha = 0.25) +
  geom_line(data = vpc_sub %>% filter(arm == "ribavirin"),
            aes(x = time, y = logV, group = sample_id),
            color = "gray60", alpha = 0.3, linewidth = 0.3) +
  geom_line(data = bands_all %>% filter(arm == "ribavirin"),
            aes(x = time, y = median_log),
            color = "#B2182B", linewidth = 1.5) +
  geom_vline(xintercept = 3, color = "red", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 3.5, y = max(bands_all$median_log[bands_all$arm == "ribavirin"]) - 0.5,
           label = "Treatment start", color = "red", size = 3.5, hjust = 0) +
  scale_y_continuous(name = "Viral load (log10 copies/mL)",
                     breaks = seq(0, 10, by = 1),
                     labels = function(x) paste0("10^", x)) +
  scale_x_continuous(name = "Time (days)", breaks = seq(0, 21, by = 3)) +
  ggtitle("Ribavirin", subtitle = "N=50 parameter sets, t_start = 3 days") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave("outputs/vpc_plot.png", plot = arrangeGrob(p_placebo, p_ribavirin, nrow = 2),
       width = 8, height = 10, dpi = 200)
cat("Saved: outputs/vpc_plot.png\n")

# Summary
cat("\n=== VPC Summary ===\n")
for (a in c("placebo", "ribavirin")) {
  cat(sprintf("\n--- %s ---\n", toupper(a)))
  ab <- bands_all %>% filter(arm == a)
  peak_idx <- which.max(ab$median)
  cat(sprintf("  Peak VL (median): %.2e copies/mL at day %.1f\n",
              ab$median[peak_idx], ab$time[peak_idx]))
  cat(sprintf("  Peak 95%% PI: [%.2e, %.2e]\n",
              ab$p5[peak_idx], ab$p95[peak_idx]))
  clearance <- min(ab$time[ab$p95 < 100], na.rm = TRUE)
  if (is.finite(clearance)) {
    cat(sprintf("  Time to clearance (p95 <100): day %.1f\n", clearance))
  } else {
    cat("  Time to clearance: not reached within 21 days\n")
  }
}

cat("\n=== VPC plot complete ===\n")
