#' Adaptive Immunity Panel by Treatment Start Day
#'
#' Generates a faceted panel plot where each panel corresponds to a specific
#' treatment start day (1-7 post-symptom-onset). Uses single-patient simulations
#' with a consistent time grid to produce smooth curves.
#'
#' Layout: rows = adaptive compartments (CD8_E, CD4, IgM, IgG)
#'         columns = treatment start days (1-7)
#'
#' Placebo shown as grey reference in all panels.

suppressPackageStartupMessages({
  library(deSolve)
  library(ggplot2)
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
arms <- c("placebo", "ribavirin", "favipiravir", "combination")
t_end <- 21
dt <- 0.5  # consistent grid for smooth lines

cat("Running single-patient adaptive immunity simulations...\n")
cat(sprintf("  Treatment days: %s\n", paste(treatment_days, collapse = ", ")))
cat(sprintf("  Arms: %s\n", paste(arms, collapse = ", ")))

# Collect all trajectory data with a consistent time grid
all_traj <- data.frame(
  time = numeric(), treatment_day = integer(),
  arm = character(),
  CD8_E = numeric(), CD4 = numeric(),
  IgM = numeric(), IgG = numeric(),
  stringsAsFactors = FALSE
)

# Helper: run simulation and interpolate to a common time grid
run_clean_sim <- function(arm, t_start, t_end, dt) {
  # Run with fine dt internally for accuracy
  sim_out <- simulate_patient(pars, arm, t_start = t_start, t_end = t_end, dt = dt)
  if (is.null(sim_out)) return(NULL)

  # Interpolate to the common output grid (0, dt, 2*dt, ..., t_end)
  common_times <- seq(0, t_end, by = dt)

  interp <- function(x) {
    approx(sim_out$time, x, xout = common_times, rule = 2)$y
  }

  data.frame(
    time  = common_times,
    CD8_E = interp(sim_out$CD8_E),
    CD4   = interp(sim_out$CD4),
    IgM   = interp(sim_out$IgM),
    IgG   = interp(sim_out$IgG),
    stringsAsFactors = FALSE
  )
}

# Run placebo once (same for all days)
cat("  Running placebo... ")
placebo_sim <- run_clean_sim("placebo", t_start = 0, t_end = t_end, dt = dt)
cat("done\n")

# Run treatment arms for each day
for (treat_day in treatment_days) {
  cat(sprintf("  Treatment day %d: ", treat_day))
  for (arm in c("ribavirin", "favipiravir", "combination")) {
    sim <- run_clean_sim(arm, t_start = treat_day, t_end = t_end, dt = dt)
    if (is.null(sim)) {
      cat(sprintf("%s FAILED; ", arm))
      next
    }
    sim$arm <- arm
    sim$treatment_day <- treat_day
    all_traj <- rbind(all_traj, sim)
    cat(sprintf("%s ok; ", arm))
  }
  cat("\n")
}

# Add placebo (replicated for each day, shown as grey reference)
placebo_traj <- placebo_sim
placebo_traj$arm <- "placebo"
placebo_traj$treatment_day <- NA_integer_

# Build combined dataset for plotting
# Replicate placebo for each panel
placebo_panels <- do.call(rbind, lapply(treatment_days, function(d) {
  p <- placebo_traj
  p$treatment_day <- d
  p
}))

plot_data <- rbind(all_traj, placebo_panels)
plot_data$arm <- factor(plot_data$arm,
  levels = c("placebo", "ribavirin", "favipiravir", "combination"))
plot_data$treatment_day <- factor(plot_data$treatment_day,
  levels = treatment_days,
  labels = paste("Day", treatment_days))

# Convert to long format for faceting
long_data <- tidyr::pivot_longer(
  plot_data,
  cols = c("CD8_E", "CD4", "IgM", "IgG"),
  names_to = "compartment",
  values_to = "value"
)

long_data$compartment <- factor(long_data$compartment,
  levels = c("CD8_E", "CD4", "IgM", "IgG"),
  labels = c(expression("CD8"^"+"~"effector T cells"),
             expression("CD4"^"+"~"helper T cells"),
             expression("IgM antibodies"),
             expression("IgG antibodies")))

arm_colors <- c(
  placebo      = "#757575",
  ribavirin    = "#1f78b4",
  favipiravir  = "#E69F00",
  combination  = "#D55E00"
)

cat("Generating panel plot...\n")

p <- ggplot(long_data, aes(x = time, y = value, colour = arm, linetype = arm)) +
  geom_line(linewidth = 0.8) +
  facet_grid(compartment ~ treatment_day,
             scales = "free_y",
             labeller = labeller(treatment_day = label_wrap_gen(width = 8))) +
  scale_colour_manual(values = arm_colors, name = "Treatment arm") +
  scale_linetype_manual(
    values = c(placebo = "dashed", ribavirin = "solid",
               favipiravir = "solid", combination = "solid"),
    name = "Treatment arm"
  ) +
  labs(
    x = "Time post symptom onset (days)",
    y = "Value (AU)",
    title = "Adaptive Immunity Trajectories by Treatment Start Day",
    subtitle = "Rows: immune compartments | Columns: treatment initiation day (post symptom onset) | Dashed = placebo"
  ) +
  theme_qsp() +
  theme(
    strip.text.y = element_text(angle = 0, size = 10),
    strip.text.x = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    panel.spacing = unit(0.3, "lines")
  )

ggsave_safe("outputs/adaptive_immunity_by_day.png", width = 16, height = 12, dpi = 300)

cat("Saved: outputs/adaptive_immunity_by_day.png\n")
cat("Done.\n")
