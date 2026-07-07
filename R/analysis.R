#' Analysis and Plotting Functions
#'
#' Publication-quality visualization and summary statistics for the
#' Hantavirus QSP virtual trial.
#'
#' @name analysis
NULL

#' String concatenation operator (helper)
#'
#' @param a First string
#' @param b Second string
#' @return Concatenated string
`%s+%` <- function(a, b) paste0(a, b)

#' Theme for publication-quality plots
#'
#' @return ggplot2 theme object
#' @export
theme_qsp <- function() {
  ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.border     = ggplot2::element_rect(fill = NA, colour = "grey80"),
      legend.position  = "bottom",
      legend.key       = ggplot2::element_blank(),
      axis.text        = ggplot2::element_text(colour = "black"),
      axis.title       = ggplot2::element_text(face = "bold"),
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5)
    )
}

#' Plot viral load trajectories by treatment arm
#'
#' @param sim_data List of simulation outputs (one per arm), each a data frame
#'   with columns: time, V, arm, patient_id
#' @param outfile Output file path (PNG)
#' @export
plot_viral_kinetics <- function(sim_data, outfile = "outputs/viral_kinetics.png") {
  arm_colors <- c(
    placebo      = "#757575",
    ribavirin    = "#1f78b4",
    favipiravir  = "#E69F00",
    combination  = "#D55E00"
  )

  sim_data$arm <- factor(sim_data$arm,
    levels = c("placebo", "ribavirin", "favipiravir", "combination"))

  # Build facet labels: "Untreated" for day 0, "Day 1" ... "Day 7"
  sim_data$start_label <- ifelse(
    sim_data$treat_start_day == 0,
    "Untreated",
    paste0("Day ", sim_data$treat_start_day)
  )
  sim_data$start_label <- factor(sim_data$start_label,
    levels = c("Untreated", paste0("Day ", 1:7)))

  ggplot2::ggplot(sim_data, ggplot2::aes(x = time, y = V, colour = arm)) +
    # IQR ribbon (drawn behind the median)
    ggplot2::stat_summary(
      fun.min = function(x) stats::quantile(x, 0.25),
      fun.max = function(x) stats::quantile(x, 0.75),
      fun     = median,
      geom    = "ribbon",
      ggplot2::aes(fill = arm),
      alpha   = 0.2,
      colour  = NA
    ) +
    # Median line
    ggplot2::stat_summary(
      fun = median, geom = "line", linewidth = 1.0,
      ggplot2::aes(group = arm)
    ) +
    ggplot2::facet_wrap(~ start_label, ncol = 4) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(values = arm_colors) +
    ggplot2::scale_fill_manual(values = arm_colors, guide = "none") +
    ggplot2::labs(
      x       = "Time (days)",
      y       = expression("Viral RNA (copies/mL)"),
      colour  = "Treatment arm",
      title   = NULL
    ) +
    theme_qsp()

  ggsave_safe(outfile, width = 14, height = 8, dpi = 300)
  message(sprintf("Saved: %s", outfile))
}

#' Plot biomarker trajectories (platelets, cytokines, permeability)
#'
#' @param sim_data Wide-format data frame with columns: time, C, P, PLT, arm
#' @param outfile Output file path (PNG)
#' @export
plot_biomarker_trajectories <- function(sim_data,
                                         outfile = "outputs/biomarker_trajectories.png") {
  # Manual color palette with high contrast for all 4 arms
  arm_colors <- c(
    placebo      = "#757575",  # grey
    ribavirin    = "#1f78b4",  # blue
    favipiravir  = "#E69F00",  # orange
    combination  = "#D55E00"   # vermillion
  )

  # Convert to long format internally
  long_data <- tidyr::pivot_longer(
    sim_data,
    cols = c("C_pro", "P", "PLT"),
    names_to = "biomarker",
    values_to = "value"
  )

  # Scale PLT to thousands for readability
  long_data$value[long_data$biomarker == "PLT"] <-
    long_data$value[long_data$biomarker == "PLT"] / 1000

  # Ensure arm is a factor with consistent order
  long_data$arm <- factor(long_data$arm,
    levels = c("placebo", "ribavirin", "favipiravir", "combination"))

  p <- ggplot2::ggplot(long_data, ggplot2::aes(
    x = time, y = value, colour = arm, group = interaction(arm, patient_id)
  )) +
    # Individual patient lines (thin, transparent, drawn first)
    ggplot2::geom_line(alpha = 0.12, linewidth = 0.2) +
    # Median lines (thick, drawn last so they sit on top)
    ggplot2::stat_summary(
      fun = median, geom = "line", linewidth = 2.0,
      ggplot2::aes(group = arm)
    ) +
    ggplot2::facet_wrap(~biomarker, scales = "free_y", ncol = 1,
                         labeller = ggplot2::as_labeller(c(
                           C_pro = "Pro-inflammatory cytokines (AU)",
                           P = "Permeability index",
                           PLT = "Platelets (K/uL)"
                         ))) +
    ggplot2::labs(
      x = "Time (days)",
      y = "Value",
      colour = "Treatment arm",
      title = "Biomarker Trajectories"
    ) +
    theme_qsp() +
    ggplot2::scale_colour_manual(values = arm_colors)

  ggsave_safe(outfile, width = 10, height = 10, dpi = 300)
  message(sprintf("Saved: %s", outfile))
}

#' Plot organ injury trajectories (renal K and lung L)
#'
#' @param sim_data Wide-format data frame with columns: time, K, L, arm
#' @param outfile Output file path (PNG)
#' @export
plot_organ_injury <- function(sim_data,
                               outfile = "outputs/organ_injury.png") {
  arm_colors <- c(
    placebo      = "#757575",
    ribavirin    = "#1f78b4",
    favipiravir  = "#E69F00",
    combination  = "#D55E00"
  )

  long_data <- tidyr::pivot_longer(
    sim_data,
    cols = c("K", "L"),
    names_to = "organ",
    values_to = "value"
  )

  long_data$arm <- factor(long_data$arm,
    levels = c("placebo", "ribavirin", "favipiravir", "combination"))

  p <- ggplot2::ggplot(long_data, ggplot2::aes(
    x = time, y = value, colour = arm, group = interaction(arm, patient_id)
  )) +
    ggplot2::geom_line(alpha = 0.12, linewidth = 0.2) +
    ggplot2::stat_summary(
      fun = median, geom = "line", linewidth = 2.0,
      ggplot2::aes(group = arm)
    ) +
    ggplot2::facet_wrap(~organ, scales = "free_y", ncol = 1,
                         labeller = ggplot2::as_labeller(c(
                           K = "Renal injury index",
                           L = "Lung injury index"
                         ))) +
    ggplot2::geom_hline(
      yintercept = 0.7, linetype = "dashed", colour = "red", alpha = 0.6
    ) +
    ggplot2::annotate(
      "text", x = 1, y = 0.72, label = "Clinical threshold",
      colour = "red", size = 3, hjust = 0
    ) +
    ggplot2::labs(
      x = "Time (days)",
      y = "Injury index (0-1)",
      colour = "Treatment arm",
      title = "Organ-Specific Injury Trajectories"
    ) +
    theme_qsp() +
    ggplot2::scale_colour_manual(values = arm_colors)

  ggsave_safe(outfile, width = 10, height = 6, dpi = 300)
  message(sprintf("Saved: %s", outfile))
}

#' Plot PK concentration profiles
#'
#' @param sim_data Data frame with columns: time, C_RBV, C_FAV, C_FAVI_RTP, arm
#' @param outfile Output file path (PNG)
#' @export
plot_pk_profiles <- function(sim_data,
                              outfile = "outputs/pk_profiles.png") {
  arm_colors <- c(
    placebo      = "#757575",
    ribavirin    = "#1f78b4",
    favipiravir  = "#E69F00",
    combination  = "#D55E00"
  )

  # Ribavirin profile
  p1 <- ggplot2::ggplot(sim_data[sim_data$arm %in% c("ribavirin", "combination"), ],
                         ggplot2::aes(x = time, y = C_RBV, colour = arm)) +
    ggplot2::geom_line(linewidth = 0.5, alpha = 0.3) +
    ggplot2::stat_summary(fun = median, geom = "line", linewidth = 2.0) +
    ggplot2::labs(
      x = "Time (days)",
      y = expression("Ribavirin (" * mu * "g/mL)"),
      title = "Ribavirin Plasma Concentration"
    ) +
    theme_qsp() +
    ggplot2::scale_colour_manual(values = arm_colors)

  # Favipiravir profile
  p2 <- ggplot2::ggplot(sim_data[sim_data$arm %in% c("favipiravir", "combination"), ],
                         ggplot2::aes(x = time, y = C_FAV, colour = arm)) +
    ggplot2::geom_line(linewidth = 0.5, alpha = 0.3) +
    ggplot2::stat_summary(fun = median, geom = "line", linewidth = 2.0) +
    ggplot2::labs(
      x = "Time (days)",
      y = expression("Favipiravir (" * mu * "g/mL)"),
      title = "Favipiravir Plasma Concentration"
    ) +
    theme_qsp() +
    ggplot2::scale_colour_manual(values = arm_colors)

  # RTP metabolite profile
  p3 <- ggplot2::ggplot(sim_data[sim_data$arm %in% c("favipiravir", "combination"), ],
                         ggplot2::aes(x = time, y = C_FAVI_RTP, colour = arm)) +
    ggplot2::geom_line(linewidth = 0.5, alpha = 0.3) +
    ggplot2::stat_summary(fun = median, geom = "line", linewidth = 2.0) +
    ggplot2::labs(
      x = "Time (days)",
      y = expression("Favipiravir-RTP (" * mu * "g/mL)"),
      title = "Favipiravir Active Metabolite (RTP)"
    ) +
    theme_qsp() +
    ggplot2::scale_colour_manual(values = arm_colors)

  # Combine
  gridExtra::grid.arrange(p1, p2, p3, ncol = 1)
  ggsave_safe(outfile, width = 10, height = 12, dpi = 300)
  message(sprintf("Saved: %s", outfile))
}

#' Plot treatment-window analysis
#'
#' @param window_data Data frame from [run_window_analysis()]
#' @param outfile Output file path (PNG)
#' @export
plot_treatment_window <- function(window_data,
                                   outfile = "outputs/treatment_window.png") {
  arm_colors <- c(
    placebo      = "#757575",
    ribavirin    = "#1f78b4",
    favipiravir  = "#E69F00",
    combination  = "#D55E00"
  )

  window_data$arm <- factor(window_data$arm,
    levels = c("placebo", "ribavirin", "favipiravir", "combination"))

  p <- ggplot2::ggplot(window_data, ggplot2::aes(
    x = factor(treatment_day), y = V_AUC, fill = arm
  )) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::labs(
      x = "Treatment Start Day",
      y = "Viral Load AUC (copies·day/mL)",
      fill = "Treatment arm",
      title = NULL
    ) +
    theme_qsp() +
    ggplot2::scale_fill_manual(values = arm_colors) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0))

  ggsave_safe(outfile, width = 10, height = 6, dpi = 300)
  message(sprintf("Saved: %s", outfile))
}

#' Plot comprehensive treatment window sensitivity analysis
#'
#' Creates a 3-panel figure showing the sensitivity of viral load,
#' organ injury, and mortality to treatment start day.
#'
#' @param window_data Data frame from [run_window_analysis()]
#' @param trial_output Output from [run_virtual_trial()] (for patient-level data)
#' @param outfile Output file path (PNG)
#' @export
plot_window_sensitivity <- function(window_data,
                                     trial_output = NULL,
                                     outfile = "outputs/treatment_window_sensitivity.png") {
  arm_colors <- c(
    placebo     = "#757575",
    ribavirin   = "#1f78b4",
    favipiravir = "#E69F00",
    combination = "#D55E00"
  )

  window_data$arm <- factor(window_data$arm,
    levels = c("ribavirin", "favipiravir", "combination"))

  # Panel A: Mortality vs treatment day (line plot with ribbons)
  p1 <- ggplot2::ggplot(window_data, ggplot2::aes(
    x = treatment_day, y = mortality_prob * 100, colour = arm
  )) +
    ggplot2::geom_line(linewidth = 1.5) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_colour_manual(values = arm_colors[c("ribavirin", "favipiravir", "combination")]) +
    ggplot2::labs(
      x = NULL,
      y = "Mortality (%)",
      colour = "Arm",
      title = "A. Treatment Start Day vs Mortality"
    ) +
    theme_qsp() +
    ggplot2::scale_x_continuous(breaks = 1:7)

  # Panel B: Peak organ injury (K and L) vs treatment day
  p2 <- ggplot2::ggplot(window_data, ggplot2::aes(
    x = treatment_day, colour = arm
  )) +
    ggplot2::geom_line(ggplot2::aes(y = K_peak), linewidth = 1.5) +
    ggplot2::geom_line(ggplot2::aes(y = L_peak / 10), linewidth = 1.5, linetype = "dashed") +
    ggplot2::geom_point(ggplot2::aes(y = K_peak), size = 2.5) +
    ggplot2::geom_point(ggplot2::aes(y = L_peak / 10), size = 2.5) +
    ggplot2::scale_colour_manual(values = arm_colors[c("ribavirin", "favipiravir", "combination")]) +
    ggplot2::labs(
      x = NULL,
      y = "Organ Injury Index",
      colour = "Arm",
      title = "B. Treatment Start Day vs Peak Organ Injury\n(K = solid, L/10 = dashed)"
    ) +
    theme_qsp() +
    ggplot2::scale_x_continuous(breaks = 1:7)

  # Panel C: Viral load AUC vs treatment day (bar plot)
  p3 <- ggplot2::ggplot(window_data, ggplot2::aes(
    x = factor(treatment_day), y = V_AUC / 1e6, fill = arm
  )) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::scale_fill_manual(values = arm_colors[c("ribavirin", "favipiravir", "combination")]) +
    ggplot2::labs(
      x = "Treatment Start Day",
      y = expression("Viral Load AUC" ~ (10^6 * " copies" * cdot ~ "day/mL")),
      fill = "Arm",
      title = "C. Treatment Start Day vs Viral Exposure (AUC)"
    ) +
    theme_qsp() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0))

  # Combine panels using grid
  gridExtra::grid.arrange(p1, p2, p3, ncol = 1, heights = c(1, 1, 1))
  ggsave_safe(outfile, width = 10, height = 14, dpi = 300)
  message(sprintf("Saved: %s", outfile))
}

#' Safe ggsave wrapper with directory creation
#'
#' @param filename Output file path
#' @param ... Passed to ggplot2::ggsave
#' @export
ggsave_safe <- function(filename, ...) {
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(filename, ...)
}

#' Plot adaptive immunity trajectories
#'
#' Creates a 4-panel figure showing CD8+ T cells, CD4+ T cells, IgM, and
#' IgG antibody trajectories by treatment arm.
#'
#' @param sim_data Data frame with columns: time, CD8_E, CD4, IgM, IgG, arm,
#'   patient_id
#' @param outfile Output file path (PNG)
#' @export
plot_adaptive_immunity <- function(sim_data,
                                    outfile = "outputs/adaptive_immunity.png") {
  arm_colors <- c(
    placebo      = "#757575",
    ribavirin    = "#1f78b4",
    favipiravir  = "#E69F00",
    combination  = "#D55E00"
  )

  # Convert to long format
  long_data <- tidyr::pivot_longer(
    sim_data,
    cols = c("CD8_E", "CD4", "IgM", "IgG"),
    names_to = "compartment",
    values_to = "value"
  )

  long_data$arm <- factor(long_data$arm,
    levels = c("placebo", "ribavirin", "favipiravir", "combination"))

  p <- ggplot2::ggplot(long_data, ggplot2::aes(
    x = time, y = value, colour = arm,
    group = interaction(arm, patient_id)
  )) +
    ggplot2::geom_line(alpha = 0.12, linewidth = 0.2) +
    ggplot2::stat_summary(
      fun = median, geom = "line", linewidth = 2.0,
      ggplot2::aes(group = arm)
    ) +
    ggplot2::facet_wrap(~compartment, scales = "free_y", ncol = 2,
                         labeller = ggplot2::as_labeller(c(
                           CD8_E = "CD8+ effector T cells (AU)",
                           CD4 = "CD4+ helper T cells (AU)",
                           IgM = "IgM antibodies (AU)",
                           IgG = "IgG antibodies (AU)"
                         ))) +
    ggplot2::labs(
      x = "Time (days)",
      y = "Value (AU)",
      colour = "Treatment arm",
      title = "Adaptive Immunity Trajectories"
    ) +
    theme_qsp() +
    ggplot2::scale_colour_manual(values = arm_colors)

  ggsave_safe(outfile, width = 10, height = 8, dpi = 300)
  message(sprintf("Saved: %s", outfile))
}

#' Generate clinical endpoints summary table
#'
#' @param trial_output Output from [run_virtual_trial()]
#' @param outfile CSV output path
#' @export
write_clinical_endpoints <- function(trial_output,
                                      outfile = "outputs/clinical_endpoints.csv") {
  ep <- trial_output$endpoints

  summary_tbl <- ep |>
    dplyr::group_by(arm) |>
    dplyr::summarise(
      n              = dplyr::n(),
      V_peak_median  = median(V_peak, na.rm = TRUE),
      V_peak_IQR     = IQR(V_peak, na.rm = TRUE),
      V_AUC_median   = median(V_AUC, na.rm = TRUE),
      V_AUC_IQR      = IQR(V_AUC, na.rm = TRUE),
      clearance_rate = mean(!is.na(time_clearance), na.rm = TRUE),
      PLT_nadir_median = median(PLT_nadir, na.rm = TRUE),
      K_peak_median  = median(K_peak, na.rm = TRUE),
      L_peak_median  = median(L_peak, na.rm = TRUE),
      dialysis_risk_mean = mean(dialysis_prob, na.rm = TRUE),
      ecmo_risk_mean = mean(ecmo_prob, na.rm = TRUE),
      mortality_risk_mean = mean(mortality_prob, na.rm = TRUE),
      dialysis_median = median(dialysis_prob, na.rm = TRUE),
      ecmo_median    = median(ecmo_prob, na.rm = TRUE),
      mortality_median = median(mortality_prob, na.rm = TRUE),
      .groups = "drop"
    )

  dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(summary_tbl, outfile, row.names = FALSE)
  message(sprintf("Saved: %s", outfile))
  summary_tbl
}

#' Generate simulation summary text
#'
#' @param trial_output Output from [run_virtual_trial()]
#' @param window_data Output from [run_window_analysis()]
#' @param outfile Text output path
#' @export
write_simulation_summary <- function(trial_output, window_data,
                                      outfile = "outputs/simulation_summary.txt") {
  ep <- trial_output$endpoints

  lines <- character()
  lines <- c(lines, "=" %s+% paste(rep("=", 59), collapse = ""))
  lines <- c(lines, "  HANTAVIRUS QSP VIRTUAL TRIAL — SIMULATION SUMMARY")
  lines <- c(lines, "=" %s+% paste(rep("=", 59), collapse = ""))
  lines <- c(lines, "")
  lines <- c(lines, sprintf("Date: %s", Sys.Date()))
  lines <- c(lines, sprintf("Patients per arm: %d",
                             length(unique(ep$patient_id))))
  lines <- c(lines, sprintf("Arms: %s",
                             paste(unique(ep$arm), collapse = ", ")))
  lines <- c(lines, "")

  # Viral load summary
  lines <- c(lines, "--- VIRAL LOAD ---")
  for (arm in unique(ep$arm)) {
    sub <- ep[ep$arm == arm, ]
    lines <- c(lines, sprintf(
      "  %-15s  Peak: %.0f (IQR: %.0f)  AUC: %.1f  Clearance: %.0f%%",
      arm,
      median(sub$V_peak, na.rm = TRUE),
      IQR(sub$V_peak, na.rm = TRUE),
      median(sub$V_AUC, na.rm = TRUE),
      100 * mean(!is.na(sub$time_clearance), na.rm = TRUE)
    ))
  }

  lines <- c(lines, "")
  lines <- c(lines, "--- CLINICAL ENDPOINTS ---")
  for (arm in unique(ep$arm)) {
    sub <- ep[ep$arm == arm, ]
    lines <- c(lines, sprintf(
      "  %-15s  Dialysis: %.1f%%  ECMO: %.1f%%  Mortality: %.1f%%",
      arm,
      100 * mean(sub$dialysis_prob, na.rm = TRUE),
      100 * mean(sub$ecmo_prob, na.rm = TRUE),
      100 * mean(sub$mortality_prob, na.rm = TRUE)
    ))
  }

  lines <- c(lines, "")
  lines <- c(lines, "--- TREATMENT WINDOW ---")
  if (!is.null(window_data) && nrow(window_data) > 0) {
    for (arm in unique(window_data$arm)) {
      sub <- window_data[window_data$arm == arm, ]
      lines <- c(lines, sprintf("  %s:", arm))
      for (day in sort(unique(sub$treatment_day))) {
        sub_day <- sub[sub$treatment_day == day, ]
        lines <- c(lines, sprintf(
          "    Day %d: V_AUC=%.1f  Mortality=%.1f%%",
          day, sub_day$V_AUC[1], 100 * sub_day$mortality_prob[1]
        ))
      }
    }
  }

  lines <- c(lines, "")
  lines <- c(lines, "String concatenation helper: %s% <- function(a, b) paste0(a, b)")
  lines <- c(lines, "=" %s+% paste(rep("=", 59), collapse = ""))

  dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
  writeLines(lines, outfile)
  message(sprintf("Saved: %s", outfile))
}
