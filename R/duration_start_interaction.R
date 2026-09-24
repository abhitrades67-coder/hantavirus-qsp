#' Interaction Between Treatment Start Day and Course Duration
#'
#' The treatment-window analysis (R/simulate_trial.R, run_window_analysis) holds
#' the course length fixed (ribavirin 10 d, favipiravir 15 d) and varies only the
#' start day. Under that design, day-1 favipiravir monotherapy looks WORSE than
#' day-2 -- the "favipiravir day-1 paradox". This script shows the paradox is an
#' artefact of holding duration fixed.
#'
#' Mechanism: early suppression prevents the virus from consuming the
#' target-cell pool (61.9% of T_0 survives to the end of a day-1 course, vs
#' 23.0% for day-2 and ~0% untreated). Because the shipped humoral module
#' contributes only ~0.1-0.4% of viral clearance, nothing terminates the
#' infection when dosing stops, so the virus rebounds into the spared pool and
#' produces a full-magnitude second peak AFTER the drug is withdrawn. The
#' resulting organ-injury peak lands outside the 21-day observation window, so
#' the window analysis both mis-attributes and understates it.
#'
#' Prediction, tested here: extending the course past the viral course abolishes
#' the rebound and restores the expected "earlier is better" ordering.
#'
#' Runs on a 42-day horizon (the rebound peak for a day-1 start falls at
#' ~day 25, so a 21-day horizon truncates it). Data cached to
#' outputs/duration_start_data.csv.

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
cache_path <- "outputs/duration_start_data.csv"

T_END       <- 42
START_DAYS  <- 1:4
DURATIONS   <- c(10, 15, 21, 28, 35)
ARMS        <- c("ribavirin", "favipiravir", "combination")

#' Simulate one patient with an explicit course duration
#'
#' Mirrors [simulate_patient()] but threads `duration_days` into both dosing
#' schedules instead of using their defaults.
#'
#' @param pars Parameter list
#' @param arm One of placebo, ribavirin, favipiravir, combination
#' @param t_start Treatment start day post-symptom-onset
#' @param duration_days Course length of the arm's longest-running drug. For
#'   ribavirin monotherapy this is the ribavirin course; for favipiravir
#'   monotherapy and for combination it is the favipiravir course, and in
#'   combination ribavirin runs 5 days fewer, as in the main trial design
#'   (favipiravir 15 d / ribavirin 10 d), floored at 5 days. Defined this way so
#'   that a column labelled "n days" means the same thing in every arm.
#' @param t_end Simulation end day
#' @param dt Output time step
#' @return Data frame of the trajectory, or NULL if the solver did not complete
simulate_patient_duration <- function(pars, arm, t_start, duration_days,
                                       t_end = T_END, dt = 0.1) {
  y0    <- build_symptom_onset_state(pars)
  times <- seq(0, t_end, by = dt)
  bw    <- if ("body_weight_kg" %in% names(pars)) pars$body_weight_kg else 75

  events <- data.frame(time = numeric(), var = character(), value = numeric(),
                       method = character(), stringsAsFactors = FALSE)
  if (arm %in% c("ribavirin", "combination")) {
    rbv_days <- if (arm == "ribavirin") {
      duration_days
    } else {
      max(duration_days - 5, 5)
    }
    events <- rbind(events, build_deSolve_events(
      ribavirin_dosing_schedule(body_weight_kg = bw, start_day = t_start,
                                duration_days = rbv_days),
      "C_RBV", body_weight_kg = bw, Vd = pars$Vd_RBV))
  }
  if (arm %in% c("favipiravir", "combination")) {
    events <- rbind(events, build_deSolve_events(
      favipiravir_dosing_schedule(regimen = "standard", start_day = t_start,
                                  duration_days = duration_days),
      "C_FAV_gut", body_weight_kg = bw, Vd = pars$Vd_FAV))
  }
  if (nrow(events) > 0) {
    events <- events[events$time >= 0 & events$time <= t_end, ]
    events <- events[order(events$time), ]
    rownames(events) <- NULL
    if (nrow(events) > 0) times <- sort(unique(c(times, events$time)))
  }

  args <- list(y = y0, times = times, func = hantavirus_qsp_ode, parms = pars,
               atol = 1e-10, rtol = 1e-8, maxsteps = 100000, hmax = 0.5)
  if (nrow(events) > 0) args$events <- list(data = events)

  out <- tryCatch(as.data.frame(do.call(deSolve::lsoda, args)),
                  error = function(e) NULL)
  # Same truncation guard as the other simulators: an incomplete trajectory
  # summarised with max()/min() is indistinguishable from a mild one.
  if (is.null(out) || nrow(out) < length(times)) return(NULL)
  if (any(!is.finite(as.matrix(out)))) return(NULL)
  out
}

if (file.exists(cache_path)) {
  cat("Loading cached duration/start-day data...\n")
  results <- read.csv(cache_path, stringsAsFactors = FALSE)
} else {
  cat("Running duration x start-day analysis (42-day horizon)...\n")
  rows <- list()

  sim_pl <- simulate_patient_duration(pars, "placebo", 0, 15)
  stopifnot(!is.null(sim_pl))
  ep_pl <- extract_endpoints_from_simulation(sim_pl, pars, syndrome = "HFRS")
  cat(sprintf("  placebo: mortality %.3f%%, K_peak %.2f at day %.1f\n",
              100 * ep_pl$mortality_prob, max(sim_pl$K),
              sim_pl$time[which.max(sim_pl$K)]))

  for (arm in ARMS) {
    for (sd in START_DAYS) {
      for (dur in DURATIONS) {
        sim <- simulate_patient_duration(pars, arm, sd, dur)
        if (is.null(sim)) {
          warning(sprintf("solver did not complete: %s, start %d, duration %d",
                          arm, sd, dur))
          rows[[length(rows) + 1]] <- data.frame(
            arm = arm, start_day = sd, duration_days = dur,
            V_peak = NA_real_, t_V_peak = NA_real_, K_peak = NA_real_,
            t_K_peak = NA_real_, L_peak = NA_real_, T_end_frac = NA_real_,
            mortality_prob = NA_real_, rebound = NA)
          next
        }
        ep <- extract_endpoints_from_simulation(sim, pars, syndrome = "HFRS")
        # "Rebound" = the viral maximum occurs after dosing stops, i.e. the
        # second wave exceeds the on-treatment peak.
        t_dose_end <- sd + dur
        rows[[length(rows) + 1]] <- data.frame(
          arm            = arm,
          start_day      = sd,
          duration_days  = dur,
          V_peak         = max(sim$V),
          t_V_peak       = sim$time[which.max(sim$V)],
          K_peak         = max(sim$K),
          t_K_peak       = sim$time[which.max(sim$K)],
          L_peak         = max(sim$L),
          T_end_frac     = sim$T[which.min(abs(sim$time - t_dose_end))] / pars$T_0,
          mortality_prob = ep$mortality_prob,
          rebound        = sim$time[which.max(sim$V)] > t_dose_end - 1)
      }
      cat(sprintf("  %s start day %d: done\n", arm, sd))
    }
  }

  results <- do.call(rbind, rows)
  results$placebo_mortality <- ep_pl$mortality_prob
  results$placebo_K_peak    <- max(sim_pl$K)
  results$RRR <- 100 * (results$placebo_mortality - results$mortality_prob) /
    results$placebo_mortality
  write.csv(results, cache_path, row.names = FALSE)
  cat(sprintf("  Cached: %s\n", cache_path))
}

# --- Report ------------------------------------------------------------------
cat("\n=== Mortality (%) by start day and course duration, 42-day horizon ===\n")
for (arm in ARMS) {
  cat(sprintf("\n  %s (placebo %.2f%%)\n", arm,
              100 * results$placebo_mortality[1]))
  sub <- results[results$arm == arm, ]
  wide <- reshape(sub[, c("start_day", "duration_days", "mortality_prob")],
                  idvar = "start_day", timevar = "duration_days",
                  direction = "wide")
  names(wide) <- sub("mortality_prob.", "d", names(wide), fixed = TRUE)
  wide[, -1] <- round(100 * wide[, -1], 3)
  print(wide, row.names = FALSE)
}

cat("\n=== Does the earliest start win, at each duration? ===\n")
for (arm in ARMS) {
  for (dur in DURATIONS) {
    sub <- results[results$arm == arm & results$duration_days == dur, ]
    sub <- sub[order(sub$start_day), ]
    best <- sub$start_day[which.min(sub$mortality_prob)]
    cat(sprintf("  %-12s %2d-day course: best start = day %d  (day1 %.3f%%, day2 %.3f%%) %s\n",
                arm, dur, best, 100 * sub$mortality_prob[1],
                100 * sub$mortality_prob[2],
                if (best == 1) "" else "<- day 1 NOT optimal"))
  }
}

# --- Figure ------------------------------------------------------------------
plot_dat <- results
plot_dat$arm <- factor(plot_dat$arm, levels = ARMS,
                       labels = c("Ribavirin", "Favipiravir", "Combination"))
plot_dat$start_lbl <- factor(sprintf("Day %d", plot_dat$start_day))

# Colourblind-safe (Okabe-Ito), matching the rest of the figure set
start_cols <- c("Day 1" = "#D55E00", "Day 2" = "#E69F00",
                "Day 3" = "#0072B2", "Day 4" = "#56B4E9")

p1 <- ggplot2::ggplot(plot_dat, ggplot2::aes(
        x = duration_days, y = 100 * mortality_prob,
        colour = start_lbl, group = start_lbl)) +
  ggplot2::geom_hline(ggplot2::aes(yintercept = 100 * placebo_mortality),
                      linetype = "dashed", colour = "grey40", linewidth = 0.6) +
  ggplot2::geom_line(linewidth = 1.2) +
  ggplot2::geom_point(size = 2.4) +
  ggplot2::facet_wrap(~arm, nrow = 1) +
  ggplot2::scale_colour_manual(values = start_cols) +
  ggplot2::scale_x_continuous(breaks = DURATIONS) +
  ggplot2::labs(
    x = "Course duration (days)", y = "Mortality (%)",
    colour = "Treatment start",
    title = "(a) Mortality depends jointly on when treatment starts and how long it runs",
    subtitle = "Dashed line = placebo. 42-day horizon. Earlier starts need longer courses.") +
  theme_qsp()

p2 <- ggplot2::ggplot(plot_dat, ggplot2::aes(
        x = duration_days, y = t_K_peak,
        colour = start_lbl, group = start_lbl)) +
  ggplot2::geom_hline(yintercept = 21, linetype = "dotted",
                      colour = "grey30", linewidth = 0.6) +
  ggplot2::geom_line(linewidth = 1.2) +
  ggplot2::geom_point(size = 2.4) +
  ggplot2::facet_wrap(~arm, nrow = 1) +
  ggplot2::scale_colour_manual(values = start_cols) +
  ggplot2::scale_x_continuous(breaks = DURATIONS) +
  ggplot2::labs(
    x = "Course duration (days)", y = "Day of peak renal injury",
    colour = "Treatment start",
    title = "(b) The injury peak tracks the end of dosing",
    subtitle = "Dotted line = the 21-day observation window used in the main analysis") +
  theme_qsp()

g <- gridExtra::arrangeGrob(p1, p2, ncol = 1)
ggsave_safe("outputs/duration_start_interaction.png", plot = g,
            width = 12, height = 10, dpi = 300)
cat("\nSaved: outputs/duration_start_interaction.png\n")
