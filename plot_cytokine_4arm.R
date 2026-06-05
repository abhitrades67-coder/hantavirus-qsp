# 4-arm cytokine trajectory comparison
# Uses pure R RK4 since deSolve lsoda crashes on this machine

rk4_pure <- function(y, times, func, parms) {
  n <- length(times)
  ny <- length(y)
  state_names <- names(y)
  out <- matrix(NA, nrow = n, ncol = ny + 1)
  colnames(out) <- c("time", state_names)
  out[1, ] <- c(times[1], y)
  for (i in 2:n) {
    t <- times[i - 1]
    h <- times[i] - times[i - 1]
    y_curr <- setNames(out[i - 1, -1], state_names)
    k1 <- func(t, y_curr, parms)[[1]]
    k2 <- func(t + h / 2, y_curr + h / 2 * k1, parms)[[1]]
    k3 <- func(t + h / 2, y_curr + h / 2 * k2, parms)[[1]]
    k4 <- func(t + h, y_curr + h * k3, parms)[[1]]
    y_next <- y_curr + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    names(y_next) <- state_names
    out[i, ] <- c(times[i], y_next)
  }
  out
}

# RK4 with dosing events (IV bolus = add to compartment)
rk4_with_events <- function(y, times, func, parms, events = NULL) {
  n <- length(times)
  ny <- length(y)
  state_names <- names(y)
  out <- matrix(NA, nrow = n, ncol = ny + 1)
  colnames(out) <- c("time", state_names)

  # Apply any events at time 0
  y0 <- y
  if (!is.null(events) && nrow(events) > 0) {
    t0_events <- events[events$time == times[1], ]
    for (r in seq_len(nrow(t0_events))) {
      if (t0_events$var[r] %in% names(y0)) {
        y0[t0_events$var[r]] <- y0[t0_events$var[r]] + t0_events$value[r]
      }
    }
  }
  out[1, ] <- c(times[1], y0)

  for (i in 2:n) {
    t <- times[i - 1]
    h <- times[i] - times[i - 1]
    y_curr <- setNames(out[i - 1, -1], state_names)

    # Apply events that fall within (t, t+h] at the step start
    if (!is.null(events) && nrow(events) > 0) {
      step_events <- events[events$time > t & events$time <= t + h, ]
      if (nrow(step_events) > 0) {
        for (r in seq_len(nrow(step_events))) {
          if (step_events$var[r] %in% names(y_curr)) {
            y_curr[step_events$var[r]] <- y_curr[step_events$var[r]] +
                                            step_events$value[r]
          }
        }
      }
    }

    k1 <- func(t, y_curr, parms)[[1]]
    k2 <- func(t + h / 2, y_curr + h / 2 * k1, parms)[[1]]
    k3 <- func(t + h / 2, y_curr + h / 2 * k2, parms)[[1]]
    k4 <- func(t + h, y_curr + h * k3, parms)[[1]]
    y_next <- y_curr + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    names(y_next) <- state_names
    out[i, ] <- c(times[i], y_next)
  }
  out
}

library(deSolve)
source("model/hantavirus_qsp.R")
source("model/parameters.R")
source("R/pk_models.R")
source("R/virtual_population.R")
source("R/simulate_trial.R")

pop <- read.csv("outputs/virtual_population.csv")
pars <- get_parameters()

set.seed(42)
n_pts <- 20
sample_ids <- sample(pop$patient_id, n_pts)
dt <- 0.25
times <- seq(0, 21, by = dt)

arms <- c("placebo", "ribavirin", "favipiravir", "combination")
arm_colors <- c(
  placebo     = "grey40",
  ribavirin   = "#D55E00",
  favipiravir = "#0072B2",
  combination = "#009E73"
)
arm_labels <- c(
  placebo     = "Placebo",
  ribavirin   = "Ribavirin",
  favipiravir = "Favipiravir",
  combination = "Combo (RBV+FAV)"
)

# Pre-allocate: list of arms, each a list of C_pro, C_anti, P matrices
results <- setNames(vector("list", length(arms)), arms)
for (a in arms) {
  results[[a]] <- list(
    C_pro  = matrix(NA, nrow = n_pts, ncol = length(times)),
    C_anti = matrix(NA, nrow = n_pts, ncol = length(times)),
    P      = matrix(NA, nrow = n_pts, ncol = length(times))
  )
}

cat("Running", n_pts, "patients x 4 arms with RK4 (dt =", dt, ")...\n")

for (i in seq_along(sample_ids)) {
  pid <- sample_ids[i]
  row <- pop[pop$patient_id == pid, ]
  pars_i <- apply_vpop_to_model(row, pars)
  y0 <- build_symptom_onset_state(pars_i)
  bw <- if ("body_weight_kg" %in% names(pars_i)) pars_i$body_weight_kg else 75

  # Build event lists for each arm
  for (arm in arms) {
    if (arm == "placebo") {
      evts <- NULL
    } else if (arm == "ribavirin") {
      sched <- ribavirin_dosing_schedule(body_weight_kg = bw, start_day = 3,
                                          duration_days = 10)
      evts <- build_deSolve_events(sched, "C_RBV", body_weight_kg = bw,
                                    Vd = pars_i$Vd_RBV)
    } else if (arm == "favipiravir") {
      sched <- favipiravir_dosing_schedule(regimen = "standard", start_day = 3)
      evts <- build_deSolve_events(sched, "C_FAV_gut", body_weight_kg = bw,
                                    Vd = pars_i$Vd_FAV)
    } else if (arm == "combination") {
      sched_rbv <- ribavirin_dosing_schedule(body_weight_kg = bw, start_day = 3,
                                              duration_days = 10)
      sched_fav <- favipiravir_dosing_schedule(regimen = "standard", start_day = 3)
      evts_rbv <- build_deSolve_events(sched_rbv, "C_RBV", body_weight_kg = bw,
                                        Vd = pars_i$Vd_RBV)
      evts_fav <- build_deSolve_events(sched_fav, "C_FAV_gut", body_weight_kg = bw,
                                        Vd = pars_i$Vd_FAV)
      evts <- rbind(evts_rbv, evts_fav)
      evts <- evts[order(evts$time), ]
      rownames(evts) <- NULL
    }

    sim <- rk4_with_events(y0, times, hantavirus_qsp_ode, pars_i, evts)
    results[[arm]]$C_pro[i, ]  <- sim[, "C_pro"]
    results[[arm]]$C_anti[i, ] <- sim[, "C_anti"]
    results[[arm]]$P[i, ]      <- sim[, "P"]
  }

  if (i %% 5 == 0) cat("  ", i, "/", n_pts, "\n")
}

# Compute medians and IQRs
stats <- list()
for (arm in arms) {
  stats[[arm]] <- list(
    med_C_pro  = apply(results[[arm]]$C_pro, 2, median, na.rm = TRUE),
    q25_C_pro  = apply(results[[arm]]$C_pro, 2, quantile, 0.25, na.rm = TRUE),
    q75_C_pro  = apply(results[[arm]]$C_pro, 2, quantile, 0.75, na.rm = TRUE),
    med_C_anti = apply(results[[arm]]$C_anti, 2, median, na.rm = TRUE),
    q25_C_anti = apply(results[[arm]]$C_anti, 2, quantile, 0.25, na.rm = TRUE),
    q75_C_anti = apply(results[[arm]]$C_anti, 2, quantile, 0.75, na.rm = TRUE),
    med_P      = apply(results[[arm]]$P, 2, median, na.rm = TRUE),
    q25_P      = apply(results[[arm]]$P, 2, quantile, 0.25, na.rm = TRUE),
    q75_P      = apply(results[[arm]]$P, 2, quantile, 0.75, na.rm = TRUE)
  )
}

# Compute overall y-limits (shared across arms in each panel)
global_max_C_pro  <- max(sapply(arms, function(a) max(stats[[a]]$q75_C_pro)))
global_max_C_anti <- max(sapply(arms, function(a) max(stats[[a]]$q75_C_anti)))
global_max_P      <- max(sapply(arms, function(a) max(stats[[a]]$q75_P)))

png("outputs/cytokine_trajectories_4arm.png", width = 1800, height = 1200, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), mgp = c(2.5, 0.8, 0))

# Panel A: C_pro — all 4 arms
ylimA <- c(0, global_max_C_pro * 1.15)
plot(times, stats$placebo$med_C_pro, type = "n",
     xlab = "Days post symptom onset", ylab = "C_pro (AU)",
     main = "Pro-inflammatory (C_pro)", ylim = ylimA, cex.main = 1.1)
for (a in arms) {
  s <- stats[[a]]
  col <- arm_colors[a]
  polygon(c(times, rev(times)), c(s$q25_C_pro, rev(s$q75_C_pro)),
          col = adjustcolor(col, alpha.f = 0.15), border = NA)
  lines(times, s$med_C_pro, col = col, lwd = 2.2)
}
abline(v = 3, lty = 3, col = "grey50")
text(3.2, ylimA[2] * 0.97, "Tx start", col = "grey50", cex = 0.7)
legend("topright", arm_labels, col = arm_colors, lwd = 2.2, bty = "n", cex = 0.75)

# Panel B: C_anti — all 4 arms
ylimB <- c(0, global_max_C_anti * 1.15)
plot(times, stats$placebo$med_C_anti, type = "n",
     xlab = "Days post symptom onset", ylab = "C_anti (AU)",
     main = "Anti-inflammatory (C_anti)", ylim = ylimB, cex.main = 1.1)
for (a in arms) {
  s <- stats[[a]]
  col <- arm_colors[a]
  polygon(c(times, rev(times)), c(s$q25_C_anti, rev(s$q75_C_anti)),
          col = adjustcolor(col, alpha.f = 0.15), border = NA)
  lines(times, s$med_C_anti, col = col, lwd = 2.2)
}
abline(v = 3, lty = 3, col = "grey50")

# Panel C: Normalized C_pro vs C_anti (just placebo for clarity of offset)
plot(times, stats$placebo$med_C_pro / max(stats$placebo$med_C_pro), type = "n",
     xlab = "Days post symptom onset", ylab = "Normalized (0-1)",
     main = "Normalized C_pro vs C_anti", cex.main = 1.1)
for (a in arms) {
  s <- stats[[a]]
  col <- arm_colors[a]
  np <- s$med_C_pro / max(s$med_C_pro)
  na <- s$med_C_anti / max(s$med_C_anti)
  lines(times, np, col = col, lwd = 2.2, lty = 1)
  lines(times, na, col = col, lwd = 2.2, lty = 2)
}
legend("topright",
       c("C_pro (solid)", "C_anti (dashed)"),
       lty = c(1, 2), lwd = 1.5, bty = "n", cex = 0.7)
legend("bottomright", arm_labels, col = arm_colors, lwd = 2, bty = "n", cex = 0.65)

# Panel D: Permeability — all 4 arms
ylimD <- c(0, global_max_P * 1.15)
plot(times, stats$placebo$med_P, type = "n",
     xlab = "Days post symptom onset", ylab = "P (AU)",
     main = "Vascular permeability (P)", ylim = ylimD, cex.main = 1.1)
for (a in arms) {
  s <- stats[[a]]
  col <- arm_colors[a]
  polygon(c(times, rev(times)), c(s$q25_P, rev(s$q75_P)),
          col = adjustcolor(col, alpha.f = 0.15), border = NA)
  lines(times, s$med_P, col = col, lwd = 2.2)
}
abline(v = 3, lty = 3, col = "grey50")

dev.off()

cat("\nPlot saved to outputs/cytokine_trajectories_4arm.png\n")

cat("\nPeak C_pro by arm (median of", n_pts, "patients):\n")
for (a in arms) {
  peak_day <- times[which.max(stats[[a]]$med_C_pro)]
  peak_val <- max(stats[[a]]$med_C_pro)
  cat(sprintf("  %-12s  day %5.1f  value = %8.4f\n", a, peak_day, peak_val))
}

cat("\nPeak C_anti by arm:\n")
for (a in arms) {
  peak_day <- times[which.max(stats[[a]]$med_C_anti)]
  peak_val <- max(stats[[a]]$med_C_anti)
  cat(sprintf("  %-12s  day %5.1f  value = %8.4f\n", a, peak_day, peak_val))
}

cat("\nPeak permeability by arm:\n")
for (a in arms) {
  peak_day <- times[which.max(stats[[a]]$med_P)]
  peak_val <- max(stats[[a]]$med_P)
  cat(sprintf("  %-12s  day %5.1f  value = %8.2f\n", a, peak_day, peak_val))
}

cat("\nPro-anti offset (C_anti peak day - C_pro peak day):\n")
for (a in arms) {
  off <- times[which.max(stats[[a]]$med_C_anti)] - times[which.max(stats[[a]]$med_C_pro)]
  cat(sprintf("  %-12s  %5.1f days\n", a, off))
}
