# Clinical summary: efficacy broken down by cytokine storm offset
# Uses pure R RK4 since deSolve lsoda crashes on this machine
#
# Key question: does the pro-inflammatory / anti-inflammatory offset
# (C_anti peak time - C_pro peak time) predict mortality? Does treatment
# shrink the offset?

# ── RK4 solver with dosing events ──────────────────────────────────────

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

rk4_with_events <- function(y, times, func, parms, events = NULL) {
  n <- length(times)
  ny <- length(y)
  state_names <- names(y)
  out <- matrix(NA, nrow = n, ncol = ny + 1)
  colnames(out) <- c("time", state_names)

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

# ── Load model ─────────────────────────────────────────────────────────

library(deSolve)
source("model/hantavirus_qsp.R")
source("model/parameters.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/virtual_population.R")
source("R/simulate_trial.R")

pop <- read.csv("outputs/virtual_population.csv")
pars <- get_parameters()

# ── Configuration ──────────────────────────────────────────────────────

set.seed(42)
n_pts <- 60
sample_ids <- sample(pop$patient_id, n_pts)
dt <- 0.25
times <- seq(0, 21, by = dt)

arms <- c("placebo", "ribavirin", "favipiravir", "combination")
arm_labels <- c(
  placebo     = "Placebo",
  ribavirin   = "Ribavirin",
  favipiravir = "Favipiravir",
  combination = "Combo"
)

# ── Run simulations ────────────────────────────────────────────────────

results <- list()
for (arm in arms) {
  results[[arm]] <- data.frame(
    patient_id       = integer(),
    C_pro_peak       = numeric(),
    C_pro_peak_day   = numeric(),
    C_anti_peak      = numeric(),
    C_anti_peak_day  = numeric(),
    offset_days      = numeric(),
    P_peak           = numeric(),
    K_peak           = numeric(),
    L_peak           = numeric(),
    C_pro_AUC        = numeric(),
    K_AUC            = numeric(),
    L_AUC            = numeric(),
    V_peak           = numeric(),
    V_AUC            = numeric(),
    dialysis_prob    = numeric(),
    ecmo_prob        = numeric(),
    mortality_prob   = numeric(),
    stringsAsFactors = FALSE
  )
}

cat("Running", n_pts, "patients x 4 arms with RK4 (dt =", dt, ")...\n")

clamp_nan <- function(x) { x[!is.finite(x)] <- 0; x }

for (i in seq_along(sample_ids)) {
  pid <- sample_ids[i]
  row <- pop[pop$patient_id == pid, ]
  pars_i <- apply_vpop_to_model(row, pars)
  y0 <- build_symptom_onset_state(pars_i)
  bw <- if ("body_weight_kg" %in% names(pars_i)) pars_i$body_weight_kg else 75

  for (arm in arms) {
    # Build events
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
    sim <- as.data.frame(sim)

    # Extract cytokine peaks and offset
    C_pro_idx  <- which.max(sim$C_pro)
    C_anti_idx <- which.max(sim$C_anti)
    offset <- sim$time[C_anti_idx] - sim$time[C_pro_idx]

    # Organ injury
    K_peak <- max(sim$K, na.rm = TRUE)
    L_peak <- max(sim$L, na.rm = TRUE)
    P_peak <- max(sim$P, na.rm = TRUE)
    V_peak <- max(sim$V, na.rm = TRUE)

    # AUCs (trapezoidal)
    C_pro_AUC <- auc_trapezoidal(sim$time, clamp_nan(sim$C_pro))
    K_AUC     <- auc_trapezoidal(sim$time, clamp_nan(sim$K))
    L_AUC     <- auc_trapezoidal(sim$time, clamp_nan(sim$L))
    V_AUC     <- auc_trapezoidal(sim$time, clamp_nan(sim$V))

    # Clinical endpoints
    ep <- compute_clinical_endpoints(
      K_peak, L_peak, max(sim$C_pro, na.rm = TRUE),
      max(sim$C_pro, na.rm = TRUE), "HFRS", pars_i,
      K_auc = K_AUC, L_auc = L_AUC, C_auc = C_pro_AUC
    )

    results[[arm]][i, ] <- data.frame(
      patient_id       = pid,
      C_pro_peak       = max(sim$C_pro, na.rm = TRUE),
      C_pro_peak_day   = sim$time[C_pro_idx],
      C_anti_peak      = max(sim$C_anti, na.rm = TRUE),
      C_anti_peak_day  = sim$time[C_anti_idx],
      offset_days      = offset,
      P_peak           = P_peak,
      K_peak           = K_peak,
      L_peak           = L_peak,
      C_pro_AUC        = C_pro_AUC,
      K_AUC            = K_AUC,
      L_AUC            = L_AUC,
      V_peak           = V_peak,
      V_AUC            = V_AUC,
      dialysis_prob    = ep$dialysis_prob,
      ecmo_prob        = ep$ecmo_prob,
      mortality_prob   = ep$mortality_prob,
      stringsAsFactors = FALSE
    )
  }

  if (i %% 10 == 0) cat("  ", i, "/", n_pts, "\n")
}

# Combine all arms into one data frame
all_data <- do.call(rbind, lapply(arms, function(a) {
  df <- results[[a]]
  df$arm <- a
  df
}))

# ── Stratify by cytokine offset ────────────────────────────────────────

# Use clinically meaningful cut-points:
#   <1d: treatment-attenuated (anti-inflammatory catches up fast)
#   1-3d: intermediate
#   >3d: unopposed cytokine storm (delayed resolution)
offset_breaks <- c(-Inf, 1, 3, Inf)
offset_labs <- c("Short (<1d)", "Intermediate (1-3d)", "Long (>3d)")
all_data$offset_stratum <- cut(all_data$offset_days,
  breaks = offset_breaks,
  labels = offset_labs,
  include.lowest = TRUE)

# ── Generate summary tables ────────────────────────────────────────────

cat("\n", paste(rep("=", 68), collapse = ""), "\n")
cat("  CLINICAL SUMMARY — EFFICACY BY CYTOKINE STORM OFFSET\n")
cat(paste(rep("=", 68), collapse = ""), "\n\n")

cat(sprintf("Date: %s\n", Sys.Date()))
cat(sprintf("Patients: %d per arm\n", n_pts))
cat(sprintf("Solver: Pure R RK4, dt = %.2f days\n\n", dt))

# Overall arm-level summary
cat("─── ARM-LEVEL CYTOKINE SUMMARY ───\n\n")
cat(sprintf("%-14s %8s %8s %8s %10s %12s %12s\n",
            "Arm", "C_pro_pk", "C_anti_pk", "Offset",
            "P_peak", "Mort_prob", "Dialysis"))
cat(sprintf("%-14s %8s %8s %8s %10s %12s %12s\n",
            "────", "───────", "────────", "──────",
            "──────", "──────────", "────────"))
for (a in arms) {
  d <- results[[a]]
  cat(sprintf("%-14s %8.3f %8.4f %7.1fd %10.1f %11.1f%% %11.1f%%\n",
    arm_labels[a],
    median(d$C_pro_peak),
    median(d$C_anti_peak),
    median(d$offset_days),
    median(d$P_peak),
    median(d$mortality_prob) * 100,
    median(d$dialysis_prob) * 100))
}

# Offset distribution by arm
cat("\n─── CYTOKINE OFFSET DISTRIBUTION BY ARM ───\n\n")
for (a in arms) {
  d <- results[[a]]
  ofs <- d$offset_days[is.finite(d$offset_days)]
  cat(sprintf("  %-12s: median=%5.1fd  IQR=[%5.1f, %5.1f]  range=[%5.1f, %5.1f]\n",
    arm_labels[a],
    median(ofs), quantile(ofs, 0.25), quantile(ofs, 0.75),
    min(ofs), max(ofs)))
}

# Efficacy broken down by offset stratum (pooled across arms)
cat("\n─── EFFICACY BY OFFSET STRATUM (POOLED ARMS) ───\n\n")
offset_levels <- offset_labs
cat(sprintf("%-18s %8s %12s %12s %12s %12s\n",
            "Stratum", "N", "Mort_prob", "Dialysis", "K_peak", "P_peak"))
cat(sprintf("%-18s %8s %12s %12s %12s %12s\n",
            "──────", "──", "──────────", "────────", "──────", "──────"))
for (os in offset_levels) {
  sub <- all_data[all_data$offset_stratum == os, ]
  n <- nrow(sub)
  if (n > 0) {
    cat(sprintf("%-18s %8d %11.1f%% %11.1f%% %11.3f %11.1f\n",
      os, n,
      median(sub$mortality_prob) * 100,
      median(sub$dialysis_prob) * 100,
      median(sub$K_peak),
      median(sub$P_peak)))
  }
}

# Efficacy by arm × offset stratum
cat("\n─── MORTALITY BY ARM × OFFSET STRATUM ───\n\n")
cat(sprintf("%-18s %12s %12s %12s %12s\n",
            "Stratum", "Placebo", "Ribavirin", "Favipiravir", "Combo"))
cat(sprintf("%-18s %12s %12s %12s %12s\n",
            "──────", "───────", "─────────", "───────────", "─────"))
for (os in offset_levels) {
  vals <- sapply(arms, function(a) {
    sub <- all_data[all_data$offset_stratum == os & all_data$arm == a, ]
    if (nrow(sub) > 0) median(sub$mortality_prob) * 100 else NA
  })
  cat(sprintf("%-18s %11.1f%% %11.1f%% %11.1f%% %11.1f%%\n",
    os, vals[1], vals[2], vals[3], vals[4]))
}

# Dialysis by arm × offset stratum
cat("\n─── DIALYSIS BY ARM × OFFSET STRATUM ───\n\n")
cat(sprintf("%-18s %12s %12s %12s %12s\n",
            "Stratum", "Placebo", "Ribavirin", "Favipiravir", "Combo"))
cat(sprintf("%-18s %12s %12s %12s %12s\n",
            "──────", "───────", "─────────", "───────────", "─────"))
for (os in offset_levels) {
  vals <- sapply(arms, function(a) {
    sub <- all_data[all_data$offset_stratum == os & all_data$arm == a, ]
    if (nrow(sub) > 0) median(sub$dialysis_prob) * 100 else NA
  })
  cat(sprintf("%-18s %11.1f%% %11.1f%% %11.1f%% %11.1f%%\n",
    os, vals[1], vals[2], vals[3], vals[4]))
}

# Correlation: offset vs mortality
cat("\n─── CORRELATION: OFFSET vs MORTALITY ───\n\n")
for (a in arms) {
  d <- results[[a]]
  rho <- cor(d$offset_days, d$mortality_prob, method = "spearman",
             use = "complete.obs")
  cat(sprintf("  %-12s: Spearman rho = %+.3f\n", arm_labels[a], rho))
}

# Treatment effect on offset
cat("\n─── TREATMENT EFFECT ON OFFSET ───\n\n")
placebo_offsets <- results$placebo$offset_days
for (a in c("ribavirin", "favipiravir", "combination")) {
  delta <- median(results[[a]]$offset_days) - median(placebo_offsets)
  pct_change <- delta / median(placebo_offsets) * 100
  cat(sprintf("  %-12s: delta offset = %+.1f days (%+.0f%% vs placebo)\n",
    arm_labels[a], delta, pct_change))
}

cat("\n", paste(rep("=", 68), collapse = ""), "\n")
cat("  END OF SUMMARY\n")
cat(paste(rep("=", 68), collapse = ""), "\n")

# ── Save to file ──────────────────────────────────────────────────────
sink("outputs/clinical_summary_cytokine_offset.txt")
cat("# Clinical Summary — Efficacy by Cytokine Storm Offset\n")
cat(paste(rep("=", 68), collapse = ""), "\n\n")
cat(sprintf("Date: %s\n", Sys.Date()))
cat(sprintf("Patients: %d per arm\n", n_pts))
cat(sprintf("Solver: Pure R RK4, dt = %.2f days\n\n", dt))

cat("## Arm-Level Cytokine Summary\n\n")
cat(sprintf("%-14s %8s %8s %8s %10s %12s %12s\n",
            "Arm", "C_pro_pk", "C_anti_pk", "Offset",
            "P_peak", "Mort_prob", "Dialysis"))
cat(sprintf("%-14s %8s %8s %8s %10s %12s %12s\n",
            "────", "───────", "────────", "──────",
            "──────", "──────────", "────────"))
for (a in arms) {
  d <- results[[a]]
  cat(sprintf("%-14s %8.3f %8.4f %7.1fd %10.1f %11.1f%% %11.1f%%\n",
    arm_labels[a], median(d$C_pro_peak), median(d$C_anti_peak),
    median(d$offset_days), median(d$P_peak),
    median(d$mortality_prob) * 100, median(d$dialysis_prob) * 100))
}

cat("\n## Offset Distribution by Arm\n\n")
for (a in arms) {
  d <- results[[a]]
  ofs <- d$offset_days
  cat(sprintf("* %-12s: median=%5.1fd  IQR=[%5.1f, %5.1f]  range=[%5.1f, %5.1f]\n",
    arm_labels[a], median(ofs), quantile(ofs, 0.25), quantile(ofs, 0.75),
    min(ofs), max(ofs)))
}

cat("\n## Efficacy by Offset Stratum (Pooled Arms)\n\n")
for (os in offset_levels) {
  sub <- all_data[all_data$offset_stratum == os, ]
  n <- nrow(sub)
  if (n > 0) {
    cat(sprintf("* %-20s (N=%d): Mort=%5.1f%%, Dialysis=%5.1f%%, P_peak=%6.1f, K_peak=%5.2f\n",
      os, n, median(sub$mortality_prob) * 100, median(sub$dialysis_prob) * 100,
      median(sub$P_peak), median(sub$K_peak)))
  }
}

cat("\n## Mortality by Arm × Offset Stratum\n\n")
cat(sprintf("%-20s %10s %10s %10s %10s\n", "Stratum", "Placebo", "RBV", "FAV", "Combo"))
for (os in offset_levels) {
  vals <- sapply(arms, function(a) {
    sub <- all_data[all_data$offset_stratum == os & all_data$arm == a, ]
    if (nrow(sub) > 0) median(sub$mortality_prob) * 100 else NA
  })
  if (!all(is.na(vals))) {
    cat(sprintf("%-20s %9.1f%% %9.1f%% %9.1f%% %9.1f%%\n",
      os, vals[1], vals[2], vals[3], vals[4]))
  }
}

cat("\n## Treatment Effect on Offset vs Placebo\n\n")
for (a in c("ribavirin", "favipiravir", "combination")) {
  d <- median(results[[a]]$offset_days) - median(results$placebo$offset_days)
  p <- d / median(results$placebo$offset_days) * 100
  cat(sprintf("* %-12s: delta offset = %+.1f days (%+.0f%% vs placebo)\n",
    arm_labels[a], d, p))
}

cat("\n## Correlation: Offset vs Mortality (Spearman)\n\n")
for (a in arms) {
  rho <- cor(results[[a]]$offset_days, results[[a]]$mortality_prob,
             method = "spearman", use = "complete.obs")
  cat(sprintf("* %-12s: rho = %+.3f\n", arm_labels[a], rho))
}
sink()
cat("Summary saved to outputs/clinical_summary_cytokine_offset.txt\n")
