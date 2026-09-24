# Treatment duration optimization: how long should therapy continue so
# adaptive immunity can take over viral clearance after drug withdrawal?
#
# Simulates ribavirin, favipiravir, and combination arms with treatment
# durations of 5, 7, 10, 14, 21, and 28 days. For each, observes whether
# viral load rebounds after treatment stops and whether adaptive immunity
# (CD8, IgG) has reached a protective threshold.

source("model/parameters.R")
source("model/hantavirus_qsp.R")
source("R/virtual_population.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/simulate_trial.R")

pop <- generate_virtual_population(N = 200, seed = 42)
pars <- get_parameters()

# Treatment durations to test (days)
durations <- c(5, 7, 10, 14, 21, 28)

# Use representative patient (median V0, average immune strength)
row <- pop[10, ]  # representative patient; syndrome and V0 are printed below
pars_i <- apply_vpop_to_model(row, pars)

cat(sprintf("Representative patient: %s, V0=%.0f, immune_str=%.2f\n\n",
    row$syndrome, row$V0_copies_mL, row$immune_strength))

results <- data.frame(
  arm         = character(),
  duration    = integer(),
  V_at_stop   = numeric(),
  V_nadir     = numeric(),
  V_day60     = numeric(),
  V_rebound   = logical(),
  CD8_at_stop = numeric(),
  CD8_max     = numeric(),
  IgG_at_stop = numeric(),
  IgG_max     = numeric(),
  mortality   = numeric(),
  stringsAsFactors = FALSE
)

for (arm in c("ribavirin", "favipiravir", "combination")) {
  cat(sprintf("=== %s ===\n", toupper(arm)))
  
  for (dur in durations) {
    cat(sprintf("  Duration: %d days... ", dur))
    
    # Simulate out to day 60
    sim <- tryCatch({
      # Build dosing events manually for custom duration
      source("model/hantavirus_qsp.R")
      source("model/parameters.R")
      source("R/pk_models.R")
      source("R/pd_models.R")
      
      # Time origin is SYMPTOM ONSET, matching every other therapeutic script.
      # This previously used build_initial_state(), the INFECTION-time state,
      # while still hardcoding start_day = 3; with symptom_onset_day = 5 that
      # dosed 2 days BEFORE symptom onset, making this an early-prophylaxis
      # study mislabelled as a treatment-duration study.
      # NOTE: endpoints here use a 60-day window and are therefore NOT
      # comparable to the 21-day virtual-trial numbers.
      times <- seq(0, 60, by = 0.1)
      y0 <- build_symptom_onset_state(pars_i)
      
      all_events <- data.frame(
        time = numeric(), var = character(),
        value = numeric(), method = character(),
        stringsAsFactors = FALSE
      )
      
      if (arm %in% c("ribavirin", "combination")) {
        rbv_sched <- ribavirin_dosing_schedule(
          body_weight_kg = 75,
          start_day = 3,
          duration_days = dur
        )
        rbv_events <- build_deSolve_events(
          rbv_sched, "C_RBV",
          body_weight_kg = 75,
          Vd = pars_i$Vd_RBV
        )
        # Filter to treatment duration
        rbv_events <- rbv_events[rbv_events$time <= (3 + dur), ]
        all_events <- rbind(all_events, rbv_events)
      }
      
      if (arm %in% c("favipiravir", "combination")) {
        fav_sched <- favipiravir_dosing_schedule(
          regimen = "standard",
          start_day = 3,
          duration_days = dur
        )
        fav_events <- build_deSolve_events(
          fav_sched, "C_FAV_gut",
          body_weight_kg = 75,
          Vd = pars_i$Vd_FAV
        )
        # Filter to treatment duration
        fav_events <- fav_events[fav_events$time <= (3 + dur), ]
        all_events <- rbind(all_events, fav_events)
      }
      
      if (nrow(all_events) > 0) {
        out <- deSolve::lsoda(
          y = y0, times = times, func = hantavirus_qsp_ode,
          parms = pars_i, events = list(data = all_events),
          atol = 1e-10, rtol = 1e-8
        )
      } else {
        out <- deSolve::lsoda(
          y = y0, times = times, func = hantavirus_qsp_ode,
          parms = pars_i, atol = 1e-10, rtol = 1e-8
        )
      }
      as.data.frame(out)
    }, error = function(e) NULL)
    
    if (is.null(sim)) {
      cat("FAILED\n")
      next
    }
    
    t_stop <- 3 + dur
    idx_stop <- which.min(abs(sim$time - t_stop))
    idx_60   <- which.min(abs(sim$time - 60))
    
    V_stop <- sim$V[idx_stop]
    V_nadir <- min(sim$V[sim$time >= 3 & sim$time <= t_stop])
    V_60 <- sim$V[idx_60]
    
    # Rebound: V at day 60 > V at stop by more than 2x
    rebound <- V_60 > V_stop * 2
    
    CD8_stop <- sim$CD8_E[idx_stop]
    CD8_max  <- max(sim$CD8_E)
    IgG_stop <- sim$IgG[idx_stop]
    IgG_max  <- max(sim$IgG)
    
    # Mortality endpoint
    ep <- extract_endpoints_from_simulation(sim, pars_i, syndrome = row$syndrome)
    
    results <- rbind(results, data.frame(
      arm         = arm,
      duration    = dur,
      V_at_stop   = V_stop,
      V_nadir     = V_nadir,
      V_day60     = V_60,
      V_rebound   = rebound,
      CD8_at_stop = CD8_stop,
      CD8_max     = CD8_max,
      IgG_at_stop = IgG_stop,
      IgG_max     = IgG_max,
      mortality   = ep$mortality_prob,
      stringsAsFactors = FALSE
    ))
    
    cat(sprintf("V_stop=%.2e, V_day60=%.2e, rebound=%s, CD8=%.1f, IgG=%.2f, mort=%.1f%%\n",
        V_stop, V_60, rebound, CD8_stop, IgG_stop, ep$mortality_prob * 100))
  }
  cat("\n")
}

# Summary table
cat("\n========== SUMMARY TABLE ==========\n\n")
cat(sprintf("%-15s %-8s %-10s %-10s %-8s %-6s %-6s %-8s\n",
    "Arm", "Dur(d)", "V_stop", "V_day60", "Rebound", "CD8*", "IgG*", "Mort%"))
cat(paste(rep("-", 80), collapse = ""), "\n")

for (arm in c("ribavirin", "favipiravir", "combination")) {
  r <- results[results$arm == arm, ]
  for (i in 1:nrow(r)) {
    cat(sprintf("%-15s %-8d %-10.2e %-10.2e %-8s %-6.1f %-6.2f %-8.1f\n",
        arm, r$duration[i], r$V_at_stop[i], r$V_day60[i],
        ifelse(r$V_rebound[i], "YES", "NO"),
        r$CD8_at_stop[i], r$IgG_at_stop[i], r$mortality[i] * 100))
  }
}

cat("\nKey thresholds:\n")
cat("  - Adaptive takeover: CD8 > 10 AU AND IgG > 1 AU at drug stop\n")
cat("  - No rebound: V_day60 <= V_stop * 2\n")
cat("  - Clinical success: mortality < 2%\n")

# Identify optimal durations
cat("\n========== OPTIMAL DURATIONS ==========\n\n")
for (arm in c("ribavirin", "favipiravir", "combination")) {
  r <- results[results$arm == arm, ]
  # Find shortest duration with no rebound and adaptive takeover
  ok <- r$CD8_at_stop > 10 & r$IgG_at_stop > 1 & !r$V_rebound
  if (any(ok)) {
    best <- r[which(ok)[1], ]
    cat(sprintf("%s: minimum %d days (CD8=%.1f, IgG=%.2f at stop, no rebound)\n",
        toupper(arm), best$duration, best$CD8_at_stop, best$IgG_at_stop))
  } else {
    # Find best available (lowest rebound)
    best <- r[which.min(r$V_day60 / r$V_at_stop), ]
    cat(sprintf("%s: no duration achieves full adaptive takeover (best=%d days, ratio=%.1fx)\n",
        toupper(arm), best$duration, best$V_day60 / best$V_at_stop))
  }
}
