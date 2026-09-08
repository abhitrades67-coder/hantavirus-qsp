#' Post-Exposure Prophylaxis (PEP) Analysis
#'
#' Evaluates whether starting antiviral treatment during the incubation period
#' (days 0-5 post-exposure, before symptom onset) can prevent disease onset.
#'
#' Simulates contact-tracing scenarios where exposed individuals are identified
#' and treated before symptoms appear.
#'
#' @export

suppressPackageStartupMessages({
  library(deSolve)
})

# Source model components
source("model/hantavirus_qsp.R")
source("model/parameters.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/simulate_trial.R")

# Get default parameters
pars <- get_parameters()

# Post-exposure treatment days (post-infection)
# Day 0 = immediate PEP at exposure
# Day 5 = treatment at symptom onset (equivalent to current "Day 1" post-symptom)
pep_days <- 0:5

# Treatment arms
arms <- c("placebo", "ribavirin", "favipiravir", "combination")

# Disease establishment threshold: V_peak > 1e4 copies/mL
# If V stays below this, the infection was contained by treatment
V_ESTABLISH_THRESHOLD <- 1e4

cat("============================================================\n")
cat("  POST-EXPOSURE PROPHYLAXIS (PEP) ANALYSIS\n")
cat("============================================================\n\n")
cat(sprintf("Treatment days post-exposure: %s\n", paste(pep_days, collapse = ", ")))
cat(sprintf("Treatment arms: %s\n", paste(arms, collapse = ", ")))
cat(sprintf("Disease establishment threshold: V_peak > %.0e copies/mL\n",
            V_ESTABLISH_THRESHOLD))
cat(sprintf("Simulation length: 26 days (5-day incubation + 21-day follow-up)\n\n"))

# Run analysis
results <- data.frame(
  arm               = character(),
  pep_day           = integer(),
  V_peak            = numeric(),
  V_AUC             = numeric(),
  K_peak            = numeric(),
  L_peak            = numeric(),
  PLT_nadir         = numeric(),
  C_pro_peak        = numeric(),
  mortality_prob    = numeric(),
  dialysis_prob     = numeric(),
  ecmo_prob         = numeric(),
  established       = logical(),
  time_clearance    = numeric(),
  stringsAsFactors  = FALSE
)

for (arm in arms) {
  for (pep_day in pep_days) {
    cat(sprintf("  Running %s at day %d post-exposure... ", arm, pep_day))

    sim_out <- simulate_patient_presymptomatic(
      pars       = pars,
      arm        = arm,
      t_exposure = pep_day,
      t_end      = 26,
      dt         = 0.1
    )

    if (is.null(sim_out)) {
      # Record the failure as a visible NA row rather than dropping it. A
      # dropped row silently shrinks the results table (six of twenty-four rows
      # vanished the first time the solver-truncation check was added), which is
      # indistinguishable from a scenario that was never run.
      cat("FAILED (solver returned an incomplete trajectory)\n")
      results <- rbind(results, data.frame(
        arm            = arm,
        pep_day        = pep_day,
        V_peak         = NA_real_,
        V_AUC          = NA_real_,
        K_peak         = NA_real_,
        L_peak         = NA_real_,
        PLT_nadir      = NA_real_,
        C_pro_peak     = NA_real_,
        mortality_prob = NA_real_,
        dialysis_prob  = NA_real_,
        ecmo_prob      = NA_real_,
        established    = NA,
        time_clearance = NA_real_,
        stringsAsFactors = FALSE
      ))
      next
    }

    # Extract endpoints
    V_peak     <- max(sim_out$V)
    # Trapezoidal integration on the ACTUAL time vector. The output grid is not
    # uniform: dose times are merged into it, so a fixed-width rectangle sum
    # (sum(V * 0.1)) over-counted drug-containing arms by up to 32% while
    # leaving placebo untouched -- an arm-dependent bias. Use the shared helper.
    V_AUC      <- compute_viral_AUC(sim_out)
    K_peak     <- max(sim_out$K)
    L_peak     <- max(sim_out$L)
    PLT_nadir  <- min(sim_out$PLT)
    C_pro_peak <- max(sim_out$C_pro)

    # Clinical endpoints come from the SINGLE canonical implementation in
    # R/pd_models.R. This script previously re-implemented them as a pair of
    # logistic functions of peak injury (breakpoints K=100, L=300), which is a
    # different model: it understated dialysis risk ~3-fold and, for arms where
    # infection was prevented (K ~ 0), returned a ~0.7% floor instead of ~0,
    # flattening the whole post-exposure gradient.
    endpoints      <- extract_endpoints_from_simulation(sim_out, pars,
                                                        syndrome = "HFRS")
    mortality_prob <- endpoints$mortality_prob
    dialysis_prob  <- endpoints$dialysis_prob
    ecmo_prob      <- endpoints$ecmo_prob

    # Disease establishment
    established <- V_peak > V_ESTABLISH_THRESHOLD

    # Time to clearance (V < 100 copies/mL after peak)
    clearance_idx <- which(sim_out$V < 100 &
                             sim_out$time > pars$symptom_onset_day)
    time_clearance <- if (length(clearance_idx) > 0 && established) {
      sim_out$time[min(clearance_idx)]
    } else if (!established) {
      0  # Never established
    } else {
      NA  # Did not clear
    }

    results <- rbind(results, data.frame(
      arm            = arm,
      pep_day        = pep_day,
      V_peak         = V_peak,
      V_AUC          = V_AUC,
      K_peak         = K_peak,
      L_peak         = L_peak,
      PLT_nadir      = PLT_nadir,
      C_pro_peak     = C_pro_peak,
      mortality_prob = mortality_prob,
      dialysis_prob  = dialysis_prob,
      ecmo_prob      = ecmo_prob,
      established    = established,
      time_clearance = time_clearance,
      stringsAsFactors = FALSE
    ))

    status <- if (established) "ESTABLISHED" else "PREVENTED"
    cat(sprintf("V_peak=%.1e, mortality=%.1f%% [%s]\n",
                V_peak, mortality_prob * 100, status))
  }
}

# Summary table
cat("\n============================================================\n")
cat("  POST-EXPOSURE PROPHYLAXIS RESULTS SUMMARY\n")
cat("============================================================\n\n")

for (arm in arms) {
  arm_data <- results[results$arm == arm, ]
  cat(sprintf("--- %s ---\n", toupper(arm)))
  cat(sprintf("  %-8s  %-12s  %-10s  %-10s  %-10s  %-10s  %s\n",
              "PEP Day", "V_peak", "K_peak", "L_peak", "Mortality", "Dialysis", "Established?"))
  cat(sprintf("  %-8s  %-12s  %-10s  %-10s  %-10s  %-10s  %s\n",
              "-------", "----------", "--------", "--------", "---------", "--------", "-----------"))
  for (i in seq_len(nrow(arm_data))) {
    r <- arm_data[i, ]
    cat(sprintf("  %-8d  %-12s  %-10.1f  %-10.1f  %-10s  %-10s  %s\n",
                r$pep_day,
                formatC(r$V_peak, format = "e", digits = 2),
                r$K_peak,
                r$L_peak,
                sprintf("%.1f%%", r$mortality_prob * 100),
                sprintf("%.1f%%", r$dialysis_prob * 100),
                ifelse(r$established, "YES", "NO")))
  }
  cat("\n")
}

# Disease prevention summary
cat("============================================================\n")
cat("  DISEASE PREVENTION SUMMARY\n")
cat("============================================================\n\n")
cat("Can single drugs prevent disease onset (V_peak < 1e4)?\n\n")

prevention_table <- reshape(
  results[, c("arm", "pep_day", "established")],
  idvar     = "pep_day",
  timevar   = "arm",
  direction = "wide"
)
names(prevention_table) <- gsub("established\\.", "", names(prevention_table))
prevention_table <- prevention_table[order(prevention_table$pep_day), ]

cat(sprintf("  %-8s  %-10s  %-10s  %-10s  %s\n",
            "PEP Day", "Placebo", "Ribavirin", "Favipiravir", "Combination"))
cat(sprintf("  %-8s  %-10s  %-10s  %-10s  %s\n",
            "-------", "---------", "---------", "-----------", "-----------"))
for (i in seq_len(nrow(prevention_table))) {
  r <- prevention_table[i, ]
  cat(sprintf("  %-8d  %-10s  %-10s  %-10s  %s\n",
              r$pep_day,
              ifelse(is.na(r$placebo), "-", ifelse(r$placebo, "FAIL", "PREVENT")),
              ifelse(is.na(r$ribavirin), "-", ifelse(r$ribavirin, "FAIL", "PREVENT")),
              ifelse(is.na(r$favipiravir), "-", ifelse(r$favipiravir, "FAIL", "PREVENT")),
              ifelse(is.na(r$combination), "-", ifelse(r$combination, "FAIL", "PREVENT"))))
}

# Save results
dir.create("outputs", showWarnings = FALSE)
write.csv(results, "outputs/preexposure_results.csv", row.names = FALSE)

# Write summary text
sink("outputs/preexposure_summary.txt")
cat("============================================================\n")
cat("  HANTAVIRUS QSP — POST-EXPOSURE PROPHYLAXIS ANALYSIS\n")
cat("============================================================\n\n")
cat(sprintf("Date: %s\n", Sys.Date()))
cat(sprintf("Treatment days: %s (post-exposure)\n", paste(pep_days, collapse = ", ")))
cat(sprintf("Arms: %s\n", paste(arms, collapse = ", ")))
cat(sprintf("Disease establishment threshold: V_peak > %.0e copies/mL\n",
            V_ESTABLISH_THRESHOLD))
cat("\n--- DETAILED RESULTS ---\n\n")

for (arm in arms) {
  arm_data <- results[results$arm == arm, ]
  cat(sprintf("ARM: %s\n", toupper(arm)))
  cat(sprintf("  %-8s  %-12s  %-10s  %-10s  %-10s  %-10s  %s\n",
              "PEP Day", "V_peak", "K_peak", "L_peak", "Mortality", "Dialysis", "Established?"))
  for (i in seq_len(nrow(arm_data))) {
    r <- arm_data[i, ]
    cat(sprintf("  %-8d  %-12s  %-10.1f  %-10.1f  %-10s  %-10s  %s\n",
                r$pep_day,
                formatC(r$V_peak, format = "e", digits = 2),
                r$K_peak, r$L_peak,
                sprintf("%.1f%%", r$mortality_prob * 100),
                sprintf("%.1f%%", r$dialysis_prob * 100),
                ifelse(r$established, "YES", "NO")))
  }
  cat("\n")
}

cat("\n--- DISEASE PREVENTION ---\n\n")
cat("PREVENT = treatment contained the inoculum (V_peak < threshold)\n")
cat("FAIL    = infection broke through despite treatment\n\n")
for (i in seq_len(nrow(prevention_table))) {
  r <- prevention_table[i, ]
  cat(sprintf("Day %d: Placebo=%s, Ribavirin=%s, Favipiravir=%s, Combination=%s\n",
              r$pep_day,
              ifelse(is.na(r$placebo), "-", ifelse(r$placebo, "FAIL", "PREVENT")),
              ifelse(is.na(r$ribavirin), "-", ifelse(r$ribavirin, "FAIL", "PREVENT")),
              ifelse(is.na(r$favipiravir), "-", ifelse(r$favipiravir, "FAIL", "PREVENT")),
              ifelse(is.na(r$combination), "-", ifelse(r$combination, "FAIL", "PREVENT"))))
}
sink()

cat("\n\nResults saved to:\n")
cat("  outputs/preexposure_results.csv\n")
cat("  outputs/preexposure_summary.txt\n")
cat("\nDone.\n")
