#' Hantavirus QSP Virtual Trial — Main Pipeline
#'
#' Runs the complete simulation pipeline:
#'   1. Load parameters
#'   2. Generate virtual population
#'   3. Run 4-arm virtual trial
#'   4. Run treatment-window analysis
#'   5. Generate all output figures and tables
#'
#' Usage:
#'   Rscript R/run_pipeline.R
#'
#' Or from R:
#'   source("R/run_pipeline.R")
#'
#' @export

#' String concatenation helper
`%+%` <- function(a, b) paste0(a, b)

main <- function() {
  start_time <- Sys.time()
  message("=" %+% paste(rep("=", 59), collapse = ""))
  message("  HANTAVIRUS QSP VIRTUAL TRIAL PIPELINE")
  message("=" %+% paste(rep("=", 59), collapse = ""))

  # ---- 0. Ensure dependencies ----
  required_packages <- c("deSolve", "ggplot2", "dplyr", "tidyr", "gridExtra",
                          "foreach", "doParallel")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf(
        "Package '%s' is required. Install with: install.packages('%s')",
        pkg, pkg
      ))
    }
  }

  # ---- 1. Load model and parameters ----
  message("\n[1/7] Loading model and parameters...")
  source("model/hantavirus_qsp.R")
  source("model/parameters.R")
  source("R/pk_models.R")
  source("R/pd_models.R")
  source("R/virtual_population.R")
  source("R/simulate_trial.R")
  source("R/analysis.R")

  pars <- get_parameters()
  message(sprintf("  Loaded %d model parameters", length(pars)))

  # ---- 2. Generate virtual population ----
  message("\n[2/7] Generating virtual population (N=1000)...")
  pop <- generate_virtual_population(N = 1000, seed = 42)
  pop_summary <- summarise_virtual_population(pop)
  print(pop_summary)

  # Save population
  dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(pop, "outputs/virtual_population.csv", row.names = FALSE)
  message("  Saved: outputs/virtual_population.csv")

  # ---- 3. Save parameter table ----
  message("\n[3/7] Saving parameter table...")
  param_tbl <- build_parameter_table(pars)
  utils::write.csv(param_tbl, "outputs/parameter_table.csv", row.names = FALSE)
  message("  Saved: outputs/parameter_table.csv")
  message(sprintf("  %d parameters documented", nrow(param_tbl)))

  # ---- 4. Run virtual trial (200 patients per arm for tractability) ----
  message("\n[4/7] Running virtual trial (4 arms x 200 patients)...")
  trial_output <- run_virtual_trial(
    pop        = pop,
    pars       = pars,
    arms       = c("placebo", "ribavirin", "favipiravir", "combination"),
    n_patients = 200,
    t_end      = 21,
    dt         = 0.5,
    seed       = 123
  )
  message(sprintf("  Completed: %d patient-simulations",
                   nrow(trial_output$endpoints)))

  # Diagnostic: check ribavirin V_peak distribution
  message("\n[DIAG] Ribavirin V_peak distribution:")
  rbv_idx <- trial_output$endpoints$arm == "ribavirin"
  rbv_vpeak <- trial_output$endpoints$V_peak[rbv_idx]
  message(sprintf("  Median: %.0f", median(rbv_vpeak)))
  message(sprintf("  Mean:   %.0f", mean(rbv_vpeak)))
  message(sprintf("  Min:    %.0f", min(rbv_vpeak)))
  message(sprintf("  Max:    %.0f", max(rbv_vpeak)))
  message(sprintf("  Q1:     %.0f", quantile(rbv_vpeak, 0.25)))
  message(sprintf("  Q3:     %.0f", quantile(rbv_vpeak, 0.75)))
  message(sprintf("  Fraction with V_peak < 1000: %.3f",
                   mean(rbv_vpeak < 1000)))

  # ---- 5. Run treatment-window analysis ----
  message("\n[5/7] Running treatment-window analysis...")
  window_data <- run_window_analysis(
    pars        = pars,
    arms        = c("ribavirin", "favipiravir", "combination"),
    window_days = 1:7,
    t_end       = 21
  )
  message(sprintf("  Window analysis: %d scenarios", nrow(window_data)))

  # Full-precision record of the treatment-window analysis. The text summary
  # rounds mortality to one decimal place, which is not enough for the
  # published table, so the same numbers are written here as a CSV and the
  # manuscript table is generated from this file.
  pl_sim <- simulate_patient(pars, "placebo", t_start = 0, t_end = 21)
  if (is.null(pl_sim)) stop("placebo reference simulation failed")
  pl_mort <- extract_endpoints_from_simulation(pl_sim, pars, "HFRS")$mortality_prob
  window_out <- window_data
  window_out$placebo_mortality <- pl_mort
  # Section 3.3 states each arm's viral exposure as a percentage of placebo.
  # Record the denominator here: it is this patient's own 21-day placebo
  # exposure, not the placebo row of the prophylaxis analysis, which runs on
  # a different horizon and gives a different figure.
  window_out$placebo_V_AUC <- compute_viral_AUC(pl_sim)
  window_out$RRR <- 100 * (pl_mort - window_out$mortality_prob) / pl_mort
  utils::write.csv(window_out, "outputs/treatment_window.csv", row.names = FALSE)
  message("  Saved: outputs/treatment_window.csv")

  # ---- 6. Generate figures ----
  message("\n[6/7] Generating figures...")

  # 6a. Viral kinetics — use individual simulation data
  # Re-run a small set for trajectory plots (5 patients per arm)
  # Simulate at fixed treatment start days (1–7) so panels are comparable
  message("  Generating trajectory data for plots...")
  traj_data <- data.frame(
    time            = numeric(),
    V               = numeric(),
    C_pro           = numeric(),
    P               = numeric(),
    PLT             = numeric(),
    K               = numeric(),
    L               = numeric(),
    CD8_N           = numeric(),
    CD8_E           = numeric(),
    CD4             = numeric(),
    IgM             = numeric(),
    IgG             = numeric(),
    C_RBV           = numeric(),
    C_FAV           = numeric(),
    C_FAVI_RTP      = numeric(),
    arm             = character(),
    patient_id      = integer(),
    treat_start_day = integer(),
    stringsAsFactors = FALSE
  )

  set.seed(99)
  plot_patients <- sample(pop$patient_id, 5)

  # Placebo — no treatment start; labeled as day 0 for faceting
  for (pid in plot_patients) {
    row    <- pop[pop$patient_id == pid, ]
    pars_i <- apply_vpop_to_model(row, pars)
    sim_out <- simulate_patient(pars_i, "placebo", t_start = 0, t_end = 21, dt = 0.5)
    if (is.null(sim_out)) next
    sim_out$arm             <- "placebo"
    sim_out$patient_id      <- pid
    sim_out$treat_start_day <- 0L
    traj_data <- rbind(traj_data, sim_out)
  }

  # Active arms at fixed treatment start days 1–7
  active_arms <- c("ribavirin", "favipiravir", "combination")
  for (start_day in 1:7) {
    for (arm in active_arms) {
      for (pid in plot_patients) {
        row    <- pop[pop$patient_id == pid, ]
        pars_i <- apply_vpop_to_model(row, pars)
        sim_out <- simulate_patient(pars_i, arm, t_start = start_day, t_end = 21, dt = 0.5)
        if (is.null(sim_out)) next
        sim_out$arm             <- arm
        sim_out$patient_id      <- pid
        sim_out$treat_start_day <- as.integer(start_day)
        traj_data <- rbind(traj_data, sim_out)
      }
    }
  }

  # Viral kinetics plot
  plot_viral_kinetics(
    traj_data[, c("time", "V", "arm", "patient_id", "treat_start_day")],
    outfile = "outputs/viral_kinetics.png"
  )

  # Biomarker trajectories (wide format)
  plot_biomarker_trajectories(
    traj_data[, c("time", "C_pro", "P", "PLT", "arm", "patient_id")],
    outfile = "outputs/biomarker_trajectories.png"
  )

  # Organ injury (wide format)
  plot_organ_injury(
    traj_data[, c("time", "K", "L", "arm", "patient_id")],
    outfile = "outputs/organ_injury.png"
  )

  # PK profiles
  pk_data <- traj_data[, c("time", "C_RBV", "C_FAV", "C_FAVI_RTP", "arm")]
  plot_pk_profiles(pk_data, outfile = "outputs/pk_profiles.png")

  # Treatment window
  plot_treatment_window(window_data,
                         outfile = "outputs/treatment_window.png")

  # Treatment window sensitivity analysis (multi-panel)
  plot_window_sensitivity(window_data,
                          trial_output = trial_output,
                          outfile = "outputs/treatment_window_sensitivity.png")

  # Adaptive immunity trajectories
  adaptive_data <- traj_data[, c("time", "CD8_E", "CD4", "IgM", "IgG",
                                  "arm", "patient_id")]
  plot_adaptive_immunity(adaptive_data,
                          outfile = "outputs/adaptive_immunity.png")

  # ---- 7. Generate summary outputs ----
  message("\n[7/7] Generating summary outputs...")

  # Clinical endpoints CSV
  ep_summary <- write_clinical_endpoints(
    trial_output,
    outfile = "outputs/clinical_endpoints.csv"
  )
  print(ep_summary)

  # Simulation summary text
  write_simulation_summary(trial_output, window_data,
                            outfile = "outputs/simulation_summary.txt")

  # ---- Done ----
  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  message("\n" %+% "=" %+% paste(rep("=", 59), collapse = ""))
  message(sprintf("  Pipeline complete in %.1f minutes", elapsed))
  message("=" %+% paste(rep("=", 59), collapse = ""))

  message("\nOutput files:")
  out_files <- list.files("outputs", full.names = TRUE)
  for (f in out_files) {
    message(sprintf("  %s (%.1f KB)", f, file.info(f)$size / 1024))
  }

  invisible(list(
    trial   = trial_output,
    window  = window_data,
    pop     = pop,
    params  = param_tbl
  ))
}

# Run if sourced directly
if (!interactive() && sys.nframe() == 0) {
  main()
}
