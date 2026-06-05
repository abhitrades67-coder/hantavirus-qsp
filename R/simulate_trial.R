#' Virtual Trial Simulation
#'
#' Runs the 4-arm virtual trial (Placebo, Ribavirin, Favipiravir, Combination)
#' across early/mid/late treatment windows for each virtual patient.
#'
#' @name simulate_trial
NULL

#' Build the infection-time state vector and initial conditions
#'
#' @param pars Parameter list (may be modified by vpop)
#' @return Named numeric vector of initial conditions
#' @export
build_initial_state <- function(pars) {
  # ===========================================================================
  # Adaptive immunity initial conditions
  #
  # The simulation starts at symptom onset (t=0), but infection began ~10 days
  # earlier. By symptom onset, the adaptive immune system has already been
  # priming for several days. We set small non-zero initial conditions for
  # CD4, CD8, IgM, IgG to reflect this pre-priming.
  #
  # Literature: Hantavirus incubation period = 7-42 days (median ~14 days).
  # T cell responses detectable by day 5-7 post-infection.
  # IgM detectable by day 7-10 post-infection.
  # (PMID: 19072554, 2902106, 36590594)
  #
  # These values are ~5-10% of peak levels seen in placebo simulations,
  # representing early effector cell activation during the incubation period.
  # ===========================================================================
  # Adaptive strength modulates pre-existing immunity (from virtual population)
  adapt_strength <- if ("adaptive_strength" %in% names(pars)) {
    pars$adaptive_strength
  } else {
    1.0
  }

  # Baseline adaptive levels at symptom onset.
  # CD8_N (naive) starts with a small resident population — these are the
  # antigen-inexperienced precursors. CD8_E (effector) starts at zero because
  # effectors only appear after CD8_N are primed and mature (k_mat = 0.14/day,
  # ~7-day delay). Placing the baseline in CD8_N instead of CD8_E prevents
  # premature killing of infected cells before de novo priming occurs.
  # CD4, IgM, IgG follow the same logic: small naive CD4, trace IgM from
  # innate signals, no IgG until class switch.
  # (PMID: 19072554, 2902106, 36590594)
  CD8_N_init <- 0.5 * adapt_strength  # AU, small naive CD8 precursor pool
  CD8_E_init <- 0                      # AU, no pre-existing effectors
  CD4_init   <- 0.5 * adapt_strength   # AU, small naive CD4 population
  IgM_init   <- 0.1 * adapt_strength   # AU, trace IgM from innate signals
  IgG_init   <- 0                      # AU, no IgG (requires class switch)

  c(
    T          = pars$T_0,
    I          = pars$I_0,
    V          = pars$V_0,
    NSs        = 0,
    F_I        = 0,
    NK         = 0,
    F_II       = 0,
    C_pro      = 0,
    C_anti     = 0,
    P          = 0,
    PLT        = pars$PLT_0,
    K          = 0,
    L          = 0,
    CD8_N      = CD8_N_init,
    CD8_E      = CD8_E_init,
    CD4        = CD4_init,
    IgM        = IgM_init,
    IgG        = IgG_init,
    C_RBV      = 0,
    C_FAV_gut  = 0,
    C_FAV      = 0,
    C_FAVI_RTP = 0,
    Hgb_drop   = 0
  )
}

#' Generate symptom-onset initial conditions via untreated burn-in
#'
#' Therapeutic simulations use t = 0 as symptom onset. This advances the
#' infection-time model without drug exposure for pars$symptom_onset_day days
#' and uses that state as the starting condition for post-symptom simulations.
#'
#' @param pars Parameter list
#' @param pre_dt Burn-in output step in days
#' @return Named numeric state vector at symptom onset
#' @export
build_symptom_onset_state <- function(pars, pre_dt = 0.1) {
  symptom_onset_day <- if ("symptom_onset_day" %in% names(pars)) {
    pars$symptom_onset_day
  } else {
    10
  }

  y0 <- build_initial_state(pars)
  if (symptom_onset_day <= 0) return(y0)

  times_pre <- seq(0, symptom_onset_day, by = pre_dt)
  if (tail(times_pre, 1) < symptom_onset_day) {
    times_pre <- c(times_pre, symptom_onset_day)
  }

  out_pre <- deSolve::lsoda(
    y     = y0,
    times = times_pre,
    func  = hantavirus_qsp_ode,
    parms = pars,
    atol  = 1e-10,
    rtol  = 1e-8,
    maxsteps = 20000
  )

  state <- as.numeric(out_pre[nrow(out_pre), -1])
  names(state) <- colnames(out_pre)[-1]

  for (nm in c("C_RBV", "C_FAV_gut", "C_FAV", "C_FAVI_RTP", "Hgb_drop")) {
    if (nm %in% names(state)) state[nm] <- 0
  }

  state
}

#' Run a single patient simulation
#'
#' @param pars Parameter list (patient-specific)
#' @param arm Treatment arm: "placebo", "ribavirin", "favipiravir", "combination"
#' @param t_start Treatment start day (relative to symptom onset)
#' @param t_end Simulation end day (relative to symptom onset)
#' @param dt Output time step (days)
#' @return Data frame with time-series output, or NULL on error
#' @export
simulate_patient <- function(pars, arm, t_start = 3, t_end = 21, dt = 0.1) {
  # Source the ODE model
  source("model/hantavirus_qsp.R")
  source("model/parameters.R")
  source("R/pk_models.R")

  # Simulate from symptom onset to t_end (treatment starts at t_start)
  times <- seq(0, t_end, by = dt)
  y0    <- build_symptom_onset_state(pars)

  body_weight_kg <- if ("body_weight_kg" %in% names(pars)) {
    pars$body_weight_kg
  } else {
    75
  }

  # Build dosing events based on arm
  all_events <- data.frame(
    time = numeric(), var = character(),
    value = numeric(), method = character(),
    stringsAsFactors = FALSE
  )

  if (arm %in% c("ribavirin", "combination")) {
    rbv_sched <- ribavirin_dosing_schedule(
      body_weight_kg = body_weight_kg,
      start_day      = t_start
    )
    rbv_events <- build_deSolve_events(
      rbv_sched, "C_RBV",
      body_weight_kg = body_weight_kg,
      Vd = pars$Vd_RBV
    )
    if (t_start == 0 && arm == "ribavirin" && nrow(all_events) == 0) {
      message(sprintf("[DEBUG] Ribavirin events: %d events, first dose at day %.1f, dose=%.0f mg",
                       nrow(rbv_events), rbv_events$time[1], rbv_events$value[1]))
    }
    all_events <- rbind(all_events, rbv_events)
  }

  if (arm %in% c("favipiravir", "combination")) {
    fav_sched <- favipiravir_dosing_schedule(
      regimen   = "standard",
      start_day = t_start
    )
    fav_events <- build_deSolve_events(
      fav_sched, "C_FAV_gut",
      body_weight_kg = body_weight_kg,
      Vd = pars$Vd_FAV
    )
    all_events <- rbind(all_events, fav_events)
  }

  # Ensure events are within simulation window
  if (nrow(all_events) > 0) {
    all_events <- all_events[all_events$time >= 0 &
                              all_events$time <= t_end, ]
    all_events <- all_events[order(all_events$time), ]
    rownames(all_events) <- NULL

    # deSolve applies event doses at exact output times. Dosing schedules use
    # hourly intervals, so merge event times into the requested output grid.
    if (nrow(all_events) > 0) {
      times <- sort(unique(c(times, all_events$time)))
    }
  }

  # Run ODE solver
  tryCatch({
    # --- Viral establishment pre-check ---
    # Some patients (~20%) sit at the bistable infection threshold where the
    # solver takes 25-230s hitting maxsteps. Run a fast coarse check first:
    # if I never exceeds 1e5 by day 5, the infection burns out → mortality ≈ 0.
    pre_times <- seq(0, min(5, t_end), by = 1)
    pre_out <- deSolve::lsoda(
      y     = y0,
      times = pre_times,
      func  = hantavirus_qsp_ode,
      parms = pars,
      atol   = 1e-6,
      rtol   = 1e-4,
      maxsteps = 5000
    )
    pre_out <- as.data.frame(pre_out)
    I_max_pre <- max(pre_out$I, na.rm = TRUE)
    V_max_pre <- max(pre_out$V, na.rm = TRUE)

    if (I_max_pre < 1e5 && V_max_pre < 1e6) {
      # Infection did not establish — return burn-out result immediately
      n_out <- length(times)
      burn_out <- data.frame(
        time = times,
        T    = rep(pars$T_0, n_out),
        I    = rep(0, n_out),
        V    = rep(0, n_out),
        F_I  = rep(0, n_out),
        NK   = rep(0, n_out),
        F_II = rep(0, n_out),
        C_pro = rep(0, n_out),
        P    = rep(0, n_out),
        K    = rep(0, n_out),
        L    = rep(0, n_out),
        PLT  = rep(2.5e5, n_out),
        CD4  = rep(0, n_out),
        CD8_N = rep(0, n_out),
        CD8_E = rep(0, n_out),
        IgM  = rep(0, n_out),
        IgG  = rep(0, n_out),
        NSs  = rep(0, n_out),
        T_reg = rep(0, n_out),
        C_RBV = rep(0, n_out),
        C_FAV_gut = rep(0, n_out),
        C_FAV_plasma = rep(0, n_out),
        C_FAVI_RTP = rep(0, n_out),
        Hgb  = rep(15, n_out),
        CD8_M = rep(0, n_out),
        F_III = rep(0, n_out),
        B_cell = rep(0, n_out),
        IL6  = rep(0, n_out),
        IL10 = rep(0, n_out),
        TNFa = rep(0, n_out),
        D_dimer = rep(0, n_out),
        check.names = FALSE
      )
      return(burn_out)
    }

    if (nrow(all_events) > 0) {
      out <- deSolve::lsoda(
        y     = y0,
        times = times,
        func  = hantavirus_qsp_ode,
        parms = pars,
        events = list(data = all_events),
        atol   = 1e-10,
        rtol   = 1e-8,
        maxsteps = 100000,
        hmax    = 0.5
      )
    } else {
      out <- deSolve::lsoda(
        y     = y0,
        times = times,
        func  = hantavirus_qsp_ode,
        parms = pars,
        atol   = 1e-10,
        rtol   = 1e-8,
        maxsteps = 100000,
        hmax    = 0.5
      )
    }

    out <- as.data.frame(out)

    # Validate: no NaN/Inf in key states
    if (any(is.nan(out$V)) || any(is.infinite(out$V))) {
      warning("Simulation produced NaN/Inf in viral load")
      return(NULL)
    }

    out
  }, error = function(e) {
    warning(sprintf("Simulation failed for arm=%s: %s", arm, e$message))
    NULL
  })
}

#' Run the full virtual trial
#'
#' Simulates all 4 arms across all patients. For efficiency, runs a
#' representative subset when N is large. Uses parallel processing
#' (foreach + doParallel) to speed up patient-level simulations.
#'
#' @param pop Virtual population data frame
#' @param pars Base parameter list
#' @param arms Character vector of treatment arms
#' @param n_patients Number of patients per arm (samples from pop)
#' @param t_start Simulation start (default 0 = symptom onset)
#' @param t_end Simulation end day (default 21)
#' @param dt Output time step (days)
#' @param seed Random seed
#' @param n_cores Number of CPU cores for parallel processing (default:
#'   parallel::detectCores() - 1, minimum 1)
#' @return List with elements:
#'   - results: list of data frames (one per arm)
#'   - endpoints: data frame of summary endpoints
#'   - pop_used: the patient subset used
#' @export
run_virtual_trial <- function(pop, pars,
                               arms = c("placebo", "ribavirin",
                                        "favipiravir", "combination"),
                               n_patients = 200,
                               t_start = 0, t_end = 21, dt = 0.1,
                               seed = 123,
                               n_cores = NULL) {
  set.seed(seed)

  # Source dependencies
  source("model/hantavirus_qsp.R")
  source("model/parameters.R")
  source("R/pk_models.R")
  source("R/pd_models.R")
  source("R/virtual_population.R")

  # Sample patients
  if (nrow(pop) > n_patients) {
    idx <- sample(nrow(pop), n_patients)
    pop_used <- pop[idx, ]
  } else {
    pop_used <- pop
  }

  # Set up parallel backend
  available_cores <- parallel::detectCores()
  if (is.null(n_cores)) {
    n_cores <- max(1, available_cores - 1)
  } else {
    n_cores <- min(n_cores, available_cores)
  }
  message(sprintf("Using %d cores for parallel processing", n_cores))

  # Build list of all (arm, patient_index) tasks
  task_grid <- expand.grid(
    arm = arms,
    idx = seq_len(nrow(pop_used)),
    stringsAsFactors = FALSE
  )

  # Run all simulations in parallel using parLapply
  message(sprintf("Running %d patient-simulations across %d arms...",
                   nrow(task_grid), length(arms)))

  cl <- parallel::makeCluster(n_cores, type = "PSOCK")
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # Export everything needed to workers
  parallel::clusterExport(cl,
    varlist = c("task_grid", "pop_used", "pars", "t_end", "dt"),
    envir = environment()
  )

  endpoint_list <- parallel::parLapply(cl, seq_len(nrow(task_grid)),
    function(i) {
      source("model/hantavirus_qsp.R", local = TRUE)
      source("model/parameters.R", local = TRUE)
      source("R/pk_models.R", local = TRUE)
      source("R/pd_models.R", local = TRUE)
      source("R/virtual_population.R", local = TRUE)
      source("R/simulate_trial.R", local = TRUE)

      arm <- task_grid$arm[i]
      pidx <- task_grid$idx[i]
      row <- pop_used[pidx, ]
      pars_i <- apply_vpop_to_model(row, pars)

      t_treat <- if (arm == "placebo") {
        0
      } else if ("time_to_treatment_days" %in% names(row)) {
        row$time_to_treatment_days
      } else {
        t_start
      }

      sim_out <- simulate_patient(
        pars    = pars_i,
        arm     = arm,
        t_start = t_treat,
        t_end   = t_end,
        dt      = dt
      )

      if (is.null(sim_out)) return(NULL)

      pk  <- extract_peak_biomarkers(sim_out)
      auc <- compute_viral_AUC(sim_out)
      ttc <- time_to_clearance(sim_out)
      ep  <- extract_endpoints_from_simulation(sim_out, pars_i,
                                                syndrome = row$syndrome)

      data.frame(
        patient_id     = row$patient_id,
        arm            = arm,
        syndrome       = row$syndrome,
        treatment_start_day = t_treat,
        V_peak         = pk$V_peak,
        V_AUC          = auc,
        time_clearance = ttc,
        PLT_nadir      = pk$PLT_nadir,
        C_pro_peak     = pk$C_pro_peak,
        K_peak         = pk$K_peak,
        L_peak         = pk$L_peak,
        CD8_E_peak     = pk$CD8_E_peak,
        CD4_peak       = pk$CD4_peak,
        IgM_peak       = pk$IgM_peak,
        IgG_peak       = pk$IgG_peak,
        dialysis_prob  = ep$dialysis_prob,
        ecmo_prob      = ep$ecmo_prob,
        mortality_prob = ep$mortality_prob,
        Hgb_drop_max   = pk$Hgb_drop_max,
        stringsAsFactors = FALSE
      )
    }
  )

  # Combine results
  valid_results <- Filter(function(x) is.data.frame(x), endpoint_list)
  if (length(valid_results) == 0) {
    stop("All simulations failed!")
  }
  endpoints <- do.call(rbind, valid_results)

  # Build per-arm results
  results <- list()
  for (arm in arms) {
    results[[arm]] <- endpoints[endpoints$arm == arm, ]
    message(sprintf("  %s: %d/%d patients completed",
                     arm, nrow(results[[arm]]), nrow(pop_used)))
  }

  list(
    results   = results,
    endpoints = endpoints,
    pop_used  = pop_used
  )
}

#' Run a single patient simulation with pre-symptomatic (post-exposure) treatment
#'
#' Starts at t=0 (infection, V_0 seeded) and runs the full ODE with drug events
#' starting at t_exposure days post-infection. No drug-free burn-in period.
#' This simulates contact-tracing-initiated prophylaxis during the incubation
#' period (days 0-5 post-exposure, before symptom onset).
#'
#' @param pars Parameter list (patient-specific)
#' @param arm Treatment arm: "placebo", "ribavirin", "favipiravir", "combination"
#' @param t_exposure Treatment start day post-infection (0-5; 0 = immediate PEP)
#' @param t_end Simulation end day (default 26 = 5-day incubation + 21-day follow-up)
#' @param dt Output time step (days)
#' @return Data frame with time-series output, or NULL on error
#' @export
simulate_patient_presymptomatic <- function(pars, arm, t_exposure = 0,
                                             t_end = 26, dt = 0.1) {
  source("model/hantavirus_qsp.R")
  source("model/parameters.R")
  source("R/pk_models.R")

  # Start directly from infection state — no drug-free burn-in
  times <- seq(0, t_end, by = dt)
  y0    <- build_initial_state(pars)

  body_weight_kg <- if ("body_weight_kg" %in% names(pars)) {
    pars$body_weight_kg
  } else {
    75
  }

  # Build dosing events based on arm
  all_events <- data.frame(
    time = numeric(), var = character(),
    value = numeric(), method = character(),
    stringsAsFactors = FALSE
  )

  if (arm %in% c("ribavirin", "combination")) {
    rbv_sched <- ribavirin_dosing_schedule(
      body_weight_kg = body_weight_kg,
      start_day      = t_exposure
    )
    rbv_events <- build_deSolve_events(
      rbv_sched, "C_RBV",
      body_weight_kg = body_weight_kg,
      Vd = pars$Vd_RBV
    )
    all_events <- rbind(all_events, rbv_events)
  }

  if (arm %in% c("favipiravir", "combination")) {
    fav_sched <- favipiravir_dosing_schedule(
      regimen   = "standard",
      start_day = t_exposure
    )
    fav_events <- build_deSolve_events(
      fav_sched, "C_FAV_gut",
      body_weight_kg = body_weight_kg,
      Vd = pars$Vd_FAV
    )
    all_events <- rbind(all_events, fav_events)
  }

  # Ensure events are within simulation window
  if (nrow(all_events) > 0) {
    all_events <- all_events[all_events$time >= 0 &
                              all_events$time <= t_end, ]
    all_events <- all_events[order(all_events$time), ]
    rownames(all_events) <- NULL

    if (nrow(all_events) > 0) {
      times <- sort(unique(c(times, all_events$time)))
    }
  }

  # Run ODE solver
  tryCatch({
    if (nrow(all_events) > 0) {
      out <- deSolve::lsoda(
        y     = y0,
        times = times,
        func  = hantavirus_qsp_ode,
        parms = pars,
        events = list(data = all_events),
        atol   = 1e-10,
        rtol   = 1e-8,
        maxsteps = 100000,
        hmax    = 0.5
      )
    } else {
      out <- deSolve::lsoda(
        y     = y0,
        times = times,
        func  = hantavirus_qsp_ode,
        parms = pars,
        atol   = 1e-10,
        rtol   = 1e-8,
        maxsteps = 100000,
        hmax    = 0.5
      )
    }

    out <- as.data.frame(out)

    if (any(is.nan(out$V)) || any(is.infinite(out$V))) {
      warning("Pre-symptomatic simulation produced NaN/Inf in viral load")
      return(NULL)
    }

    out
  }, error = function(e) {
    warning(sprintf("Pre-symptomatic simulation failed for arm=%s, t_exposure=%d: %s",
                     arm, t_exposure, e$message))
    NULL
  })
}


#' Run treatment-window analysis
#'
#' For a representative patient, simulate treatment at days 1-7 and
#' compute benefit metrics.
#'
#' @param pars Parameter list
#' @param arms Arms to compare
#' @param window_days Vector of treatment start days to test
#' @param t_end Simulation end day
#' @return Data frame with window analysis results
#' @export
run_window_analysis <- function(pars,
                                 arms = c("ribavirin", "favipiravir", "combination"),
                                 window_days = 1:7,
                                 t_end = 21) {
  source("model/hantavirus_qsp.R")
  source("R/pk_models.R")
  source("R/pd_models.R")

  results <- data.frame(
    arm            = character(),
    treatment_day  = integer(),
    V_peak         = numeric(),
    V_AUC          = numeric(),
    time_clearance = numeric(),
    K_peak         = numeric(),
    L_peak         = numeric(),
    mortality_prob = numeric(),
    stringsAsFactors = FALSE
  )

  for (arm in arms) {
    for (day in window_days) {
      sim_out <- simulate_patient(pars, arm, t_start = day, t_end = t_end)

      if (is.null(sim_out)) next

      pk <- extract_peak_biomarkers(sim_out)
      auc <- compute_viral_AUC(sim_out)
      ttc <- time_to_clearance(sim_out)
      ep  <- extract_endpoints_from_simulation(sim_out, pars, syndrome = "HFRS")

      results <- rbind(results, data.frame(
        arm            = arm,
        treatment_day  = day,
        V_peak         = pk$V_peak,
        V_AUC          = auc,
        time_clearance = ttc,
        K_peak         = pk$K_peak,
        L_peak         = pk$L_peak,
        mortality_prob = ep$mortality_prob,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Add placebo baseline
  sim_placebo <- simulate_patient(pars, "placebo", t_start = 0, t_end = t_end)
  if (!is.null(sim_placebo)) {
    pk_p <- extract_peak_biomarkers(sim_placebo)
    ep_p <- extract_endpoints_from_simulation(sim_placebo, pars, syndrome = "HFRS")
    for (day in window_days) {
      results <- rbind(results, data.frame(
        arm            = "placebo",
        treatment_day  = day,
        V_peak         = pk_p$V_peak,
        V_AUC          = compute_viral_AUC(sim_placebo),
        time_clearance = time_to_clearance(sim_placebo),
        K_peak         = pk_p$K_peak,
        L_peak         = pk_p$L_peak,
        mortality_prob = ep_p$mortality_prob,
        stringsAsFactors = FALSE
      ))
    }
  }

  results
}
