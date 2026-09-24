# Guards on the claims the Royal Society manuscript makes about mechanism.
#
# test-regression.R pins the calibrated endpoints. This file pins the claims
# that go beyond them: that the solver crosses the memory threshold, that the
# favipiravir start-day reversal is a duration effect, which mechanisms carry
# the outcome, and that the identifiability ridge is real. If one of these
# fails, a sentence in the manuscript has become false.

with_project_root <- function(code) {
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(qsp_root)
  force(code)
}

trapz <- function(t, y) sum(diff(t) * (head(y, -1) + tail(y, -1)) / 2)

test_that("every post-exposure prophylaxis scenario integrates to completion", {
  # The hard if/else memory switch was an attracting discontinuity that pinned
  # LSODA at CD8_E = 1, silently truncating 6 of these 24 scenarios -- including
  # both days on which prophylaxis is reported to prevent establishment.
  # simulate_patient_presymptomatic() returns NULL on an incomplete trajectory.
  with_project_root({
    p <- get_parameters()
    for (arm in c("ribavirin", "favipiravir", "combination")) {
      for (d in 0:2) {
        sim <- simulate_patient_presymptomatic(p, arm, t_exposure = d,
                                               t_end = 26, dt = 0.1)
        expect_false(is.null(sim),
                     info = sprintf("%s, exposure day %d, did not complete",
                                    arm, d))
      }
    }
  })
})

test_that("the memory decay switch is continuous across its threshold", {
  with_project_root({
    p <- get_parameters()
    expect_equal(p$n_mem_switch, 50)

    # Outside +/-9% of the threshold the smoothed switch must still behave like
    # the step function it replaces (0.1 * d below, 1.0 * d above): it agrees to
    # within 1% of d outside CD8_E in [0.914, 1.094].
    sw <- function(x) 0.1 + 0.9 / (1 + x^(-p$n_mem_switch))
    expect_lt(abs(sw(0.91) - 0.1), 0.01)
    expect_lt(abs(sw(1.10) - 1.0), 0.01)
    # and the band is no wider than claimed
    expect_gt(abs(sw(0.93) - 0.1), 0.01)
    expect_gt(abs(sw(1.07) - 1.0), 0.01)
    # ... and be continuous at it, which the step function was not.
    expect_lt(abs(sw(1 + 1e-9) - sw(1 - 1e-9)), 1e-6)
  })
})

test_that("prophylaxis on days 0-1 prevents establishment in every active arm", {
  with_project_root({
    p <- get_parameters()
    for (arm in c("ribavirin", "favipiravir", "combination")) {
      for (d in 0:1) {
        sim <- simulate_patient_presymptomatic(p, arm, t_exposure = d,
                                               t_end = 26, dt = 0.1)
        expect_false(is.null(sim))
        expect_lt(max(sim$V), 1e4)      # the establishment threshold
      }
    }
    # Day 2 is the documented breakthrough point; if this stops being true the
    # "days 0-1" claim needs rewording.
    sim2 <- simulate_patient_presymptomatic(p, "combination", t_exposure = 2,
                                            t_end = 26, dt = 0.1)
    expect_false(is.null(sim2))
    expect_gt(max(sim2$V), 1e4)
  })
})

test_that("the favipiravir day-1 penalty is a duration effect, not a timing one", {
  # Under a fixed 15-day course, day 1 is worse than day 2 (the published
  # "paradox"); with a course long enough to outlast the viral course, day 1 is
  # better, as it is for every other regimen.
  with_project_root({
    p <- get_parameters()
    run <- function(start, dur) {
      y0 <- build_symptom_onset_state(p)
      times <- seq(0, 42, by = 0.1)
      sched <- favipiravir_dosing_schedule(regimen = "standard",
                                           start_day = start,
                                           duration_days = dur)
      ev <- build_deSolve_events(sched, "C_FAV_gut",
                                 body_weight_kg = 75, Vd = p$Vd_FAV)
      ev <- ev[ev$time >= 0 & ev$time <= 42, ]
      times <- sort(unique(c(times, ev$time)))
      out <- as.data.frame(deSolve::lsoda(
        y = y0, times = times, func = hantavirus_qsp_ode, parms = p,
        events = list(data = ev), atol = 1e-10, rtol = 1e-8,
        maxsteps = 1e5, hmax = 0.5))
      expect_equal(nrow(out), length(times))
      t_dose_end <- start + dur
      list(m = extract_endpoints_from_simulation(out, p, "HFRS")$mortality_prob,
           t_K = out$time[which.max(out$K)],
           # Target cells surviving to the end of dosing, the quantity the
           # target-cell-sparing claim rests on. Same definition as
           # R/duration_start_interaction.R.
           T_frac = out$T[which.min(abs(out$time - t_dose_end))] / p$T_0,
           # The rebound is the viral maximum AFTER dosing stops. For a day-2
           # start the global maximum is the pre-treatment peak instead, so the
           # two must not be conflated.
           V_rebound = max(out$V[out$time > t_dose_end - 1]))
    }

    short_d1 <- run(1, 15)
    short_d2 <- run(2, 15)
    long_d1  <- run(1, 35)
    long_d2  <- run(2, 35)

    expect_gt(short_d1$m, short_d2$m)   # reversal present under a short course
    expect_lt(long_d1$m,  long_d2$m)    # and absent once the course is long
    expect_lt(long_d1$m,  0.04)         # 2.9% in the manuscript
    expect_gt(short_d1$m, 0.08)         # 9.2% in the manuscript

    # Target-cell sparing: 63.9% of the pool survives a day-1 course against
    # 28.7% for day-2 (abstract and 3.4). An earlier draft reported the day-2
    # figure as 23.0%, which the model has never given.
    expect_equal(short_d1$T_frac, 0.639, tolerance = 2e-3)
    expect_equal(short_d2$T_frac, 0.287, tolerance = 4e-3)
    expect_gt(short_d1$T_frac, short_d2$T_frac)

    # Rebound magnitudes: 3.33e7 for day-1 and 1.64e7 for day-2 (3.4).
    expect_equal(short_d1$V_rebound, 3.33e7, tolerance = 5e-3)
    expect_equal(short_d2$V_rebound, 1.64e7, tolerance = 5e-3)

    # The injury peak under a short day-1 course falls outside the 21-day
    # window used for the primary analysis: this is why that analysis
    # understated the effect roughly sixfold (2.29 vs 0.36 percentage points).
    expect_gt(short_d1$t_K, 21)
    expect_lt(long_d1$t_K, 21)
  })
})

test_that("mortality is renal and virus-driven, not cytokine-driven", {
  # The manuscript states that the cytokine term carries 0.09% of mortality,
  # NK cells 0.002% of infected-cell clearance and antibodies 0.14% of viral
  # clearance. These are the numbers that forbid reading the model as evidence
  # about cytokine-driven pathology or humoral protection.
  with_project_root({
    p  <- get_parameters()
    pl <- simulate_patient(p, "placebo", t_start = 0, t_end = 21)
    t  <- pl$time

    K <- max(pl$K); L <- max(pl$L); C <- max(pl$C_pro)
    Ka <- trapz(t, pl$K); La <- trapz(t, pl$L); Ca <- trapz(t, pl$C_pro)
    peak <- c(p$w_K_HFRS * K / (K + p$K50),
              p$w_L_HFRS * L / (L + p$L50),
              p$w_C_HFRS * C / (C + p$C50))
    auc  <- c(p$w_K_HFRS * Ka / (Ka + p$K_auc50),
              p$w_L_HFRS * La / (La + p$L_auc50),
              p$w_C_HFRS * Ca / (Ca + p$C_auc50))
    share <- ((1 - p$w_auc) * peak + p$w_auc * auc) /
      ((1 - p$w_auc) * sum(peak) + p$w_auc * sum(auc))

    expect_gt(share[1], 0.80)           # renal, 85.3% in the manuscript
    expect_lt(share[3], 0.005)          # cytokine, 0.090%

    nk  <- trapz(t, p$delta_NK * pl$NK * pl$I)
    nat <- trapz(t, p$delta_natural * pl$I)
    wm  <- pl$I / (200 + pl$I)
    cd8 <- trapz(t, wm * (p$delta_CD8 * pl$CD8_E * pl$I) +
                   (1 - wm) * (5.0 * pl$CD8_E * pl$I))
    expect_lt(nk / (nk + nat + cd8), 1e-4)          # NK, 0.0019%

    vc  <- trapz(t, p$c * pl$V)
    vab <- trapz(t, p$k_neut_IgM * pl$IgM * pl$V) +
      trapz(t, p$k_neut_IgG * pl$IgG * pl$V)
    expect_lt(vab / (vc + vab), 0.01)               # antibodies, 0.14%

    # Peak viraemia is target-cell-limited, not interferon-limited. Two forms,
    # both stated in the manuscript: against I_peak -- which occurs half a day
    # before the viral peak -- the agreement is 1.2%; against I at the time of
    # the viral peak, the true quasi-steady-state comparison, it is 0.03%. The
    # earlier 5% tolerance here was loose enough to pass while the manuscript
    # claimed 1%, which it was not.
    expect_lt(abs(max(pl$V) - p$p * max(pl$I) / p$c) / max(pl$V), 0.015)
    i_at_Vpeak <- pl$I[which.max(pl$V)]
    expect_lt(abs(max(pl$V) - p$p * i_at_Vpeak / p$c) / max(pl$V), 0.001)
    expect_lt(max(p$E_max_IFN * pl$F_I / (p$EC50_F + pl$F_I)), 1e-4)
  })
})

test_that("ribavirin's organ protection is endothelial, not immunomodulatory", {
  with_project_root({
    p <- get_parameters()
    m <- function(pp) {
      s <- simulate_patient(pp, "ribavirin", t_start = 1, t_end = 21)
      expect_false(is.null(s))
      extract_endpoints_from_simulation(s, pp, "HFRS")$mortality_prob
    }
    base <- m(p)

    p_imm <- p; p_imm$rbv_immuno_scale <- 0
    p_end <- p; p_end$rbv_endothelial_scale <- 0
    p_all <- p
    p_all$rbv_immuno_scale <- 0
    p_all$rbv_endothelial_scale <- 0
    p_all$rbv_recovery_scale <- 0

    # Removing immunomodulation does essentially nothing ...
    expect_lt(abs(m(p_imm) - base), 5e-4)
    # ... while removing endothelial protection accounts for the whole effect.
    expect_gt(m(p_end) - base, 0.03)
    expect_lt(abs(m(p_end) - m(p_all)), 1e-3)
  })
})

test_that("the identifiability ridge leaves calibration targets invariant", {
  # p -> p/f, beta -> beta*f, V_ref -> V_ref/f. Every calibration-visible
  # quantity must hold while peak viral load and the treatment effect move.
  with_project_root({
    p <- get_parameters()
    f <- 10
    q <- p
    q$p <- p$p / f; q$beta <- p$beta * f; q$V_ref <- p$V_ref / f

    pl0 <- simulate_patient(p, "placebo", t_start = 0, t_end = 21)
    pl1 <- simulate_patient(q, "placebo", t_start = 0, t_end = 21)
    e0 <- extract_endpoints_from_simulation(pl0, p, "HFRS")
    e1 <- extract_endpoints_from_simulation(pl1, q, "HFRS")

    # Invariant to within a fraction of a percent ...
    expect_equal(e1$mortality_prob, e0$mortality_prob, tolerance = 5e-3)
    expect_equal(e1$dialysis_prob,  e0$dialysis_prob,  tolerance = 5e-3)
    expect_equal(max(pl1$K), max(pl0$K), tolerance = 5e-3)
    expect_equal(min(pl1$PLT), min(pl0$PLT), tolerance = 5e-3)

    # ... while peak viraemia moves by about the full factor f.
    expect_lt(max(pl1$V) / max(pl0$V), 0.2)

    # ... and so does the predicted day-1 combination benefit (90.0% -> 60.8%).
    cb0 <- simulate_patient(p, "combination", t_start = 1, t_end = 21)
    cb1 <- simulate_patient(q, "combination", t_start = 1, t_end = 21)
    rrr <- function(pp, plc, trt) {
      a <- extract_endpoints_from_simulation(plc, pp, "HFRS")$mortality_prob
      b <- extract_endpoints_from_simulation(trt, pp, "HFRS")$mortality_prob
      100 * (a - b) / a
    }
    expect_equal(rrr(p, pl0, cb0), 90.0, tolerance = 1e-2)
    expect_lt(rrr(q, pl1, cb1), 70)
  })
})
