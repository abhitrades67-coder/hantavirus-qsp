# Regression guards on the numbers the manuscript reports.
#
# These pin the published results so that a change which silently moves a
# headline number fails a test instead of shipping. Every value below was
# measured from the model itself; if one of them changes, either the change was
# intended (update the constant AND the manuscript) or it is a bug.
#
# simulate_patient() source()s the model with paths relative to the project
# root, so these tests run from there and restore the working directory after.

with_project_root <- function(code) {
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(qsp_root)
  force(code)
}

test_that("the burn-out placeholder has exactly the model's state columns", {
  # Regression guard for a real defect: the non-establishment shortcut used to
  # return a hand-written frame that omitted C_anti, C_FAV and Hgb_drop and
  # invented ten states the model does not have, so max(sim$Hgb_drop) returned
  # -Inf and rbind() against real solver output failed.
  with_project_root({
    p <- get_parameters()
    p$beta <- 1e-9                      # guarantees non-establishment
    burn <- simulate_patient(p, "placebo", t_start = 0, t_end = 21)
    real <- names(build_initial_state(p))

    expect_setequal(names(burn), c("time", real))
    expect_true(all(vapply(burn, function(x) all(is.finite(x)), logical(1))))
    expect_true(is.finite(max(burn$Hgb_drop)))
    expect_equal(unname(burn$T[1]), p$T_0)
    expect_equal(unname(burn$PLT[1]), p$PLT_0)
  })
})

test_that("placebo HFRS endpoints match the published calibration", {
  with_project_root({
    p  <- get_parameters()
    pl <- simulate_patient(p, "placebo", t_start = 0, t_end = 21)
    expect_false(is.null(pl))

    ep <- extract_endpoints_from_simulation(pl, p, syndrome = "HFRS")
    expect_equal(ep$mortality_prob, 0.09584, tolerance = 1e-3)
    expect_equal(ep$dialysis_prob,  0.43110, tolerance = 1e-3)
    expect_equal(max(pl$V),         3.4183e7, tolerance = 1e-3)
    expect_equal(min(pl$PLT),       40236,    tolerance = 1e-3)
    expect_equal(max(pl$K),         53.044,   tolerance = 1e-3)
  })
})

test_that("day-1 combination gives the published 90% relative risk reduction", {
  with_project_root({
    p  <- get_parameters()
    pl <- simulate_patient(p, "placebo",     t_start = 0, t_end = 21)
    cb <- simulate_patient(p, "combination", t_start = 1, t_end = 21)
    expect_false(is.null(cb))

    m_pl <- extract_endpoints_from_simulation(pl, p, "HFRS")$mortality_prob
    m_cb <- extract_endpoints_from_simulation(cb, p, "HFRS")$mortality_prob
    expect_equal(100 * (m_pl - m_cb) / m_pl, 90.0, tolerance = 1e-2)
  })
})

test_that("placebo receives no drug and treated arms do", {
  with_project_root({
    p  <- get_parameters()
    pl <- simulate_patient(p, "placebo",     t_start = 0, t_end = 21)
    rb <- simulate_patient(p, "ribavirin",   t_start = 1, t_end = 21)
    fv <- simulate_patient(p, "favipiravir", t_start = 1, t_end = 21)

    expect_equal(max(pl$C_RBV), 0)
    expect_equal(max(pl$C_FAV), 0)
    expect_gt(max(rb$C_RBV), 0)
    expect_equal(max(rb$C_FAV), 0)
    expect_gt(max(fv$C_FAVI_RTP), 0)
    expect_equal(max(fv$C_RBV), 0)
  })
})

test_that("the dosing schedules are unchanged at their default duration", {
  # duration_days became a real control for the optimal-duration analysis;
  # the defaults must still reproduce the regimen used for the trial.
  r <- ribavirin_dosing_schedule()
  f <- favipiravir_dosing_schedule()
  expect_equal(nrow(r), 35L)
  expect_equal(max(r$time), 13.0)
  expect_equal(nrow(f), 30L)
  expect_equal(max(f$time), 17.5)

  # ... and a longer course must actually emit more doses.
  expect_gt(nrow(ribavirin_dosing_schedule(duration_days = 21)), nrow(r))
  expect_gt(nrow(favipiravir_dosing_schedule(duration_days = 21)), nrow(f))
})

test_that("synergy_combination() matches the rule used inside the ODE", {
  E1 <- 0.7; E2 <- 0.98; psi <- 0.1
  ode_rule <- 1 - (1 - E1) * (1 - E2) * (1 - psi * E1 * E2)
  expect_equal(synergy_combination(E1, E2, psi), ode_rule)
})
