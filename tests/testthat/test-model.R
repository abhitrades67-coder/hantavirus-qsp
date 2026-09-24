make_state <- function(p, ...) {
  base <- c(T = p$T_0, I = p$I_0, V = p$V_0, NSs = 0, F_I = 0, NK = 0, F_II = 0,
            C_pro = 0, C_anti = 0, P = 0, PLT = p$PLT_0, K = 0, L = 0,
            CD8_N = 0.5, CD8_E = 0, CD4 = 0.5, IgM = 0.1, IgG = 0,
            C_RBV = 0, C_FAV_gut = 0, C_FAV = 0, C_FAVI_RTP = 0, Hgb_drop = 0)
  ov <- list(...)
  for (nm in names(ov)) base[[nm]] <- ov[[nm]]
  base
}

test_that("ODE returns 23 finite derivatives at the initial state", {
  p <- get_parameters()
  d <- hantavirus_qsp_ode(0, make_state(p), p)
  expect_type(d, "list")
  dydt <- d[[1]]
  expect_length(dydt, 23L)
  expect_true(all(is.finite(dydt)))
})

test_that("ODE stays finite under an active infection + drug state", {
  p <- get_parameters()
  y <- make_state(p, T = 5e5, I = 1e4, V = 1e6, NSs = 20, F_I = 80, NK = 5,
                  F_II = 10, C_pro = 300, C_anti = 50, P = 0.4, PLT = 1.2e5,
                  K = 20, L = 30, CD8_N = 5, CD8_E = 10, CD4 = 8, IgM = 5,
                  IgG = 3, C_RBV = 10, C_FAV_gut = 50, C_FAV = 5,
                  C_FAVI_RTP = 2, Hgb_drop = 1)
  dydt <- hantavirus_qsp_ode(0, y, p)[[1]]
  expect_length(dydt, 23L)
  expect_true(all(is.finite(dydt)))
})

test_that("clinical endpoints are valid probabilities and increase with injury", {
  p  <- get_parameters()
  lo <- compute_clinical_endpoints(K = 1,   L = 1,   C_pro = 10,   C_pro_max = 10,
                                   syndrome = "HFRS", pars = p)
  hi <- compute_clinical_endpoints(K = 100, L = 100, C_pro = 1000, C_pro_max = 1000,
                                   syndrome = "HFRS", pars = p)
  for (e in c("dialysis_prob", "ecmo_prob", "mortality_prob")) {
    expect_gte(lo[[e]], 0); expect_lte(lo[[e]], 1)
    expect_gte(hi[[e]], 0); expect_lte(hi[[e]], 1)
    expect_gte(hi[[e]], lo[[e]])
  }
})

test_that("the endpoint mapping is HFRS-only and refuses any other syndrome", {
  # The lung-weighted HCPS mapping was removed rather than shipped
  # uncalibrated. A request for it must fail loudly, not fall through to the
  # renal-weighted mapping and return a number that looks like an HCPS result.
  p <- get_parameters()
  expect_error(
    compute_clinical_endpoints(K = 0, L = 200, C_pro = 0, C_pro_max = 0,
                               syndrome = "HCPS", pars = p),
    "HFRS only"
  )
  expect_false(any(grepl("HCPS", names(p))))
})

test_that("mortality is renal-weighted, as an HFRS model requires", {
  p <- get_parameters()
  renal <- compute_clinical_endpoints(K = 200, L = 0, C_pro = 0, C_pro_max = 0,
                                      syndrome = "HFRS", pars = p)
  lung  <- compute_clinical_endpoints(K = 0, L = 200, C_pro = 0, C_pro_max = 0,
                                      syndrome = "HFRS", pars = p)
  expect_gt(renal$mortality_prob, lung$mortality_prob)
})
