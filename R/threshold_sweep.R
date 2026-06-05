# Quick test: check what fraction of patients cross the threshold
source("model/hantavirus_qsp.R")
source("model/parameters.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/virtual_population.R")
source("R/simulate_trial.R")

for (beta_val in c(7e-8, 8e-8, 1e-7)) {
  pars <- get_parameters()
  pars$beta <- beta_val
  pop <- generate_virtual_population(N = 100, seed = 42)
  
  n_est <- 0
  n_burn <- 0
  for (i in 1:nrow(pop)) {
    row <- pop[i, ]
    pars_i <- apply_vpop_to_model(row, pars)
    
    # Run pre-check
    y0 <- build_symptom_onset_state(pars_i)
    pre_times <- seq(0, 5, by = 1)
    pre_out <- deSolve::lsoda(y = y0, times = pre_times, func = hantavirus_qsp_ode,
                               parms = pars_i, atol = 1e-6, rtol = 1e-4, maxsteps = 5000)
    pre_out <- as.data.frame(pre_out)
    I_max <- max(pre_out$I, na.rm = TRUE)
    V_max <- max(pre_out$V, na.rm = TRUE)
    
    if (I_max < 1e5 && V_max < 1e6) {
      n_burn <- n_burn + 1
    } else {
      n_est <- n_est + 1
    }
  }
  cat(sprintf("beta=%.1e (R0=%.2f): established=%d (%.0f%%), burn-out=%d (%.0f%%)\n",
              beta_val, beta_val*pars$T_0*pars$p/pars$c,
              n_est, 100*n_est/100, n_burn, 100*n_burn/100))
}
