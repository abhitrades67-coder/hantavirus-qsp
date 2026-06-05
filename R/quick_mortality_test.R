# Quick single-patient mortality validation at beta=7e-8
source("model/hantavirus_qsp.R")
source("model/parameters.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/simulate_trial.R")

pars <- get_parameters()
pars$beta <- 7.0e-8
cat(sprintf("Testing with beta=%.1e, R0=%.2f\n", pars$beta, 
            pars$beta * pars$T_0 * pars$p / pars$c))

sim <- simulate_patient(pars, "placebo", t_start = 0, t_end = 60)

if (is.null(sim)) {
  cat("Simulation failed!\n")
} else {
  ep <- extract_endpoints_from_simulation(sim, pars, syndrome = "HFRS")
  cat(sprintf("beta=%.1e mortal=%.4f V_peak=%.0f I_peak=%.0f\n",
              pars$beta, ep$mortality_prob, max(sim$V), max(sim$I)))
  cat(sprintf("C_pro=%.4f K_max=%.2f L_max=%.2f P_max=%.2f\n",
              max(sim$C_pro), max(sim$K), max(sim$L), max(sim$P)))
}
