source('R/simulate_trial.R')
source('R/pd_models.R')
source('R/pk_models.R')
source('model/hantavirus_qsp.R')
source('model/parameters.R')

pars <- get_parameters()

# Test simulate_patient directly
cat("Testing simulate_patient directly:\n\n")

# Test placebo
cat("=== Placebo ===\n")
s_p <- simulate_patient(pars, "placebo", t_start = 0, t_end = 10)
cat("V_max:", max(s_p$V), "\n")
cat("I_max:", max(s_p$I), "\n\n")

# Test ribavirin with t_start = 0
cat("=== Ribavirin (t_start = 0) ===\n")
s_r0 <- simulate_patient(pars, "ribavirin", t_start = 0, t_end = 10)
cat("V_max:", max(s_r0$V), "\n")
cat("I_max:", max(s_r0$I), "\n")
cat("C_RBV max:", max(s_r0$C_RBV), "\n\n")

# Test ribavirin with t_start = 1
cat("=== Ribavirin (t_start = 1) ===\n")
s_r <- simulate_patient(pars, "ribavirin", t_start = 1, t_end = 10)
cat("V_max:", max(s_r$V), "\n")
cat("I_max:", max(s_r$I), "\n")
cat("C_RBV max:", max(s_r$C_RBV), "\n\n")

# Check if the simulation was successful
cat("Simulation lengths:\n")
cat("Placebo:", nrow(s_p), "rows\n")
cat("Ribavirin D0:", nrow(s_r0), "rows\n")
cat("Ribavirin D1:", nrow(s_r), "rows\n")
