# Small serial test with 10 patients per arm
source("model/hantavirus_qsp.R")
source("model/parameters.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/virtual_population.R")
source("R/simulate_trial.R")
source("R/analysis.R")

pars <- get_parameters()
pop <- generate_virtual_population(N = 50, seed = 42)
set.seed(123)
idx <- sample(nrow(pop), 10)
pop_used <- pop[idx, ]

arms <- c("placebo", "ribavirin", "favipiravir", "combination")
endpoints <- list()

t0 <- Sys.time()
task_num <- 0
total_tasks <- length(arms) * nrow(pop_used)

for (arm in arms) {
  for (i in seq_len(nrow(pop_used))) {
    task_num <- task_num + 1
    row <- pop_used[i, ]
    pars_i <- apply_vpop_to_model(row, pars)
    t_treat <- if (arm == "placebo") 0 else row$time_to_treatment_days
    
    t1 <- Sys.time()
    sim <- tryCatch({
      simulate_patient(pars_i, arm, t_start = t_treat, t_end = 21, dt = 0.5)
    }, error = function(e) {
      cat(sprintf("  ERROR: %s\n", e$message))
      return(NULL)
    })
    elapsed <- difftime(Sys.time(), t1, units = "secs")
    
    if (is.null(sim)) {
      cat(sprintf("[%d/%d] %s pt %d: FAILED (%.1fs)\n", task_num, total_tasks, arm, i, elapsed))
      endpoints[[length(endpoints)+1]] <- data.frame(
        patient_id = row$patient_id, arm = arm, mortality_prob = NA,
        V_peak = NA, I_peak = NA, elapsed = elapsed, stringsAsFactors = FALSE
      )
    } else {
      ep <- extract_endpoints_from_simulation(sim, pars_i, syndrome = "HFRS")
      endpoints[[length(endpoints)+1]] <- data.frame(
        patient_id = row$patient_id, arm = arm, 
        mortality_prob = ep$mortality_prob,
        V_peak = max(sim$V), I_peak = max(sim$I), elapsed = elapsed,
        stringsAsFactors = FALSE
      )
      if (elapsed > 5) {
        cat(sprintf("[%d/%d] %s pt %d: mort=%.4f V=%.0f SLOW=%.1fs\n", 
                    task_num, total_tasks, arm, i, ep$mortality_prob, max(sim$V), elapsed))
      }
    }
  }
}
total <- difftime(Sys.time(), t0, units = "secs")

cat(sprintf("\nTotal: %.1fs for %d tasks\n", total, total_tasks))

ep_df <- do.call(rbind, endpoints)
for (arm in arms) {
  idx <- ep_df$arm == arm
  mort <- ep_df$mortality_prob[idx]
  n_na <- sum(is.na(mort))
  cat(sprintf("%-15s n=%d na=%d mean_mort=%.4f median=%.4f\n",
              arm, sum(idx), n_na, mean(mort, na.rm=TRUE), median(mort, na.rm=TRUE)))
}
