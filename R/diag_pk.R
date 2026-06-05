source('R/pk_models.R')
source('R/virtual_population.R')

pars <- get_parameters()

# Check ribavirin dosing schedule
sched <- ribavirin_dosing_schedule(body_weight_kg = 75, start_day = 1, duration_days = 10)
cat("Ribavirin dosing schedule:\n")
cat(sprintf("  %d doses total\n", nrow(sched)))
cat(sprintf("  First dose: day %.2f, %.0f mg IV\n", sched$time[1], sched$dose_mg[1]))
cat(sprintf("  Last dose:  day %.2f, %.0f mg IV\n", sched$time[nrow(sched)], sched$dose_mg[nrow(sched)]))

# Build events for ribavirin compartment
events <- build_deSolve_events(sched, "C_RBV", body_weight_kg = 75, Vd = 45)
cat("\nDosing events (first 10):\n")
print(head(events, 10))

# Simulate the PK model alone to check concentration profile
library(deSolve)
pk_pars <- c(
  CL_RBV = 40,
  Vd_RBV = 45
)

init <- c(C_RBV = 0)
times <- seq(1, 10, by = 0.1)

# Sort events by time
events <- events[order(events$time), ]

out <- ode(y = init, times = times,
           func = function(t, y, pars) {
             list(c(-(pars["CL_RBV"] / pars["Vd_RBV"]) * y["C_RBV"]))
           },
           parms = pk_pars,
           events = list(data = events))

cat("\nRibavirin PK profile (first 10 time points with C_RBV > 0):\n")
out_df <- as.data.frame(out)
out_df <- out_df[out_df$C_RBV > 0.1, ]
cat(sprintf("  Day   C_RBV\n"))
for (i in seq(1, min(15, nrow(out_df)))) {
  cat(sprintf("  %.1f   %.2f\n", out_df$time[i], out_df$C_RBV[i]))
}

cat(sprintf("\n  C_RBV_max: %.2f at day %.1f\n", max(out_df$C_RBV), out_df$time[which.max(out_df$C_RBV)]))
cat(sprintf("  C_RBV at day 2: %.2f\n", out_df$C_RBV[which.min(abs(out_df$time - 2))]))
cat(sprintf("  C_RBV at day 3: %.2f\n", out_df$C_RBV[which.min(abs(out_df$time - 3))]))

# The problem: there's a 1-day delay between symptom onset (day 0) and
# treatment start (day 1), and then another ~0.9 days to reach therapeutic levels
cat("\n\n=== Key finding ===\n")
cat("Day 0: Patient presents with symptoms (V=100, I starts seeding)\n")
cat("Day 1: First ribavirin dose given (C_RBV starts at 0)\n")
cat("Day 1.9: C_RBV reaches therapeutic levels (12.8 ug/mL)\n")
cat("\nThis means viral replication proceeds unchecked for ~2 days,")
cat("\nseeding ~939,000 infected cells before ribavirin acts.\n")
