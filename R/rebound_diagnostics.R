#' Trajectory diagnostics behind the withdrawal-rebound account in section 3.4.
#'
#' outputs/duration_start_data.csv records one summary row per cell of the
#' start-day by duration grid, and its V_peak column is the maximum over the
#' whole 42 days. For a day-1 start that maximum is the rebound peak, because
#' treatment suppresses the first one; for a day-2 start it is the first peak,
#' so the rebound has to be measured separately. The plateau values reached
#' while dosing continues are not in that file at all.
#'
#' This script measures those quantities directly, splitting each trajectory at
#' the end of dosing:
#'
#'   on treatment  the lowest viral load and infected-cell count reached while
#'                 dosing continues, which is what "suppresses without clearing"
#'                 refers to
#'   after dosing  the rebound peak of viral load, permeability and renal
#'                 injury, and the day each is reached
#'
#' Run from the project root: Rscript R/rebound_diagnostics.R

suppressPackageStartupMessages(library(deSolve))

source("model/parameters.R")
source("model/hantavirus_qsp.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/simulate_trial.R")
source("R/duration_start_interaction.R")   # simulate_patient_duration()

pars <- get_parameters()
T_END <- 42
DURATION <- 15          # the course used in the primary analysis

peak_after <- function(sim, col, from) {
  w <- sim$time >= from
  i <- which.max(sim[[col]][w])
  c(value = sim[[col]][w][i], day = sim$time[w][i])
}

rows <- list()
for (sd in c(1, 2)) {
  sim <- simulate_patient_duration(pars, "favipiravir", sd, DURATION,
                                   t_end = T_END, dt = 0.1)
  if (is.null(sim)) stop(sprintf("solver did not complete for a day-%d start", sd))

  stop_day <- sd + DURATION
  on <- sim$time >= sd & sim$time <= stop_day
  V_reb <- peak_after(sim, "V", stop_day)
  P_reb <- peak_after(sim, "P", stop_day)
  K_reb <- peak_after(sim, "K", stop_day)

  rows[[length(rows) + 1]] <- data.frame(
    arm = "favipiravir", start_day = sd, duration_days = DURATION,
    dosing_stops_day = stop_day,
    V_on_treatment_min = min(sim$V[on]),
    I_on_treatment_min = min(sim$I[on]),
    V_rebound_peak = V_reb[["value"]], t_V_rebound = V_reb[["day"]],
    P_rebound_peak = P_reb[["value"]], t_P_rebound = P_reb[["day"]],
    K_rebound_peak = K_reb[["value"]], t_K_rebound = K_reb[["day"]],
    T_end_frac = sim$T[nrow(sim)] / pars$T_0,
    stringsAsFactors = FALSE)
}

# Untreated reference: the first peak the rebound is compared against.
sim_pl <- simulate_patient_duration(pars, "placebo", 0, DURATION,
                                    t_end = T_END, dt = 0.1)
if (is.null(sim_pl)) stop("placebo simulation failed")
rows[[length(rows) + 1]] <- data.frame(
  arm = "placebo", start_day = NA_integer_, duration_days = NA_integer_,
  dosing_stops_day = NA_real_,
  V_on_treatment_min = NA_real_, I_on_treatment_min = NA_real_,
  V_rebound_peak = max(sim_pl$V), t_V_rebound = sim_pl$time[which.max(sim_pl$V)],
  P_rebound_peak = max(sim_pl$P), t_P_rebound = sim_pl$time[which.max(sim_pl$P)],
  K_rebound_peak = max(sim_pl$K), t_K_rebound = sim_pl$time[which.max(sim_pl$K)],
  T_end_frac = sim_pl$T[nrow(sim_pl)] / pars$T_0,
  stringsAsFactors = FALSE)

out <- do.call(rbind, rows)
dir.create("outputs", showWarnings = FALSE)
utils::write.csv(out, "outputs/rebound_diagnostics.csv", row.names = FALSE)

cat("\nFavipiravir, 15-day course, 42-day horizon\n")
for (i in seq_len(nrow(out))) {
  r <- out[i, ]
  cat(sprintf("  %-12s start %-3s  rebound V %.3g on day %.1f, P %.1f on day %.1f, K %.1f on day %.1f\n",
              r$arm, ifelse(is.na(r$start_day), "-", r$start_day),
              r$V_rebound_peak, r$t_V_rebound, r$P_rebound_peak,
              r$t_P_rebound, r$K_rebound_peak, r$t_K_rebound))
  if (!is.na(r$V_on_treatment_min)) {
    cat(sprintf("               on treatment, viral load held at or above %.3g and infected cells at or above %.3g\n",
                r$V_on_treatment_min, r$I_on_treatment_min))
  }
}
reb <- out$V_rebound_peak[out$arm == "favipiravir" & out$start_day == 1]
pl <- out$V_rebound_peak[out$arm == "placebo"]
cat(sprintf("\n  day-1 rebound peak is %.1f%% of the untreated peak, %.1f days later\n",
            100 * reb / pl,
            out$t_V_rebound[out$arm == "favipiravir" & out$start_day == 1] -
              out$t_V_rebound[out$arm == "placebo"]))
cat("\nSaved: outputs/rebound_diagnostics.csv\n")
