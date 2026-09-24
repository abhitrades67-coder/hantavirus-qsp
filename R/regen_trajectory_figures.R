#' Regenerate the trajectory figures without rerunning the whole pipeline.
#'
#' Figures 1 and S1, and the auxiliary trajectory plots, are drawn from the
#' same five sampled virtual patients simulated in every arm at every start
#' day. run_pipeline.R builds that data on its way through a much longer run,
#' so changing a plot label used to mean rerunning the trial. This reproduces
#' the trajectory set exactly: the population is read back from
#' outputs/virtual_population.csv and the same seed selects the same five
#' patients, so the figures come out identical apart from the change being made.
#'
#' Run from the project root: Rscript R/regen_trajectory_figures.R

suppressPackageStartupMessages({
  library(deSolve)
  library(ggplot2)
  library(gridExtra)
})

source("model/parameters.R")
source("model/hantavirus_qsp.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/simulate_trial.R")
source("R/virtual_population.R")
source("R/analysis.R")

pars <- get_parameters()

pop_path <- "outputs/virtual_population.csv"
if (!file.exists(pop_path)) {
  stop(pop_path, " is missing; run R/run_pipeline.R first")
}
pop <- read.csv(pop_path, stringsAsFactors = FALSE)

# Same seed and the same call as run_pipeline.R, so the same five patients.
set.seed(99)
plot_patients <- sample(pop$patient_id, 5)
cat("patients:", paste(plot_patients, collapse = ", "), "\n")

rows <- list()
add <- function(sim, arm, pid, day) {
  if (is.null(sim)) return(invisible(NULL))
  sim$arm <- arm
  sim$patient_id <- pid
  sim$treat_start_day <- as.integer(day)
  rows[[length(rows) + 1]] <<- sim
}

for (pid in plot_patients) {
  pars_i <- apply_vpop_to_model(pop[pop$patient_id == pid, ], pars)
  add(simulate_patient(pars_i, "placebo", t_start = 0, t_end = 21, dt = 0.5),
      "placebo", pid, 0L)
}
for (start_day in 1:7) {
  for (arm in c("ribavirin", "favipiravir", "combination")) {
    for (pid in plot_patients) {
      pars_i <- apply_vpop_to_model(pop[pop$patient_id == pid, ], pars)
      add(simulate_patient(pars_i, arm, t_start = start_day, t_end = 21,
                           dt = 0.5), arm, pid, start_day)
    }
  }
  cat("  start day", start_day, "done\n")
}
traj_data <- do.call(rbind, rows)
cat(nrow(traj_data), "rows from", length(rows), "simulations\n")

plot_viral_kinetics(
  traj_data[, c("time", "V", "arm", "patient_id", "treat_start_day")],
  outfile = "outputs/viral_kinetics.png")

plot_pk_profiles(
  traj_data[, c("time", "C_RBV", "C_FAV", "C_FAVI_RTP", "arm")],
  outfile = "outputs/pk_profiles.png")

plot_organ_injury(
  traj_data[, c("time", "K", "L", "arm", "patient_id")],
  outfile = "outputs/organ_injury.png")

plot_adaptive_immunity(
  traj_data[, c("time", "CD8_E", "CD4", "IgM", "IgG", "arm", "patient_id")],
  outfile = "outputs/adaptive_immunity.png")

cat("\nRegenerated: viral_kinetics.png (figure 1), pk_profiles.png (figure S1),\n")
cat("organ_injury.png and adaptive_immunity.png.\n")
