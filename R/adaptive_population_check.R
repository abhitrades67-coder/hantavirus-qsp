#' Population-level robustness check for adaptive-immunity treatment effects.
#'
#' Confirms (at the 200-patient population level) the representative-patient
#' finding that CD8+/CD4+ peak magnitudes are treatment-invariant while IgM is
#' modulated. The same 200 patients are simulated in every arm (paired design),
#' so each treated patient's peak is expressed as % of that SAME patient's
#' placebo peak. Output: outputs/adaptive_population_check.csv
#'
#' Usage: Rscript R/adaptive_population_check.R   (run from project root)

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

source("model/hantavirus_qsp.R")
source("model/parameters.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/virtual_population.R")
source("R/simulate_trial.R")

pars  <- get_parameters()
pop   <- generate_virtual_population(N = 1000, seed = 42)
trial <- run_virtual_trial(
  pop, pars,
  arms = c("placebo", "ribavirin", "favipiravir", "combination"),
  n_patients = 200, t_end = 21, dt = 0.5, seed = 123
)
ep <- trial$endpoints

comps <- c("CD8_E_peak", "CD4_peak", "IgM_peak", "IgG_peak")

# Placebo peak for each patient (paired reference)
placebo <- ep |> dplyr::filter(arm == "placebo") |>
  dplyr::select(patient_id, dplyr::all_of(comps))
names(placebo)[-1] <- paste0(names(placebo)[-1], "_pl")

m <- ep |> dplyr::filter(arm != "placebo") |>
  dplyr::left_join(placebo, by = "patient_id")
m$arm <- factor(m$arm, levels = c("ribavirin", "favipiravir", "combination"))

summ <- list()
for (cp in comps) {
  ratio <- 100 * m[[cp]] / m[[paste0(cp, "_pl")]]
  ok <- is.finite(ratio)
  d  <- data.frame(arm = m$arm[ok], ratio = ratio[ok])
  s  <- d |> dplyr::group_by(arm) |>
    dplyr::summarise(
      n           = dplyr::n(),
      median_pct  = median(ratio),
      p05_pct     = quantile(ratio, 0.05),
      p95_pct     = quantile(ratio, 0.95),
      max_abs_dev = max(abs(ratio - 100)),
      .groups = "drop"
    )
  s$compartment <- cp
  summ[[cp]] <- s
}
out <- do.call(rbind, summ)[, c("compartment", "arm", "n", "median_pct",
                                 "p05_pct", "p95_pct", "max_abs_dev")]
write.csv(out, "outputs/adaptive_population_check.csv", row.names = FALSE)
print(out, row.names = FALSE)

tcell <- out[out$compartment %in% c("CD8_E_peak", "CD4_peak"), ]
cat(sprintf(
  "\nT-cell peaks: median %% of placebo across treated arms = %.1f-%.1f%%; max per-patient deviation = %.2f%%\n",
  min(tcell$median_pct), max(tcell$median_pct), max(tcell$max_abs_dev)))
igm <- out[out$compartment == "IgM_peak", ]
cat("IgM median %% of placebo by arm: ",
    paste(sprintf("%s=%.0f%%", igm$arm, igm$median_pct), collapse = ", "), "\n")
