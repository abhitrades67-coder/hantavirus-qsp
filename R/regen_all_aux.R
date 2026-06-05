#' Regenerate all auxiliary analyses after the 2026-06-03 recalibration.
#' Each script re-sources the model and re-simulates, so it picks up the new
#' parameters (k_PLT_prod=25000, k_PLT_cons=0.002, K50_dialysis=70, L50_ecmo=400).
#' Run from the project root: Rscript R/regen_all_aux.R

scripts <- c(
  "R/plot_organ_heatmap.R",        # Figure 3 + organ_peak_data.csv (organ % of placebo)
  "R/plot_adaptive_heatmap.R",     # Figure 4 + adaptive_peak_data.csv (adaptive % of placebo)
  "R/plot_adaptive_peak_by_day.R", # Figure 5 (adaptive peak by day)
  "R/plot_biomarker_declutter.R",  # Figure 6 (biomarker delta — PLT changed)
  "R/preexposure_analysis.R",      # Figure 7 + preexposure_results.csv / summary (PEP)
  "R/ablation_analysis.R",         # Figures S3/S4 + ablation_summary.csv
  "R/sensitivity_analysis.R",      # tornado plots + sensitivity_summary.txt
  "R/uncertainty_quantification.R" # Figure S2 + uncertainty_intervals.csv
)

for (s in scripts) {
  cat("\n\n==================== RUNNING:", s, "====================\n")
  t0 <- Sys.time()
  ok <- tryCatch({ source(s, local = new.env()); TRUE },
                 error = function(e) { cat("ERROR in", s, ":", conditionMessage(e), "\n"); FALSE })
  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
  cat(sprintf("==================== %s: %s (%.2f min) ====================\n",
              s, if (ok) "OK" else "FAILED", dt))
}
cat("\n\nALL AUXILIARY ANALYSES COMPLETE.\n")
