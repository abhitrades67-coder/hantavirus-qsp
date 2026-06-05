#' Cache-busting regeneration of the figure scripts that load cached CSVs.
#' These three scripts use `if (file.exists(cache)) load else simulate`, so the
#' 2026-06-03 recalibration was NOT picked up. Delete the stale caches, then
#' rerun so they re-simulate with the new parameters and rewrite the caches.
#' Run from project root: Rscript R/regen_cached_figs.R

stale_caches <- c(
  "outputs/organ_peak_data.csv",
  "outputs/organ_peak_placebo.csv",
  "outputs/biomarker_traj_data.csv",
  "outputs/adaptive_peak_data.csv",
  "outputs/adaptive_peak_placebo.csv"
)
for (f in stale_caches) {
  if (file.exists(f)) { file.remove(f); cat("Deleted stale cache:", f, "\n") }
}

scripts <- c(
  "R/plot_organ_heatmap.R",     # Figure 3 + fresh organ_peak_data.csv (organ % of placebo)
  "R/plot_biomarker_declutter.R", # Figure 6 (biomarker delta — fresh PLT)
  "R/plot_adaptive_heatmap.R"   # Figure 4 + fresh adaptive_peak_data.csv (adaptive % of placebo)
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
cat("\n\nCACHE-BUSTING REGENERATION COMPLETE.\n")
