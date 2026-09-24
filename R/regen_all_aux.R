#' Regenerate all auxiliary analyses after the 2026-06-03 recalibration.
#' Each script re-sources the model and re-simulates, so it picks up the new
#' parameters (k_PLT_prod=25000, k_PLT_cons=0.002, K50_dialysis=70, L50_ecmo=400).
#' Run from the project root: Rscript R/regen_all_aux.R

scripts <- c(
  "R/plot_organ_heatmap.R",        # Figure 3 + organ_peak_data.csv (organ % of placebo)
  "R/plot_adaptive_heatmap.R",     # Figure 4 + adaptive_peak_data.csv (adaptive % of placebo)
  "R/plot_adaptive_peak_by_day.R", # Figure 5 (adaptive peak by day)
  "R/plot_adaptive_by_day.R",      # adaptive_immunity_by_day.png (trajectories by start day)
  "R/plot_biomarker_declutter.R",  # Figure 6 (biomarker delta — PLT changed)
  "R/preexposure_analysis.R",      # Figure 7 + preexposure_results.csv / summary (PEP)
  "R/ablation_analysis.R",         # Figures S3/S4 + ablation_summary.csv
  "R/sensitivity_analysis.R",      # tornado plots + sensitivity_summary.txt
  "R/uncertainty_quantification.R", # Figure S2 + uncertainty_intervals.csv
  # Population-level check behind the figure S4 robustness paragraph. This was
  # driven by nothing until 2026-09-19 and its output had gone three months
  # stale while the supplement quoted from it.
  "R/adaptive_population_check.R"  # adaptive_population_check.csv
)

# These two use sink() and fail with "invalid connection" when sourced inside
# another script's environment; README.md states they must be run as standalone
# Rscript invocations. Run them as separate processes instead of source().
standalone_scripts <- c(
  "R/sensitivity_analysis.R",
  "R/uncertainty_quantification.R"
)

# Rscript may not be on PATH; resolve it from the running R installation.
rscript_bin <- file.path(R.home("bin"), "Rscript")
if (.Platform$OS.type == "windows") rscript_bin <- paste0(rscript_bin, ".exe")
if (!file.exists(rscript_bin)) rscript_bin <- "Rscript"  # last resort: rely on PATH

# The cache-backed figure scripts take their `if (file.exists(cache))` branch and
# simply re-emit the stale cached numbers. Deleting the caches first is what makes
# "regenerate" actually regenerate. Same file list as R/regen_cached_figs.R.
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

ok_flags <- logical(0)

for (s in scripts) {
  cat("\n\n==================== RUNNING:", s, "====================\n")
  t0 <- Sys.time()
  if (s %in% standalone_scripts) {
    ok <- tryCatch({
      rc <- system2(rscript_bin, c(s))
      if (!identical(as.integer(rc), 0L)) {
        cat("ERROR in", s, ": Rscript exited with status", rc, "\n")
      }
      identical(as.integer(rc), 0L)
    }, error = function(e) { cat("ERROR in", s, ":", conditionMessage(e), "\n"); FALSE })
  } else {
    ok <- tryCatch({ source(s, local = new.env()); TRUE },
                   error = function(e) { cat("ERROR in", s, ":", conditionMessage(e), "\n"); FALSE })
  }
  ok_flags[s] <- ok
  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
  cat(sprintf("==================== %s: %s (%.2f min) ====================\n",
              s, if (ok) "OK" else "FAILED", dt))
}

cat("\n\n==================== SUMMARY ====================\n")
for (s in names(ok_flags)) {
  cat(sprintf("  %-32s %s\n", s, if (ok_flags[[s]]) "OK" else "FAILED"))
}

failed <- names(ok_flags)[!ok_flags]
if (length(failed) > 0) {
  cat(sprintf("\n%d of %d auxiliary analyses FAILED: %s\n",
              length(failed), length(ok_flags), paste(failed, collapse = ", ")))
  quit(status = 1)
}
cat("\nALL AUXILIARY ANALYSES COMPLETE.\n")
