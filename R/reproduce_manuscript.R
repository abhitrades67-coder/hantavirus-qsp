#' Regenerate every output the manuscript and its supplement report.
#'
#' One entry point for the whole article, so that a reader who downloads the
#' archive can reproduce the figures and tables without working out which
#' script writes which file. Run from the project root:
#'
#'     Rscript R/reproduce_manuscript.R
#'
#' It runs for several hours. Each step is run in a separate R process so that
#' a failure in one does not leave state behind for the next, and the exit
#' status is non-zero if any step fails.
#'
#' After it finishes, build the documents:
#'
#'     python manuscript/Royal/_src/build_tables.py
#'     python manuscript/Royal/_src/build.py
#'     python manuscript/Royal/_src/check_numbers.py
#'     python manuscript/Royal/_src/finalise_counts.py
#'     python manuscript/Royal/_src/check_package.py
#'     python manuscript/Royal/_src/check_provenance.py
#'
#' To redraw the trajectory figures alone, after a label or axis change, run
#' R/regen_trajectory_figures.R instead. It reads the saved population back and
#' reproduces the same five patients, so figure 1 comes out byte for byte the
#' same and only the intended change appears.
#'
#' Several scripts cache their results and take the cached branch when the file
#' already exists. Set REGENERATE=TRUE below to delete those caches first; that
#' is what makes "reproduce" mean recompute rather than reprint.

REGENERATE <- TRUE

# Ordered by dependency. The comment on each line names what it writes that the
# article uses, so a reader can find the provenance of any number.
steps <- c(
  # Core pipeline: virtual population, parameter table, virtual trial,
  # treatment window, and figures 1 and 2.
  "R/run_pipeline.R",

  # Section 3.2 and supplementary table S5: what carries the output, plus the
  # section 3.1 description of the representative placebo patient.
  "R/mechanism_weights.R",

  # Section 3.5 and supplementary table S5: the ribavirin non-antiviral terms,
  # removed together and one at a time. Also supplementary figure S3.
  "R/ablation_analysis.R",

  # Figure 3 and the organ-injury percentages in section 3.3.
  "R/plot_organ_heatmap.R",

  # Supplementary figure S2 and the adaptive percentages in sections 3.4/3.6.
  "R/plot_adaptive_heatmap.R",

  # Figure 4 and supplementary table S6: start day crossed with duration.
  "R/duration_start_interaction.R",

  # Section 3.4: the on-treatment plateau and the rebound peaks, which the
  # summary grid above does not separate.
  "R/rebound_diagnostics.R",

  # Section 3.4 and supplementary table S7: how the result depends on humoral
  # strength, and supplementary table S8, the analytic invariance.
  "R/identifiability_and_humoral.R",

  # Figure 5 and table 3: what the calibration datum identifies.
  "R/identifiability_profile.R",

  # Section 3.7, supplementary table S2 and figure S4: prophylaxis.
  "R/preexposure_analysis.R",

  # Section 3.9: local sensitivity, then the global analysis behind
  # supplementary table S4 and figure S5.
  "R/sensitivity_analysis.R",
  "R/gsa_prcc.R",

  # Section 3.8: the joint version of the identifiability question, taken
  # from the designs the global analysis above has already run.
  "R/identifiability_joint.R",

  # Supplementary figure S6: model trajectories against independent data.
  "R/make_s7.R"
)

# Caches that would otherwise short-circuit the scripts above.
caches <- c(
  "outputs/organ_peak_data.csv",
  "outputs/organ_peak_placebo.csv",
  "outputs/adaptive_peak_data.csv",
  "outputs/adaptive_peak_placebo.csv",
  "outputs/duration_start_data.csv",
  "outputs/biomarker_traj_data.csv"
)

rscript <- file.path(R.home("bin"), "Rscript")
if (.Platform$OS.type == "windows") rscript <- paste0(rscript, ".exe")
if (!file.exists(rscript)) rscript <- "Rscript"

missing <- steps[!file.exists(steps)]
if (length(missing) > 0) {
  stop("these scripts are named above but are not in the repository: ",
       paste(missing, collapse = ", "))
}

if (REGENERATE) {
  for (f in caches) {
    if (file.exists(f)) {
      file.remove(f)
      cat("deleted cache:", f, "\n")
    }
  }
}

ok <- logical(length(steps))
names(ok) <- steps
t_all <- Sys.time()

for (i in seq_along(steps)) {
  s <- steps[i]
  cat(sprintf("\n=============== [%d/%d] %s ===============\n",
              i, length(steps), s))
  t0 <- Sys.time()
  rc <- system2(rscript, c(s))
  ok[s] <- identical(as.integer(rc), 0L)
  cat(sprintf("=============== %s: %s (%.1f min) ===============\n", s,
              if (ok[s]) "OK" else sprintf("FAILED, status %s", rc),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

cat("\n=============== SUMMARY ===============\n")
for (s in steps) cat(sprintf("  %-36s %s\n", s, if (ok[s]) "OK" else "FAILED"))
cat(sprintf("\ntotal %.1f min\n",
            as.numeric(difftime(Sys.time(), t_all, units = "mins"))))

if (any(!ok)) {
  cat(sprintf("\n%d of %d steps failed: %s\n", sum(!ok), length(steps),
              paste(steps[!ok], collapse = ", ")))
  quit(status = 1)
}
cat("\nAll steps complete. Now run the six Python scripts named at the top of\n")
cat("this file, in that order.\n")
