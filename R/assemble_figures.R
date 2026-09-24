#' Assemble journal-ready main figure files from outputs/
#'
#' Copies the five canonical model outputs into the manuscript submission
#' Figures/ folder, renaming to Figure_1..5.png. This is the final step of the
#' figure pipeline and is what makes the *submitted* figure files fully
#' reproducible from code (no manual cropping / renaming / re-tagging):
#'
#'   Rscript R/run_pipeline.R                # Figures 1 & 2 (viral kinetics, treatment window)
#'   Rscript R/plot_organ_heatmap.R          # Figure 3 (organ injury heatmap)
#'   Rscript R/plot_adaptive_heatmap.R       # Figure 4 (adaptive immunity heatmap)
#'   Rscript R/duration_start_interaction.R  # Figure 5 (start day x course duration)
#'   Rscript R/assemble_figures.R            # copy the five PNGs into the submission folder
#'
#' This script covers the MAIN-TEXT figures only. The complete submission
#' package -- these five plus the seven electronic supplementary figures, and
#' the .docx conversions -- is assembled by
#'   python "manuscript/Royal Society/_src/build.py"
#' which is the canonical step and carries the full figure mapping. Use this
#' script when only the main figures need refreshing, or when Python is
#' unavailable.
#'
#' Run from the project root.

# Submission Figures/ folder (current target journal: J. R. Soc. Interface).
# Change `dest_dir` if the manuscript is retargeted to a different journal
# folder, and keep it in step with FIGURES in _src/build.py.
dest_dir <- file.path("manuscript", "Royal Society", "Figures")

# outputs/ source  ->  submission file name
fig_map <- rbind(
  c("outputs/viral_kinetics.png",             "Figure_1.png"),
  c("outputs/treatment_window.png",           "Figure_2.png"),
  c("outputs/organ_injury_heatmap.png",       "Figure_3.png"),
  c("outputs/adaptive_immunity_heatmap.png",  "Figure_4.png"),
  c("outputs/duration_start_interaction.png", "Figure_5.png")
)

if (!dir.exists(dest_dir)) {
  dir.create(dest_dir, recursive = TRUE)
  cat("Created:", dest_dir, "\n")
}

missing <- fig_map[!file.exists(fig_map[, 1]), 1]
if (length(missing) > 0) {
  stop("Missing source figure(s) — run the figure scripts first:\n  ",
       paste(missing, collapse = "\n  "))
}

cat("Assembling submission figures into:\n  ", dest_dir, "\n\n", sep = "")
for (i in seq_len(nrow(fig_map))) {
  src  <- fig_map[i, 1]
  dest <- file.path(dest_dir, fig_map[i, 2])
  ok   <- file.copy(src, dest, overwrite = TRUE)
  cat(sprintf("  %-42s -> %-14s [%s]\n", src, fig_map[i, 2],
              if (ok) "OK" else "FAILED"))
}
cat("\nFigure assembly complete.\n")
