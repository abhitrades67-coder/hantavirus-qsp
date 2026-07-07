#' Assemble journal-ready main figure files from outputs/
#'
#' Copies the four canonical model outputs into the manuscript submission
#' Figures/ folder, renaming to Figure_1..4.png. This is the final step of the
#' figure pipeline and is what makes the *submitted* figure files fully
#' reproducible from code (no manual cropping / renaming / re-tagging):
#'
#'   Rscript R/run_pipeline.R            # Figures 1 & 2 (viral kinetics, treatment window)
#'   Rscript R/plot_organ_heatmap.R      # Figure 3 (organ injury heatmap)
#'   Rscript R/plot_adaptive_heatmap.R   # Figure 4 (adaptive immunity heatmap)
#'   Rscript R/assemble_figures.R        # copy the four PNGs into the submission folder
#'
#' Run from the project root.

# Submission Figures/ folder (current target journal). Change `dest_dir` if the
# manuscript is retargeted to a different journal folder.
dest_dir <- file.path("manuscript", "Antiviral Research",
                      "Journal of Infection and Chemotherapy", "Figures")

# outputs/ source  ->  submission file name
fig_map <- rbind(
  c("outputs/viral_kinetics.png",            "Figure_1.png"),
  c("outputs/treatment_window.png",          "Figure_2.png"),
  c("outputs/organ_injury_heatmap.png",      "Figure_3.png"),
  c("outputs/adaptive_immunity_heatmap.png", "Figure_4.png")
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
