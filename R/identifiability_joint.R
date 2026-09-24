#' The joint version of the identifiability question.
#'
#' R/identifiability_profile.R walks one parameter at a time with the rest held
#' at their calibrated values. That is a conditional region, and a reader is
#' entitled to ask whether varying parameters together would close it. This
#' applies the same criterion to the global sensitivity designs, in which all 32
#' uncertain parameters vary simultaneously: keep the designs whose placebo
#' mortality the calibration datum cannot reject, and report the spread of the
#' predicted treatment effect across them.
#'
#' The criterion is the same in both cases. The datum is 7 deaths among 75
#' placebo patients; a parameter set is kept when its predicted placebo
#' mortality falls inside the exact binomial 95 per cent interval of that
#' observation, which is the set of values a test on the datum cannot reject.
#'
#' Designs in which infection does not establish are excluded, because the
#' relative risk reduction is undefined there, not because the datum rejects
#' them.
#'
#' Run from the project root: Rscript R/identifiability_joint.R

lhs_path <- "outputs/gsa_lhs_samples.csv"
if (!file.exists(lhs_path)) {
  stop(lhs_path, " is missing; run R/gsa_prcc.R first")
}
lhs <- read.csv(lhs_path, stringsAsFactors = FALSE)

ci <- binom.test(7, 75)$conf.int
lo <- 100 * ci[1]
hi <- 100 * ci[2]

# The column is written as 0/1 by R/gsa_prcc.R; accept the logical and
# character spellings too rather than silently selecting nothing.
flag <- lhs$established
est <- lhs[as.character(flag) %in% c("1", "TRUE", "true"), ]
est <- est[!is.na(est$combo_RRR), ]
mort <- 100 * est$placebo_mortality
keep <- est[mort >= lo & mort <= hi, ]

if (nrow(keep) == 0) {
  stop("no designs were selected; check how outputs/gsa_lhs_samples.csv ",
       "spells its established flag and its mortality scale")
}

rrr <- 100 * keep$combo_RRR
qs <- stats::quantile(rrr, c(0.05, 0.25, 0.5, 0.75, 0.95), names = FALSE)

out <- data.frame(
  quantity = c("ci_lower_pct", "ci_upper_pct", "n_designs", "n_established",
               "n_consistent", "n_rejected",
               "rrr_min", "rrr_p05", "rrr_q1", "rrr_median", "rrr_q3",
               "rrr_p95", "rrr_max"),
  value = c(lo, hi, nrow(lhs), nrow(est), nrow(keep), nrow(est) - nrow(keep),
            min(rrr), qs[1], qs[2], qs[3], qs[4], qs[5], max(rrr)),
  stringsAsFactors = FALSE)

dir.create("outputs", showWarnings = FALSE)
utils::write.csv(out, "outputs/identifiability_joint.csv", row.names = FALSE)

cat(sprintf("\nExact binomial 95%% interval on 7/75: %.2f to %.2f per cent\n", lo, hi))
cat(sprintf("  %d designs, %d establishing, %d consistent with the datum (%.1f%%)\n",
            nrow(lhs), nrow(est), nrow(keep), 100 * nrow(keep) / nrow(est)))
cat(sprintf("  the datum rejects %d of the %d establishing designs\n",
            nrow(est) - nrow(keep), nrow(est)))
cat(sprintf("  day-1 combination RRR across the consistent designs:\n"))
cat(sprintf("    %.1f to %.1f per cent, median %.1f, interquartile range %.1f to %.1f\n",
            min(rrr), max(rrr), qs[3], qs[2], qs[4]))
cat("\nSaved: outputs/identifiability_joint.csv\n")
