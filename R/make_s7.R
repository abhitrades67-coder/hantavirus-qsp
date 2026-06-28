#' Figure S7 — External corroboration of model trajectories
#'
#' Regenerates Supplementary Figure S7 from the QSP model. Four panels overlay
#' representative single-patient model time-courses with independent published
#' reference values that were NOT used in calibration (external consistency
#' checks, not formal predictive validation; cf. Supplementary Table S4):
#'
#'   A. Placebo viral load with the reported viraemic window shaded
#'      (PUUV RNA ~first 9 d, Pettersson 2014 [27]; DOBV peak ~1e7 copies/mL,
#'      viraemia ~16-30 d, Korva 2013 [28]).
#'   B. Virus / IgM / IgG (each scaled to its own maximum) — viral decline as
#'      antibody rises (Evander 2007 [29]).
#'   C. Placebo platelets with the reported 2nd-week nadir ~40,000/uL marked
#'      (dashed line + shaded window, Rasche 2004 [25]).
#'   D. Placebo vs favipiravir (day-1 start) viral load — early suppression and
#'      later rebound (Safronetz 2013 [14]; lethal ANDV hamster model showed a
#'      survival benefit only when favipiravir began before viraemia onset).
#'
#' Independent values are shown as annotated reference ranges, not digitized
#' curves; arbitrary-unit compartments are compared on timing and shape only.
#'
#' Run from project root:  Rscript R/make_s7.R
#' Output:                 outputs/figure_S7_external_validation.png

suppressPackageStartupMessages({
  library(deSolve)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

source("model/hantavirus_qsp.R")
source("model/parameters.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/simulate_trial.R")
source("R/analysis.R")

# --- Representative single-patient simulations (default parameter set) --------
pars <- get_parameters()

t_end <- 21
dt    <- 0.1

cat("Simulating placebo (representative patient)...\n")
pl  <- simulate_patient(pars, "placebo",     t_start = 0, t_end = t_end, dt = dt)
cat("Simulating favipiravir, treatment start day 1...\n")
fav <- simulate_patient(pars, "favipiravir", t_start = 1, t_end = t_end, dt = dt)

stopifnot(!is.null(pl), !is.null(fav))

# --- Shared aesthetics -------------------------------------------------------
col_placebo <- "#757575"  # grey   (placebo / virus)
col_igm     <- "#1f78b4"  # blue   (IgM)
col_igg     <- "#33a02c"  # green  (IgG; panel B has no red, CVD-safe here)
col_fav     <- "#E69F00"  # orange (favipiravir; colorblind-safe palette)
col_window  <- "#9ECAE1"  # light blue shaded reference windows
col_nadir   <- "#E69F00"  # platelet-nadir reference line

# Panel theme: left-aligned bold title + italic annotation subtitle
theme_s7 <- function(legend = "none") {
  theme_qsp() +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0, size = 15),
      plot.subtitle = element_text(face = "italic", hjust = 0, size = 9.5,
                                   colour = "grey30"),
      legend.position = legend,
      legend.title    = element_blank()
    )
}

# --- Panel A: placebo viral load, viraemic window shaded ---------------------
# Pin the log10 axis to decade bounds so the shaded band spans the full panel
# height (ymin/ymax = -Inf/Inf is NaN under log10 and would drop the rect).
yA_lo <- 10^floor(log10(min(pl$V[pl$V > 0])))
yA_hi <- 10^ceiling(log10(max(pl$V)))
p_a <- ggplot(pl, aes(time, V)) +
  annotate("rect", xmin = 0, xmax = 9, ymin = yA_lo, ymax = yA_hi,
           fill = col_window, alpha = 0.35) +
  geom_line(colour = col_placebo, linewidth = 1.1) +
  scale_y_log10(limits = c(yA_lo, yA_hi)) +
  labs(
    x = "Days post-symptom onset", y = "Viral RNA (copies/mL)",
    title    = "A  Viral load (placebo, model)",
    subtitle = paste0(
      "Independent: DOBV peak ~1e7 copies/mL, viraemia ~16-30 d (Korva 2013 [28]);\n",
      "PUUV RNA ~first 9 d (shaded; Pettersson 2014 [27])")
  ) +
  theme_s7()

# --- Panel B: virus vs antibody, each scaled to its maximum ------------------
pb_long <- data.frame(
  time  = pl$time,
  Virus = pl$V   / max(pl$V),
  IgM   = pl$IgM / max(pl$IgM),
  IgG   = pl$IgG / max(pl$IgG)
) |>
  pivot_longer(-time, names_to = "series", values_to = "rel") |>
  mutate(series = factor(series, levels = c("Virus", "IgM", "IgG")))

p_b <- ggplot(pb_long, aes(time, rel, colour = series)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = c(Virus = col_placebo, IgM = col_igm,
                                 IgG = col_igg)) +
  labs(
    x = "Days post-symptom onset", y = "Relative (fraction of max)",
    title    = "B  Virus vs antibody (placebo, model)",
    subtitle = paste0(
      "Independent: viraemia falls as IgM/IgG rise; fatal case = sustained\n",
      "high virus + absent antibody (Evander 2007 [29])")
  ) +
  theme_s7(legend = "top")

# --- Panel C: placebo platelets, 2nd-week nadir marked ----------------------
p_c <- ggplot(pl, aes(time, PLT)) +
  annotate("rect", xmin = 7, xmax = 14, ymin = -Inf, ymax = Inf,
           fill = col_window, alpha = 0.35) +
  geom_hline(yintercept = 40000, linetype = "dashed",
             colour = col_nadir, linewidth = 0.9) +
  geom_line(colour = col_placebo, linewidth = 1.1) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    x = "Days post-symptom onset", y = "Platelets (/uL)",
    title    = "C  Platelets (placebo, model)",
    subtitle = "Independent: nadir ~40,000/uL in 2nd week (dashed line / shaded; Rasche 2004 [25])"
  ) +
  theme_s7()

# --- Panel D: placebo vs favipiravir (day-1) viral load ---------------------
pd_df <- rbind(
  data.frame(time = pl$time,  V = pl$V,  arm = "Placebo"),
  data.frame(time = fav$time, V = fav$V, arm = "Favipiravir d1")
)
pd_df$arm <- factor(pd_df$arm, levels = c("Placebo", "Favipiravir d1"))

p_d <- ggplot(pd_df, aes(time, pmax(V, 1), colour = arm)) +
  geom_line(linewidth = 1.1) +
  scale_y_log10() +
  scale_colour_manual(values = c(Placebo = col_placebo,
                                 `Favipiravir d1` = col_fav)) +
  labs(
    x = "Days post-symptom onset", y = "Viral RNA (copies/mL)",
    title    = "D  Favipiravir antiviral effect (model)",
    subtitle = paste0(
      "Independent: hamster ANDV - oral favipiravir improved survival only if\n",
      "begun before viraemia onset (Safronetz 2013 [14])")
  ) +
  theme_s7(legend = "top")

# --- Assemble 2x2 and export (13x9 in @ 300 dpi = 3900x2700 px) --------------
g <- gridExtra::arrangeGrob(p_a, p_b, p_c, p_d, ncol = 2)
ggsave_safe("outputs/figure_S7_external_validation.png",
            plot = g, width = 13, height = 9, dpi = 300)
cat("Saved: outputs/figure_S7_external_validation.png\n")
