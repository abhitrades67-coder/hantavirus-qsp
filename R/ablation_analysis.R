#' Representative ablation analyses for manuscript robustness checks.
#'
#' This script tests whether key manuscript conclusions depend on two model
#' assumptions: ribavirin non-antiviral endothelial/immunomodulatory effects and
#' favipiravir potency/adaptive timing assumptions.

suppressPackageStartupMessages({
  library(ggplot2)
})

source("model/parameters.R")
source("model/hantavirus_qsp.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/simulate_trial.R")

pars_base <- get_parameters()

extract_ablation_metrics <- function(sim, pars, syndrome = "HFRS") {
  pk <- extract_peak_biomarkers(sim)
  ep <- extract_endpoints_from_simulation(sim, pars, syndrome = syndrome)
  data.frame(
    V_AUC = compute_viral_AUC(sim),
    P_peak = pk$P_peak,
    K_peak = pk$K_peak,
    L_peak = pk$L_peak,
    C_pro_peak = pk$C_pro_peak,
    CD8_E_peak = max(sim$CD8_E, na.rm = TRUE),
    IgM_peak = max(sim$IgM, na.rm = TRUE),
    mortality_prob = ep$mortality_prob,
    stringsAsFactors = FALSE
  )
}

run_scenario <- function(label, arm, day, pars) {
  sim <- simulate_patient(pars, arm = arm, t_start = day, t_end = 21, dt = 0.1)
  if (is.null(sim)) stop(sprintf("Simulation failed: %s", label))
  out <- extract_ablation_metrics(sim, pars)
  out$scenario <- label
  out$arm <- arm
  out$treatment_day <- day
  out
}

scenarios <- list()

# Ribavirin non-antiviral effect ablation.
pars_no_rbv_nonantiviral <- pars_base
pars_no_rbv_nonantiviral$rbv_immuno_scale <- 0
pars_no_rbv_nonantiviral$rbv_endothelial_scale <- 0
pars_no_rbv_nonantiviral$rbv_recovery_scale <- 0

scenarios[["ribavirin_day1_default"]] <-
  run_scenario("ribavirin_day1_default", "ribavirin", 1, pars_base)
scenarios[["ribavirin_day1_no_nonantiviral"]] <-
  run_scenario("ribavirin_day1_no_nonantiviral", "ribavirin", 1, pars_no_rbv_nonantiviral)
scenarios[["combination_day1_default"]] <-
  run_scenario("combination_day1_default", "combination", 1, pars_base)
scenarios[["combination_day1_no_rbv_nonantiviral"]] <-
  run_scenario("combination_day1_no_rbv_nonantiviral", "combination", 1, pars_no_rbv_nonantiviral)

# Favipiravir paradox robustness checks.
fav_variants <- list(
  fav_default = pars_base,
  fav_EC50_half = within(pars_base, { EC50_FAV <- EC50_FAV * 0.5 }),
  fav_EC50_double = within(pars_base, { EC50_FAV <- EC50_FAV * 2 }),
  fav_Emax_090 = within(pars_base, { Emax_FAV <- 0.90 }),
  fav_CD8_fast = within(pars_base, { k_mat <- 0.20 }),
  fav_CD8_slow = within(pars_base, { k_mat <- 0.08 })
)

for (nm in names(fav_variants)) {
  scenarios[[paste0(nm, "_day1")]] <-
    run_scenario(paste0(nm, "_day1"), "favipiravir", 1, fav_variants[[nm]])
  scenarios[[paste0(nm, "_day2")]] <-
    run_scenario(paste0(nm, "_day2"), "favipiravir", 2, fav_variants[[nm]])
}

ablation <- do.call(rbind, scenarios)
ablation <- ablation[, c("scenario", "arm", "treatment_day", "V_AUC", "P_peak",
                         "K_peak", "L_peak", "C_pro_peak", "CD8_E_peak",
                         "IgM_peak", "mortality_prob")]

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
write.csv(ablation, "outputs/ablation_summary.csv", row.names = FALSE)

fav <- ablation[ablation$arm == "favipiravir", ]
fav$variant <- sub("_day[12]$", "", fav$scenario)

p1 <- ggplot(fav, aes(x = factor(treatment_day), y = K_peak, fill = factor(treatment_day))) +
  geom_col(width = 0.7, show.legend = FALSE) +
  facet_wrap(~ variant, scales = "free_y") +
  labs(x = "Treatment day", y = "Peak renal injury (K)",
       title = "Favipiravir day-1/day-2 paradox robustness checks") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey80"),
        plot.title = element_text(face = "bold", hjust = 0.5))

dir.create("manuscript/supplementary", showWarnings = FALSE, recursive = TRUE)
ggsave("manuscript/supplementary/FigureS3_ablation_analysis.png",
       p1, width = 9, height = 5.5, dpi = 300)

rbv <- ablation[grepl("ribavirin|combination", ablation$scenario), ]
rbv$nonantiviral <- ifelse(grepl("no", rbv$scenario), "Ablated", "Default")

p2 <- ggplot(rbv, aes(x = nonantiviral, y = K_peak, fill = arm)) +
  geom_col(position = "dodge", width = 0.7) +
  labs(x = "Ribavirin non-antiviral effects", y = "Peak renal injury (K)",
       fill = "Arm", title = "Ribavirin endothelial/immunomodulatory ablation") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey80"),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("manuscript/supplementary/FigureS4_ribavirin_ablation.png",
       p2, width = 7, height = 5, dpi = 300)

message("Saved: outputs/ablation_summary.csv")
message("Saved: manuscript/supplementary/FigureS3_ablation_analysis.png")
message("Saved: manuscript/supplementary/FigureS4_ribavirin_ablation.png")
