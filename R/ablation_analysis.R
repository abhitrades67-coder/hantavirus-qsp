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

# The three non-antiviral terms one at a time. Removing them together says how
# much they are worth in total; removing them singly says which one carries it,
# which is what section 3.5 and supplementary table S5 report.
single <- list(
  no_immuno       = within(pars_base, { rbv_immuno_scale <- 0 }),
  no_endothelial  = within(pars_base, { rbv_endothelial_scale <- 0 }),
  no_recovery     = within(pars_base, { rbv_recovery_scale <- 0 })
)
for (nm in names(single)) {
  lab <- paste0("ribavirin_day1_", nm)
  scenarios[[lab]] <- run_scenario(lab, "ribavirin", 1, single[[nm]])
}
for (nm in names(single)) {
  lab <- paste0("combination_day1_", nm)
  scenarios[[lab]] <- run_scenario(lab, "combination", 1, single[[nm]])
}

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

# Ribavirin non-antiviral terms, removed singly and together.
# An earlier version of this figure collapsed every ablated scenario onto one
# "Ablated" bar, so the three single-term ablations drew on top of each other
# and only the tallest was visible. Each condition now gets its own bar, and
# mortality is shown beside renal injury, which is what the legend describes.
LEVELS <- c(default = "Default",
            no_immuno = "No immunomodulatory",
            no_endothelial = "No endothelial",
            no_recovery = "No organ recovery",
            no_nonantiviral = "No non-antiviral terms",
            no_rbv_nonantiviral = "No non-antiviral terms")

rbv <- ablation[grepl("^(ribavirin|combination)_day1_", ablation$scenario), ]
key <- sub("^(ribavirin|combination)_day1_", "", rbv$scenario)
stopifnot(all(key %in% names(LEVELS)))
rbv$condition <- factor(unname(LEVELS[key]), levels = unique(unname(LEVELS)))

plot_df <- rbind(
  data.frame(arm = rbv$arm, condition = rbv$condition,
             metric = "Peak renal injury (AU)", value = rbv$K_peak),
  data.frame(arm = rbv$arm, condition = rbv$condition,
             metric = "Mortality (%)", value = 100 * rbv$mortality_prob))
plot_df$metric <- factor(plot_df$metric,
                         levels = c("Peak renal injury (AU)", "Mortality (%)"))

arm_colors <- c(ribavirin = "#1f78b4", combination = "#D55E00")

p2 <- ggplot(plot_df, aes(x = condition, y = value, fill = arm)) +
  geom_col(position = position_dodge(preserve = "single"), width = 0.7) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = arm_colors,
                    labels = c(ribavirin = "Ribavirin",
                               combination = "Combination"),
                    breaks = c("ribavirin", "combination")) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey80"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

ggsave("manuscript/supplementary/FigureS4_ribavirin_ablation.png",
       p2, width = 8, height = 4.6, dpi = 300)

message("Saved: outputs/ablation_summary.csv")
message("Saved: manuscript/supplementary/FigureS3_ablation_analysis.png")
message("Saved: manuscript/supplementary/FigureS4_ribavirin_ablation.png")
