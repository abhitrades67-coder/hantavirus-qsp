#' Mini-population calibration for platelet nadir + dialysis/ECMO thresholds.
#'
#' Runs a serial placebo cohort with the platelet bug-fix applied, collects the
#' K_peak / L_peak / PLT_nadir / mortality distributions, and reports the mean
#' dialysis/ECMO rate that candidate decoupled thresholds would produce
#' (computed analytically from the per-patient peaks, so no re-simulation).

suppressWarnings(suppressMessages({
  source("model/hantavirus_qsp.R")
  source("model/parameters.R")
  source("R/pk_models.R")
  source("R/pd_models.R")
  source("R/virtual_population.R")
  source("R/simulate_trial.R")
}))

OVERRIDES <- list(k_PLT_prod = 25000, k_PLT_cons = 0.002)
N_SAMPLE  <- 80

pars0 <- get_parameters()
pop <- generate_virtual_population(N = 300, seed = 42)
set.seed(123)
pop_used <- pop[sample(nrow(pop), N_SAMPLE), ]

K <- L <- PLTn <- mort <- syn <- rep(NA_real_, nrow(pop_used))
syn <- character(nrow(pop_used))
for (i in seq_len(nrow(pop_used))) {
  row <- pop_used[i, ]
  pars_i <- apply_vpop_to_model(row, pars0)
  for (nm in names(OVERRIDES)) pars_i[[nm]] <- OVERRIDES[[nm]]
  sim <- tryCatch(simulate_patient(pars_i, "placebo", t_start = 0, t_end = 21, dt = 0.5),
                  error = function(e) NULL)
  if (is.null(sim)) next
  pk <- extract_peak_biomarkers(sim)
  ep <- extract_endpoints_from_simulation(sim, pars_i, syndrome = row$syndrome)
  K[i] <- pk$K_peak; L[i] <- pk$L_peak; PLTn[i] <- pk$PLT_nadir
  mort[i] <- ep$mortality_prob; syn[i] <- row$syndrome
}
ok <- !is.na(K)
K <- K[ok]; L <- L[ok]; PLTn <- PLTn[ok]; mort <- mort[ok]

cat(sprintf("\nPlacebo cohort N=%d (%d HFRS / %d HCPS)\n", length(K),
            sum(syn[ok]=="HFRS"), sum(syn[ok]=="HCPS")))
cat(sprintf("PLT_nadir: median=%.0f  mean=%.0f  IQR=[%.0f, %.0f]\n",
            median(PLTn), mean(PLTn), quantile(PLTn,.25), quantile(PLTn,.75)))
cat(sprintf("Mortality: mean=%.1f%%  median=%.1f%%\n", 100*mean(mort), 100*median(mort)))
cat(sprintf("K_peak: median=%.1f  mean=%.1f   L_peak: median=%.1f  mean=%.1f\n",
            median(K), mean(K), median(L), mean(L)))

cat("\n=== Mean DIALYSIS rate = mean(K/(K+K50_d)) ===\n")
for (k50 in c(30, 50, 70, 90, 120)) {
  cat(sprintf("  K50_dialysis=%-4d -> %.1f%%\n", k50, 100*mean(K/(K+k50))))
}
cat("\n=== Mean ECMO rate = mean(L/(L+L50_e)) ===\n")
for (l50 in c(100, 200, 300, 400, 500)) {
  cat(sprintf("  L50_ecmo=%-4d -> %.1f%%\n", l50, 100*mean(L/(L+l50))))
}
