#!/usr/bin/env Rscript
# ============================================================================
# Global Sensitivity Analysis (LHS + PRCC) — Hantavirus QSP Model
# ============================================================================
# Latin Hypercube Sampling over the uncertain mechanistic/PD parameters
# (log-uniform, confidence-tiered ranges) with Partial Rank Correlation
# Coefficients for the key model outputs, evaluated on the representative
# HFRS treatment-window patient:
#   - placebo mortality probability
#   - day-1 combination relative risk reduction (RRR) vs placebo
#   - peak viral load (placebo)
#
# Outputs:
#   outputs/gsa_lhs_samples.csv   raw LHS design + outputs
#   outputs/gsa_prcc.csv          PRCC per parameter x output
#   outputs/gsa_prcc_tornado.png  PRCC tornado figure (Supplementary Fig S8)
#   outputs/gsa_summary.txt       ranked text summary
#
# Usage (from project root): Rscript R/gsa_prcc.R
# ============================================================================

suppressPackageStartupMessages({
  library(deSolve); library(parallel); library(ggplot2); library(dplyr)
})

source("model/parameters.R")
source("model/hantavirus_qsp.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/simulate_trial.R")

set.seed(2026)
N      <- 1000      # LHS sample size
t_end  <- 21
dt     <- 0.5
syndrome <- "HFRS"
root   <- normalizePath(getwd())

# --- Uncertain parameters and multiplicative (log-uniform) ranges -----------
# factor range applied to the nominal value; tiered by confidence/role.
ranges <- list(
  beta=c(.5,2), p=c(.5,2), c=c(.7,1.5),
  delta_natural=c(.33,3), delta_NK=c(.33,3), delta_CD8=c(.33,3),
  k_FI=c(.5,2), E_max_IFN=c(.7,1.2), EC50_F=c(.5,2),
  k_CI=c(.5,2), k_CF2=c(.5,2), k_CC=c(.6,1.6), K_CC=c(.6,1.6),
  d_C=c(.7,1.5), I_gate=c(.5,2),
  k_PV=c(.5,2), k_PC=c(.5,2), d_P=c(.6,1.6),
  k_KP=c(.5,2), k_KC=c(.5,2), d_K=c(.6,1.6),
  k_LP=c(.5,2), k_LC=c(.5,2), d_L=c(.6,1.6),
  Emax_RBV=c(.7,1.3), EC50_RBV=c(.5,2),
  Emax_FAV=c(.85,1.02), EC50_FAV=c(.5,2), psi=c(.2,3),
  k_form=c(.5,2), k_elim_RTP=c(.5,2), CL_RBV=c(.6,1.6)
)
pnames <- names(ranges)
k <- length(pnames)
emax_cap <- c(Emax_RBV=0.95, Emax_FAV=0.999, E_max_IFN=0.95)

# --- Latin Hypercube design in [0,1] ----------------------------------------
lhs01 <- function(N, k) {
  m <- matrix(0, N, k)
  for (j in seq_len(k)) m[, j] <- (sample(N) - runif(N)) / N
  m
}
U <- lhs01(N, k)
colnames(U) <- pnames

# Map to log-uniform multiplicative factors -> absolute parameter values
pars0 <- get_parameters()
design <- as.data.frame(matrix(0, N, k)); names(design) <- pnames
for (j in seq_len(k)) {
  lo <- log(ranges[[pnames[j]]][1]); hi <- log(ranges[[pnames[j]]][2])
  fac <- exp(lo + U[, j] * (hi - lo))
  val <- pars0[[pnames[j]]] * fac
  if (pnames[j] %in% names(emax_cap)) val <- pmin(val, emax_cap[[pnames[j]]])
  design[[j]] <- val
}

cat(sprintf("GSA: N=%d LHS samples x %d parameters; 4 sims/sample\n", N, k))

# --- Worker: one LHS sample -> representative-patient outputs ----------------
n_cores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterExport(cl, c("design","pnames","pars0","t_end","dt","syndrome","root"),
                        envir = environment())

rows <- parallel::parLapply(cl, seq_len(N), function(i) {
  setwd(root)
  source("model/hantavirus_qsp.R", local = TRUE)
  source("model/parameters.R", local = TRUE)
  source("R/pk_models.R", local = TRUE)
  source("R/pd_models.R", local = TRUE)
  source("R/simulate_trial.R", local = TRUE)

  pars <- pars0
  for (nm in pnames) pars[[nm]] <- design[[nm]][i]

  mort <- function(arm, t_start) {
    s <- tryCatch(simulate_patient(pars, arm, t_start = t_start, t_end = t_end, dt = dt),
                  error = function(e) NULL)
    if (is.null(s)) return(c(m = NA, vpk = NA, vauc = NA))
    ep <- extract_endpoints_from_simulation(s, pars, syndrome = syndrome)
    pk <- extract_peak_biomarkers(s)
    c(m = ep$mortality_prob, vpk = pk$V_peak, vauc = compute_viral_AUC(s))
  }
  pl <- mort("placebo", 0)
  rb <- mort("ribavirin", 1)
  fv <- mort("favipiravir", 1)
  cb <- mort("combination", 1)
  rrr <- function(mt) if (!is.na(pl["m"]) && pl["m"] > 1e-4) (pl["m"] - mt) / pl["m"] else NA

  data.frame(
    placebo_mortality = unname(pl["m"]),
    rbv_mortality     = unname(rb["m"]),
    fav_mortality     = unname(fv["m"]),
    combo_mortality   = unname(cb["m"]),
    combo_RRR         = unname(rrr(cb["m"])),
    rbv_RRR           = unname(rrr(rb["m"])),
    placebo_Vpeak     = unname(pl["vpk"]),
    placebo_VAUC      = unname(pl["vauc"]),
    combo_VAUC        = unname(cb["vauc"])
  )
})
res <- do.call(rbind, rows)
samples <- cbind(design, res)

# --- Infection establishment ------------------------------------------------
# The RRR is NA whenever placebo mortality <= 1e-4 (see rrr() above), i.e. in
# exactly those designs where the infection never establishes. Every RRR
# statistic is therefore CONDITIONAL on establishment, while the placebo
# mortality summary quoted beside it is not. Make the conditioning an explicit,
# recorded variable rather than an implicit side effect of the NA rule.
ESTABLISH_MORT_THRESHOLD <- 1e-4
samples$established <- as.numeric(is.finite(samples$placebo_mortality) &
                                    samples$placebo_mortality > ESTABLISH_MORT_THRESHOLD)

n_valid <- sum(is.finite(res$placebo_mortality))
n_established <- sum(samples$established == 1, na.rm = TRUE)

write.csv(samples, "outputs/gsa_lhs_samples.csv", row.names = FALSE)
cat(sprintf("Completed %d samples (%d valid placebo sims)\n", N, n_valid))
cat(sprintf("Infection established in %d/%d designs (%.1f%%); %d (%.1f%%) excluded from RRR statistics\n",
            n_established, N, 100 * n_established / N,
            N - n_established, 100 * (N - n_established) / N))

# --- PRCC -------------------------------------------------------------------
prcc_vec <- function(X, y) {
  ok <- is.finite(y) & apply(is.finite(X), 1, all)
  X <- X[ok, , drop = FALSE]; y <- y[ok]
  Rx <- apply(X, 2, rank); ry <- rank(y)
  kk <- ncol(X); out <- numeric(kk)
  for (j in seq_len(kk)) {
    oth <- Rx[, -j, drop = FALSE]
    rj  <- residuals(lm(Rx[, j] ~ oth))
    ryj <- residuals(lm(ry ~ oth))
    out[j] <- cor(rj, ryj)
  }
  list(prcc = out, n = sum(ok))
}

# --- Drug parameters ---------------------------------------------------------
# These cannot act in the placebo arm. Leaving them in the placebo PRCC input
# only consumes degrees of freedom and produces spurious "drivers" (e.g.
# Emax_RBV appearing as a top driver of PLACEBO viral load). Named vector so the
# exclusion list is auditable; applied ONLY to placebo-arm outcomes.
DRUG_PARAMS <- c(
  Emax_RBV                = "ribavirin PD (Emax)",
  EC50_RBV                = "ribavirin PD (EC50)",
  gamma_RBV               = "ribavirin PD (Hill)",
  Emax_FAV                = "favipiravir PD (Emax)",
  EC50_FAV                = "favipiravir PD (EC50)",
  gamma_FAV               = "favipiravir PD (Hill)",
  psi                     = "drug interaction (synergy) term",
  rbv_immuno_scale        = "ribavirin immunomodulation ablation control",
  rbv_endothelial_scale   = "ribavirin endothelial-protection ablation control",
  rbv_recovery_scale      = "ribavirin recovery-enhancement ablation control",
  Vd_RBV                  = "ribavirin PK (volume)",
  CL_RBV                  = "ribavirin PK (clearance)",
  ka_FAV                  = "favipiravir PK (absorption)",
  Vd_FAV                  = "favipiravir PK (volume)",
  CL_FAV                  = "favipiravir PK (clearance)",
  Km_FAV                  = "favipiravir PK (saturable elimination)",
  k_form                  = "ribavirin triphosphate formation",
  k_elim_RTP              = "ribavirin triphosphate elimination",
  k_hgb                   = "ribavirin haemolytic anaemia",
  d_hgb                   = "haemoglobin recovery"
)

outcomes <- c(placebo_mortality = "Placebo mortality",
              combo_RRR        = "Day-1 combination RRR",
              placebo_Vpeak    = "Peak viral load (placebo)",
              established      = "Infection establishment (binary)")
# Outcomes generated entirely by the placebo arm: drug parameters are dropped.
placebo_outcomes <- c("placebo_mortality", "placebo_Vpeak", "established")

Xmat <- as.matrix(design)
prcc_tbl <- do.call(rbind, lapply(names(outcomes), function(oc) {
  use <- if (oc %in% placebo_outcomes) setdiff(pnames, names(DRUG_PARAMS)) else pnames
  pr  <- prcc_vec(Xmat[, use, drop = FALSE], samples[[oc]])
  data.frame(parameter = use, outcome = outcomes[[oc]],
             # "established" is a 0/1 outcome: ranking it produces midranks, so
             # the same partial-rank machinery yields a partial rank
             # point-biserial correlation. Labelled distinctly so it is not
             # read as a PRCC on a continuous outcome.
             method = if (oc == "established")
                        "partial rank point-biserial (binary outcome)" else "PRCC",
             prcc = pr$prcc, n = pr$n, k_used = length(use),
             # two-sided significance threshold ~ |PRCC| with t-test
             stringsAsFactors = FALSE)
}))
# significance: t = prcc*sqrt((n-2-(k_used-1))/(1-prcc^2)); k_used is per-outcome
# because the placebo outcomes drop the drug parameters.
prcc_tbl$t  <- with(prcc_tbl, prcc * sqrt(pmax(n - 2 - (k_used - 1), 1) / (1 - prcc^2)))
prcc_tbl$p  <- 2 * pt(-abs(prcc_tbl$t), df = pmax(prcc_tbl$n - 2 - (prcc_tbl$k_used - 1), 1))
# Multiplicity: Benjamini-Hochberg applied WITHIN each outcome (not pooled
# across outcomes), since the parameter set tested differs by outcome and each
# outcome is a separate scientific question. Raw p is retained alongside p_adj.
prcc_tbl$p_adj <- ave(prcc_tbl$p, prcc_tbl$outcome,
                      FUN = function(pp) p.adjust(pp, method = "BH"))
prcc_tbl$sig <- ifelse(prcc_tbl$p_adj < 0.05, "*", "")
write.csv(prcc_tbl[, c("parameter","outcome","method","prcc","p","p_adj","sig","n","k_used")],
          "outputs/gsa_prcc.csv", row.names = FALSE)

# --- Tornado figure ---------------------------------------------------------
plotd <- prcc_tbl
plotd$outcome <- factor(plotd$outcome, levels = unname(outcomes))
# order parameters within each facet by |prcc|
plotd <- plotd %>% group_by(outcome) %>%
  mutate(parameter = factor(parameter, levels = parameter[order(abs(prcc))])) %>%
  ungroup()

p <- ggplot(plotd, aes(x = prcc, y = parameter, fill = prcc > 0)) +
  geom_col(width = 0.7, alpha = 0.9) +
  geom_vline(xintercept = 0, colour = "grey40") +
  # ncol = 2: there are now 4 outcomes (the binary establishment outcome added)
  facet_wrap(~ outcome, scales = "free", ncol = 2) +
  scale_fill_manual(values = c(`TRUE` = "#d73027", `FALSE` = "#4575b4"),
                    labels = c(`TRUE` = "positive", `FALSE` = "negative"),
                    name = "PRCC sign") +
  labs(x = "Partial Rank Correlation Coefficient (PRCC)", y = NULL,
       title = "Global Sensitivity Analysis (LHS + PRCC)",
       subtitle = sprintf("N = %d Latin Hypercube samples; log-uniform parameter ranges", N)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave("outputs/gsa_prcc_tornado.png", p, width = 13, height = 7, dpi = 200)

# --- Text summary -----------------------------------------------------------
sink("outputs/gsa_summary.txt")
cat("GLOBAL SENSITIVITY ANALYSIS (LHS + PRCC)\n")
cat(sprintf("Date: %s | N=%d | parameters=%d\n\n", Sys.Date(), N, k))

# --- Conditioning: which designs enter which statistic ----------------------
est <- samples$established == 1
cat("DESIGN CONDITIONING\n")
cat(sprintf("  n_total (LHS designs):        %d\n", N))
cat(sprintf("  n_valid (placebo sim solved): %d\n", n_valid))
cat(sprintf("  n_established (placebo mortality > %.0e): %d (%.1f%%)\n",
    ESTABLISH_MORT_THRESHOLD, n_established, 100 * n_established / N))
cat(sprintf("  excluded (never established): %d (%.1f%%)\n",
    N - n_established, 100 * (N - n_established) / N))
cat("  NOTE: RRR is undefined when the infection never establishes, so ALL RRR\n")
cat("  statistics below are CONDITIONAL on establishment (n_established).\n\n")

cat(sprintf("Day-1 combination RRR [CONDITIONAL on establishment, n=%d]: median=%.1f%%  IQR=[%.1f, %.1f]%%  (5th pct=%.1f%%)\n",
    sum(is.finite(samples$combo_RRR)),
    100*median(samples$combo_RRR, na.rm=TRUE),
    100*quantile(samples$combo_RRR, .25, na.rm=TRUE),
    100*quantile(samples$combo_RRR, .75, na.rm=TRUE),
    100*quantile(samples$combo_RRR, .05, na.rm=TRUE)))
cat(sprintf("Placebo mortality [ALL designs, n=%d]: median=%.1f%%  IQR=[%.1f, %.1f]%%\n",
    n_valid,
    100*median(samples$placebo_mortality, na.rm=TRUE),
    100*quantile(samples$placebo_mortality, .25, na.rm=TRUE),
    100*quantile(samples$placebo_mortality, .75, na.rm=TRUE)))
cat(sprintf("Placebo mortality [ESTABLISHING subset only, n=%d]: median=%.1f%%  IQR=[%.1f, %.1f]%%\n\n",
    n_established,
    100*median(samples$placebo_mortality[est], na.rm=TRUE),
    100*quantile(samples$placebo_mortality[est], .25, na.rm=TRUE),
    100*quantile(samples$placebo_mortality[est], .75, na.rm=TRUE)))

cat("Significance is based on Benjamini-Hochberg adjusted p-values (p_adj < 0.05),\n")
cat("computed within each outcome. Drug parameters are excluded from the\n")
cat("placebo-arm outcomes because they cannot act in the placebo arm.\n\n")

for (oc in names(outcomes)) {
  sub <- prcc_tbl[prcc_tbl$outcome == outcomes[[oc]], ]
  sub <- sub[order(-abs(sub$prcc)), ]
  label <- if (oc == "established")
    "Top partial rank point-biserial drivers (BINARY outcome, NOT a PRCC on a continuous outcome)"
  else "Top PRCC drivers"
  stat_lab <- if (oc == "established") "rho_pb" else "PRCC"
  cat(sprintf("--- %s: %s ---\n", label, outcomes[[oc]]))
  # Print only rows that survive the multiplicity correction (capped at 8).
  # Printing a fixed top-8 regardless of magnitude previously listed
  # noise-level coefficients (|PRCC| ~ 1/sqrt(N) = 0.032) as "drivers".
  sub <- sub[sub$sig == "*", , drop = FALSE]
  if (nrow(sub) == 0) {
    cat("  no parameters significant after BH correction\n")
  } else {
    for (i in seq_len(min(8, nrow(sub))))
      cat(sprintf("  %-14s %s=%+.3f  p_adj=%.3g %s\n",
                  sub$parameter[i], stat_lab, sub$prcc[i], sub$p_adj[i], sub$sig[i]))
  }
  cat("\n")
}
sink()
cat("Saved: gsa_lhs_samples.csv, gsa_prcc.csv, gsa_prcc_tornado.png, gsa_summary.txt\n")
