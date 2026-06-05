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
write.csv(samples, "outputs/gsa_lhs_samples.csv", row.names = FALSE)
cat(sprintf("Completed %d samples (%d valid placebo sims)\n",
            N, sum(is.finite(res$placebo_mortality))))

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

outcomes <- c(placebo_mortality = "Placebo mortality",
              combo_RRR        = "Day-1 combination RRR",
              placebo_Vpeak    = "Peak viral load (placebo)")
Xmat <- as.matrix(design)
prcc_tbl <- do.call(rbind, lapply(names(outcomes), function(oc) {
  pr <- prcc_vec(Xmat, samples[[oc]])
  data.frame(parameter = pnames, outcome = outcomes[[oc]],
             prcc = pr$prcc, n = pr$n,
             # two-sided significance threshold ~ |PRCC| with t-test
             stringsAsFactors = FALSE)
}))
# significance: t = prcc*sqrt((n-2-(k-1))/(1-prcc^2))
prcc_tbl$t  <- with(prcc_tbl, prcc * sqrt(pmax(n - 2 - (k - 1), 1) / (1 - prcc^2)))
prcc_tbl$p  <- 2 * pt(-abs(prcc_tbl$t), df = pmax(prcc_tbl$n - 2 - (k - 1), 1))
prcc_tbl$sig <- ifelse(prcc_tbl$p < 0.05, "*", "")
write.csv(prcc_tbl[, c("parameter","outcome","prcc","p","sig","n")],
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
  facet_wrap(~ outcome, scales = "free", ncol = 3) +
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
cat(sprintf("Day-1 combination RRR: median=%.1f%%  IQR=[%.1f, %.1f]%%  (5th pct=%.1f%%)\n",
    100*median(samples$combo_RRR, na.rm=TRUE),
    100*quantile(samples$combo_RRR, .25, na.rm=TRUE),
    100*quantile(samples$combo_RRR, .75, na.rm=TRUE),
    100*quantile(samples$combo_RRR, .05, na.rm=TRUE)))
cat(sprintf("Placebo mortality: median=%.1f%%  IQR=[%.1f, %.1f]%%\n\n",
    100*median(samples$placebo_mortality, na.rm=TRUE),
    100*quantile(samples$placebo_mortality, .25, na.rm=TRUE),
    100*quantile(samples$placebo_mortality, .75, na.rm=TRUE)))
for (oc in names(outcomes)) {
  sub <- prcc_tbl[prcc_tbl$outcome == outcomes[[oc]], ]
  sub <- sub[order(-abs(sub$prcc)), ]
  cat(sprintf("--- Top PRCC drivers: %s ---\n", outcomes[[oc]]))
  for (i in 1:min(8, nrow(sub)))
    cat(sprintf("  %-14s PRCC=%+.3f %s\n", sub$parameter[i], sub$prcc[i], sub$sig[i]))
  cat("\n")
}
sink()
cat("Saved: gsa_lhs_samples.csv, gsa_prcc.csv, gsa_prcc_tornado.png, gsa_summary.txt\n")
