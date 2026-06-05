library(deSolve)

source("model/hantavirus_qsp.R")
source("model/parameters.R")
source("R/pk_models.R")
source("R/pd_models.R")
source("R/virtual_population.R")
source("R/simulate_trial.R")

# Load a representative virtual patient
vpop <- read.csv("outputs/virtual_population.csv", stringsAsFactors = FALSE)
set.seed(42)
pt <- vpop[sample(nrow(vpop), 1), ]

pars <- get_parameters()
pars <- apply_vpop_to_model(pt, pars)

# Test grid: delta_CD8 x k_neut_IgG
delta_CD8_vals <- c(1e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2)
k_neut_IgG_vals <- c(0.005, 0.02, 0.05, 0.1, 0.2)

cat("=== ADAPTIVE IMMUNITY PARAMETER CALIBRATION ===\n")
cat(sprintf("Patient: %d (adaptive_strength=%.2f)\n", pt$patient_id, pt$adaptive_strength))
cat("Testing delta_CD8 x k_neut_IgG grid\n\n")

cat(sprintf("%-12s %-12s | %-10s %-10s %-10s %-10s %-8s %-6s\n",
            "delta_CD8", "k_neut_IgG", "V_peak", "V_d13", "V_d21", "Rebound?", "CD8max", "IgGmax"))
cat(paste(rep("-", 100), collapse = ""), "\n")

for (dCD8 in delta_CD8_vals) {
  for (kIgG in k_neut_IgG_vals) {
    pars_i <- pars
    pars_i$delta_CD8 <- dCD8
    pars_i$k_neut_IgG <- kIgG

    suppressWarnings({
      out <- try(simulate_patient(pars_i, "ribavirin", t_start = 0, t_end = 21, dt = 0.5), silent = TRUE)
    })

    if (inherits(out, "try-error") || is.null(out) || is.null(out$V)) {
      cat(sprintf("%-12.1e %-12.1e | %-10s\n", dCD8, kIgG, "ERROR"))
      next
    }

    V <- out$V
    times <- out$time

    V_peak <- max(V, na.rm = TRUE)
    V_d13 <- V[which.min(abs(times - 13))]
    V_d21 <- tail(V, 1)

    rebound <- V_d21 > V_d13 * 1.1

    # Extract CD8 and IgG max
    if (!is.null(out$y) && "CD8" %in% colnames(out$y)) {
      CD8_max <- max(out$y[, "CD8"], na.rm = TRUE)
      IgG_max <- max(out$y[, "IgG"], na.rm = TRUE)
    } else {
      CD8_max <- NA
      IgG_max <- NA
    }

    cat(sprintf("%-12.1e %-12.1e | %-10.2e %-10.2e %-10.2e %-10s %-8.1f %-6.1f\n",
                dCD8, kIgG, V_peak, V_d13, V_d21, ifelse(rebound, "YES", "no"), CD8_max, IgG_max))
  }
}

# Now check placebo arm with selected "good" parameter combos
cat("\n\n=== PLACEBO ARM CHECK (no drug, adaptive only) ===\n\n")

for (dCD8 in c(1e-3, 5e-3, 1e-2)) {
  for (kIgG in c(0.02, 0.05, 0.1)) {
    pars_i <- pars
    pars_i$delta_CD8 <- dCD8
    pars_i$k_neut_IgG <- kIgG

    suppressWarnings({
      out <- try(simulate_patient(pars_i, "placebo", t_start = 0, t_end = 21, dt = 0.5), silent = TRUE)
    })
    if (inherits(out, "try-error") || is.null(out) || is.null(out$V)) next

    V <- out$V
    times <- out$time
    V_peak <- max(V, na.rm = TRUE)
    V_d21 <- tail(V, 1)

    if (!is.null(out$y) && "K" %in% colnames(out$y)) {
      K <- out$y[, "K"]
      L <- out$y[, "L"]
      C_pro <- out$y[, "C_pro"]
      K_max <- max(K, na.rm = TRUE)
      L_max <- max(L, na.rm = TRUE)
      C_max <- max(C_pro, na.rm = TRUE)

      m_risk <- pars_i$w_K_HFRS * K_max / (K_max + pars_i$K50) +
                pars_i$w_L_HFRS * L_max / (L_max + pars_i$L50) +
                pars_i$w_C_HFRS * C_max / (C_max + pars_i$C50)
      mort <- ifelse(m_risk > 0.35, 1, 0)
    } else {
      mort <- NA
    }

    cat(sprintf("delta_CD8=%.0e, k_neut_IgG=%.0e | V_peak=%.2e, V_day21=%.2e, mortality=%d\n",
                dCD8, kIgG, V_peak, V_d21, mort))
  }
}
