# Plot cytokine trajectories for the placebo arm
# Uses pure R RK4 since deSolve lsoda crashes on this machine

# Pure R RK4 solver (preserves state names)
rk4_pure <- function(y, times, func, parms) {
  n <- length(times)
  ny <- length(y)
  state_names <- names(y)
  out <- matrix(NA, nrow = n, ncol = ny + 1)
  colnames(out) <- c("time", state_names)
  out[1, ] <- c(times[1], y)
  for (i in 2:n) {
    t <- times[i - 1]
    h <- times[i] - times[i - 1]
    y_curr <- setNames(out[i - 1, -1], state_names)
    k1 <- func(t, y_curr, parms)[[1]]
    k2 <- func(t + h / 2, y_curr + h / 2 * k1, parms)[[1]]
    k3 <- func(t + h / 2, y_curr + h / 2 * k2, parms)[[1]]
    k4 <- func(t + h, y_curr + h * k3, parms)[[1]]
    y_next <- y_curr + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    names(y_next) <- state_names
    out[i, ] <- c(times[i], y_next)
  }
  out
}

library(deSolve)
source("model/hantavirus_qsp.R")
source("model/parameters.R")
source("R/pk_models.R")
source("R/virtual_population.R")
source("R/simulate_trial.R")
source("R/pd_models.R")

pop <- read.csv("outputs/virtual_population.csv")
pars <- get_parameters()

set.seed(42)
n_pts <- 20
sample_ids <- sample(pop$patient_id, n_pts)
dt <- 0.25
times <- seq(0, 21, by = dt)

all_C_pro  <- matrix(NA, nrow = n_pts, ncol = length(times))
all_C_anti <- matrix(NA, nrow = n_pts, ncol = length(times))
all_P      <- matrix(NA, nrow = n_pts, ncol = length(times))

cat("Running", n_pts, "simulations with pure R RK4 (dt =", dt, ")...\n")
for (i in seq_along(sample_ids)) {
  pid <- sample_ids[i]
  row <- pop[pop$patient_id == pid, ]
  pars_i <- apply_vpop_to_model(row, pars)
  y0 <- build_symptom_onset_state(pars_i)
  sim <- rk4_pure(y0, times, hantavirus_qsp_ode, pars_i)
  all_C_pro[i, ]  <- sim[, "C_pro"]
  all_C_anti[i, ] <- sim[, "C_anti"]
  all_P[i, ]      <- sim[, "P"]
  if (i %% 5 == 0) cat("  ", i, "/", n_pts, "\n")
}

# Median and IQR
med_C_pro  <- apply(all_C_pro, 2, median, na.rm = TRUE)
q25_C_pro  <- apply(all_C_pro, 2, quantile, 0.25, na.rm = TRUE)
q75_C_pro  <- apply(all_C_pro, 2, quantile, 0.75, na.rm = TRUE)
med_C_anti <- apply(all_C_anti, 2, median, na.rm = TRUE)
q25_C_anti <- apply(all_C_anti, 2, quantile, 0.25, na.rm = TRUE)
q75_C_anti <- apply(all_C_anti, 2, quantile, 0.75, na.rm = TRUE)
med_P      <- apply(all_P, 2, median, na.rm = TRUE)
q25_P      <- apply(all_P, 2, quantile, 0.25, na.rm = TRUE)
q75_P      <- apply(all_P, 2, quantile, 0.75, na.rm = TRUE)

# PNG output
png("outputs/cytokine_trajectories_placebo.png", width = 1400, height = 1000, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), mgp = c(2.5, 0.8, 0))

# Panel A: C_pro
ylimA <- c(0, max(q75_C_pro) * 1.1)
plot(times, med_C_pro, type = "n", xlab = "Days post symptom onset",
     ylab = "C_pro (AU)", main = "Pro-inflammatory (C_pro)",
     ylim = ylimA, cex.main = 1.1)
polygon(c(times, rev(times)), c(q25_C_pro, rev(q75_C_pro)),
        col = rgb(0.9, 0.2, 0.2, 0.25), border = NA)
lines(times, med_C_pro, col = "firebrick", lwd = 2.5)
peak_day <- times[which.max(med_C_pro)]
abline(v = peak_day, lty = 2, col = "firebrick4")
text(peak_day + 1, ylimA[2] * 0.95, paste0("Peak day ", round(peak_day, 1)),
     col = "firebrick4", cex = 0.85)
legend("topright", c("Median", "IQR"),
       col = c("firebrick", rgb(0.9, 0.2, 0.2, 0.5)),
       lwd = c(2.5, NA), pch = c(NA, 15), pt.cex = 2, bty = "n", cex = 0.8)

# Panel B: C_anti
ylimB <- c(0, max(q75_C_anti) * 1.1)
plot(times, med_C_anti, type = "n", xlab = "Days post symptom onset",
     ylab = "C_anti (AU)", main = "Anti-inflammatory (C_anti)",
     ylim = ylimB, cex.main = 1.1)
polygon(c(times, rev(times)), c(q25_C_anti, rev(q75_C_anti)),
        col = rgb(0.2, 0.4, 0.9, 0.25), border = NA)
lines(times, med_C_anti, col = "steelblue", lwd = 2.5)
peak_day <- times[which.max(med_C_anti)]
abline(v = peak_day, lty = 2, col = "steelblue4")
text(peak_day + 1, ylimB[2] * 0.95, paste0("Peak day ", round(peak_day, 1)),
     col = "steelblue4", cex = 0.85)
legend("topright", c("Median", "IQR"),
       col = c("steelblue", rgb(0.2, 0.4, 0.9, 0.5)),
       lwd = c(2.5, NA), pch = c(NA, 15), pt.cex = 2, bty = "n", cex = 0.8)

# Panel C: Both normalized
plot(times, med_C_pro / max(med_C_pro), type = "n",
     xlab = "Days post symptom onset", ylab = "Normalized (0-1)",
     main = "Normalized C_pro vs C_anti", cex.main = 1.1)
polygon(c(times, rev(times)),
        c(q25_C_pro / max(med_C_pro), rev(q75_C_pro / max(med_C_pro))),
        col = rgb(0.9, 0.2, 0.2, 0.2), border = NA)
lines(times, med_C_pro / max(med_C_pro), col = "firebrick", lwd = 2.5)
polygon(c(times, rev(times)),
        c(q25_C_anti / max(med_C_anti), rev(q75_C_anti / max(med_C_anti))),
        col = rgb(0.2, 0.4, 0.9, 0.2), border = NA)
lines(times, med_C_anti / max(med_C_anti), col = "steelblue", lwd = 2.5)
legend("topright", c("C_pro (pro-inflam)", "C_anti (anti-inflam)"),
       col = c("firebrick", "steelblue"), lwd = 2.5, bty = "n", cex = 0.85)

# Panel D: Permeability
ylimD <- c(0, max(q75_P) * 1.1)
plot(times, med_P, type = "n", xlab = "Days post symptom onset",
     ylab = "P (AU)", main = "Vascular permeability (P)",
     ylim = ylimD, cex.main = 1.1)
polygon(c(times, rev(times)), c(q25_P, rev(q75_P)),
        col = rgb(0.5, 0.3, 0.7, 0.25), border = NA)
lines(times, med_P, col = "darkviolet", lwd = 2.5)
peak_day <- times[which.max(med_P)]
abline(v = peak_day, lty = 2, col = "darkviolet")
text(peak_day + 1, ylimD[2] * 0.95, paste0("Peak day ", round(peak_day, 1)),
     col = "darkviolet", cex = 0.85)

dev.off()

cat("\nPlot saved to outputs/cytokine_trajectories_placebo.png\n")
cat("\nSummary (median of", n_pts, "patients):\n")
cat("C_pro peak:  day", round(times[which.max(med_C_pro)], 1),
    "  value =", round(max(med_C_pro), 4), "\n")
cat("C_anti peak: day", round(times[which.max(med_C_anti)], 1),
    "  value =", round(max(med_C_anti), 4), "\n")
cat("P peak:      day", round(times[which.max(med_P)], 1),
    "  value =", round(max(med_P), 2), "\n")
cat("Offset (C_anti \u2013 C_pro):",
    round(times[which.max(med_C_anti)] - times[which.max(med_C_pro)], 1), "days\n")
