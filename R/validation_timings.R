# Extract representative-patient peak/nadir timings for the external-validation table.
# t = 0 is symptom onset (= day 5 post-infection). Times reported post-symptom-onset.
source("model/parameters.R"); source("model/hantavirus_qsp.R")
source("R/pk_models.R"); source("R/pd_models.R"); source("R/simulate_trial.R")
p <- get_parameters()
s <- simulate_patient(p, "placebo", t_start = 0, t_end = 21, dt = 0.1)
onset <- p$symptom_onset_day
tp <- function(col, fn) s$time[fn(s[[col]])]
cat(sprintf("V_peak:   t_post_onset=%.1f d (t_post_infection=%.1f) value=%.3e copies/mL\n",
            tp("V", which.max), tp("V", which.max) + onset, max(s$V)))
cat(sprintf("C_pro:    t_post_onset=%.1f d (value=%.1f AU)\n", tp("C_pro", which.max), max(s$C_pro)))
cat(sprintf("PLT nadir:t_post_onset=%.1f d value=%.0f /uL\n", tp("PLT", which.min), min(s$PLT)))
cat(sprintf("K_peak:   t_post_onset=%.1f d value=%.1f\n", tp("K", which.max), max(s$K)))
cat(sprintf("L_peak:   t_post_onset=%.1f d value=%.1f\n", tp("L", which.max), max(s$L)))
cat(sprintf("CD8_E:    t_post_onset=%.1f d (t_post_infection=%.1f)\n", tp("CD8_E", which.max), tp("CD8_E", which.max) + onset))
cat(sprintf("IgM:      t_post_onset=%.1f d\n", tp("IgM", which.max)))
cat(sprintf("IgG:      t_post_onset=%.1f d (end-of-sim value=%.3f)\n", tp("IgG", which.max), tail(s$IgG, 1)))
