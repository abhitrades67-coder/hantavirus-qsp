# -*- coding: utf-8 -*-
"""Re-derive every number quoted in the results and check it against the text.

The tables are generated from outputs/ by build_tables.py and cannot drift.
The numbers quoted in running prose can, so each one is recomputed here from
the same CSV files and compared with the sentence that carries it. A mismatch
fails the build.

Values are compared as they are written, that is after rounding to the number
of digits the manuscript uses, so a claim is accepted only if it reads the way
the data say it should.

Usage, from the repository root:
    python "manuscript/Royal/_src/check_numbers.py"
"""
import csv
import io
import os
import re
import statistics as st
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(PKG))
OUT = os.path.join(REPO, "outputs")

fails = []


def rd(name):
    with io.open(os.path.join(OUT, name), encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


MAN = io.open(os.path.join(HERE, "manuscript.md"), encoding="utf-8").read()
BODY = MAN.split("## References")[0]


def says(fragment, label):
    """The manuscript must contain this exact fragment."""
    if fragment in BODY:
        print("  ok    %-46s %s" % (label, fragment[:76]))
    else:
        fails.append(label)
        print("  FAIL  %-46s not in the text: %s" % (label, fragment[:76]))


def num(x):
    return float(x)


# ---------------------------------------------------------------- 3.1 placebo
print("\n== 3.1 representative placebo patient ==")
D = {r["quantity"]: num(r["value"]) for r in rd("representative_placebo.csv")}
says("peaked at %.2f × 10^7^ copies ml^−1^ on day %d"
     % (D["V_peak"] / 1e7, round(D["t_V_peak"])), "peak viral load and day")
says("permeability at %d on day %.1f" % (round(D["P_peak"]), D["t_P_peak"]),
     "permeability peak")
says("renal injury at %.1f on day %.1f" % (D["K_peak"], D["t_K_peak"]),
     "renal injury peak")
says("platelet nadir of {:,.0f}".format(D["PLT_nadir"]), "platelet nadir")
says("pulmonary injury peaking at %.1f" % D["L_peak"], "pulmonary injury peak")
says("Mortality was %.3f per cent" % D["mortality_pct"], "mortality")
says("dialysis risk %.1f per cent" % D["dialysis_pct"], "dialysis risk")
says("support risk %.1f per cent" % D["ecmo_pct"], "extracorporeal support risk")
says("infected cells reach %.2f × 10^6^" % (D["I_peak"] / 1e6),
     "infected-cell peak")
m = re.search(r"within ([\d.]+) per cent of \*pI\*/\*c\*", BODY)
if m and num(m.group(1)) >= D["V_peak_vs_pI_over_c_pct"]:
    print("  ok    %-46s stated bound %s%% holds (actual %.4f%%)"
          % ("quasi-steady-state bound", m.group(1),
             D["V_peak_vs_pI_over_c_pct"]))
else:
    fails.append("quasi-steady-state bound")
    print("  FAIL  quasi-steady-state bound: actual %.4f%%"
          % D["V_peak_vs_pI_over_c_pct"])

ce = {r["arm"]: r for r in rd("clinical_endpoints.csv")}
says("placebo mortality averaged %.2f per cent"
     % (100 * num(ce["placebo"]["mortality_risk_mean"])), "population mean")
says("median of %.2f per cent"
     % (100 * num(ce["placebo"]["mortality_median"])), "population median")

# ------------------------------------------------------------- 3.2 mechanisms
print("\n== 3.2 mechanism weights ==")
W = {(r["quantity"], r["component"]): num(r["percent"])
     for r in rd("mechanism_weights.csv")}


def w(q, c):
    return W[(q, c)]


says("Renal injury contributes %.1f per cent"
     % w("Mortality probability", "Renal injury (peak and integral)"), "renal share")
says("pulmonary injury %.1f per cent"
     % w("Mortality probability", "Pulmonary injury"), "pulmonary share")
says("cytokine term contributes %.3f per cent"
     % w("Mortality probability", "Pro-inflammatory cytokines"), "cytokine share")
says("accounts for %.4f per cent of infected-cell removal"
     % w("Infected-cell clearance", "NK-mediated"), "NK share")
says("natural turnover (%.1f per cent)"
     % w("Infected-cell clearance", "Natural turnover"), "turnover share")
says("killing (%.1f per cent)"
     % w("Infected-cell clearance", "CD8+ killing"), "CD8 share")
ab = w("Viral clearance", "IgM neutralisation") + w("Viral clearance", "IgG neutralisation")
says("accounts for %.2f per cent of viral clearance" % ab, "antibody share")
says("at most %.5f per cent"
     % w("Type I interferon", "Maximum suppression of viral production"),
     "interferon suppression")
says("runs at %.5f per cent of its capacity"
     % w("Cytokine auto-amplification", "Fraction of capacity used"),
     "cytokine auto-amplification")

# --------------------------------------------------------- 3.3 antiviral effect
print("\n== 3.3 antiviral effect and treatment window ==")
win = {(r["arm"], int(num(r["treatment_day"]))): r for r in rd("treatment_window.csv")}
pl_auc = num(win[("ribavirin", 1)]["placebo_V_AUC"])
pl_mort = 100 * num(win[("ribavirin", 1)]["placebo_mortality"])


def auc_pct(arm):
    return 100 * num(win[(arm, 1)]["V_AUC"]) / pl_auc


says("integral to %.1f per cent of placebo" % auc_pct("ribavirin"), "ribavirin exposure")
says("favipiravir to %.1f per cent" % auc_pct("favipiravir"), "favipiravir exposure")
says("combination reached %.1f per cent" % auc_pct("combination"), "combination exposure")


def mort(arm, day):
    return 100 * num(win[(arm, day)]["mortality_prob"])


def rrr(arm, day):
    return num(win[(arm, day)]["RRR"])


says("combination mortality was %.2f per cent, a %.1f per cent relative reduction"
     % (mort("combination", 1), rrr("combination", 1)), "combination day 1")
says("falling to %.2f per cent (%.1f per cent) by day 4"
     % (mort("combination", 4), rrr("combination", 4)), "combination day 4")
says("%.2f per cent (%.1f per cent) by day 7"
     % (mort("combination", 7), rrr("combination", 7)), "combination day 7")
says("Ribavirin gave %.2f per cent on day 1 (%.1f per cent)"
     % (mort("ribavirin", 1), rrr("ribavirin", 1)), "ribavirin day 1")
says("falling to %.2f per cent on day 7 (%.1f per cent)"
     % (mort("ribavirin", 7), rrr("ribavirin", 7)), "ribavirin day 7")
says("gave %.2f per cent on day 1 (%.1f per cent) and %.2f per cent on day 2"
     % (mort("favipiravir", 1), rrr("favipiravir", 1), mort("favipiravir", 2)),
     "favipiravir days 1 and 2")

org = rd("organ_peak_data.csv")
plc = rd("organ_peak_placebo.csv")
pm = {k: st.median(num(r[k]) for r in plc) for k in plc[0]}
KEYS = ["P_peak", "K_peak", "L_peak"]


def trio(arm):
    r = [x for x in org if x["arm"] == arm and int(num(x["treatment_day"])) == 1][0]
    return ", ".join("%.1f" % (100 * num(r[k]) / pm[k]) for k in KEYS)


says("injury to %s per cent of placebo" % trio("ribavirin"), "ribavirin organ")
says("combination to %s per cent" % trio("combination"), "combination organ")
says("favipiravir reached %s per cent" % trio("favipiravir"), "favipiravir organ")
for day, word in ((5, "By day 5"), (7, "by day 7")):
    lo = min(100 * num(r[k]) / pm[k] for r in org
             if int(num(r["treatment_day"])) == day for k in KEYS)
    stated = re.search(r"%s .{0,60}?at or above (\d+) per cent" % word, BODY)
    if stated and int(stated.group(1)) <= lo:
        print("  ok    %-46s day %d lower bound %s%% holds (actual %.2f%%)"
              % ("organ injury floor", day, stated.group(1), lo))
    else:
        fails.append("organ injury floor day %d" % day)
        print("  FAIL  organ injury floor day %d: actual minimum %.4f%%" % (day, lo))

# ---------------------------------------------------------------- 3.4 duration
print("\n== 3.4 duration rather than timing ==")
ad = rd("adaptive_peak_data.csv")
adp = rd("adaptive_peak_placebo.csv")
am = {k: st.median(num(r[k]) for r in adp) for k in adp[0]}


def igm(arm, day):
    return num([r for r in ad if r["arm"] == arm
                and int(num(r["treatment_day"])) == day][0]["IgM_peak"])


says("peak IgM was %.3f with a day-1 start against %.3f for day 2 and %.3f for placebo"
     % (igm("favipiravir", 1), igm("favipiravir", 2), am["IgM_peak"]), "IgM peaks")
says("is a %.0f per cent increase"
     % (100 * (igm("favipiravir", 1) / am["IgM_peak"] - 1)), "IgM increase")

dur = {(r["arm"], int(num(r["start_day"])), int(num(r["duration_days"]))): r
       for r in rd("duration_start_data.csv")}


def dm(arm, sd, d):
    return 100 * num(dur[(arm, sd, d)]["mortality_prob"])


says("%.1f per cent of the initial target cells remain, against %.1f per cent"
     % (100 * num(dur[("favipiravir", 1, 15)]["T_end_frac"]),
        100 * num(dur[("favipiravir", 2, 15)]["T_end_frac"])), "target-cell sparing")

reb = {(r["arm"], r["start_day"]): r for r in rd("rebound_diagnostics.csv")}
r1 = reb[("favipiravir", "1")]
r2 = reb[("favipiravir", "2")]
rp = reb[("placebo", "NA")]
says("viral load holds near %.1f × 10^5^" % (num(r1["V_on_treatment_min"]) / 1e5),
     "on-treatment viral load")
says("infected cells near %.1f × 10^5^" % (num(r1["I_on_treatment_min"]) / 1e5),
     "on-treatment infected cells")
says("second peak of %.2f × 10^7^ copies ml^−1^ on day %.1f"
     % (num(r1["V_rebound_peak"]) / 1e7, num(r1["t_V_rebound"])), "rebound peak")
says("peaks at %d on day %.1f" % (round(num(r1["P_rebound_peak"])),
                                  num(r1["t_P_rebound"])), "rebound permeability")
says("renal injury at %.1f on day %.1f" % (num(r1["K_rebound_peak"]),
                                           num(r1["t_K_rebound"])), "rebound renal")
says("rebounds to %.2f × 10^7^ and a renal peak of %.1f"
     % (num(r2["V_rebound_peak"]) / 1e7, num(r2["K_rebound_peak"])), "day-2 rebound")
short = 100 * (1 - num(r1["V_rebound_peak"]) / num(rp["V_rebound_peak"]))
lag = num(r1["t_V_rebound"]) - num(rp["t_V_rebound"])
m = re.search(r"within (\d+) per cent of the untreated peak and (\d+) days later", BODY)
if m and num(m.group(1)) >= short and abs(num(m.group(2)) - lag) < 0.5:
    print("  ok    %-46s %s%% and %s days (actual %.1f%%, %.1f days)"
          % ("rebound comparison", m.group(1), m.group(2), short, lag))
else:
    fails.append("rebound comparison")
    print("  FAIL  rebound comparison: actual %.2f%% short, %.2f days later"
          % (short, lag))

d1_21, d2_21 = mort("favipiravir", 1), mort("favipiravir", 2)
says("is %.2f percentage points at 21 days but %.2f points at 42 days (%.2f against %.2f per cent)"
     % (d1_21 - d2_21, dm("favipiravir", 1, 15) - dm("favipiravir", 2, 15),
        dm("favipiravir", 1, 15), dm("favipiravir", 2, 15)), "21 versus 42 days")
says("falls from %.2f per cent on a 10-day course to %.2f per cent on a 35-day course"
     % (dm("favipiravir", 1, 10), dm("favipiravir", 1, 35)), "course extension")
says("renal peak returns to day %.1f" % num(dur[("favipiravir", 1, 35)]["t_K_peak"]),
     "renal peak after extension")
says("moving from %.1f to %.1f to %.1f for 15, 21 and 28-day courses"
     % tuple(num(dur[("favipiravir", 1, d)]["t_K_peak"]) for d in (15, 21, 28)),
     "renal peak day by duration")
says("%.2f against %.2f per cent for favipiravir and %.2f against %.2f per cent for the combination"
     % (dm("favipiravir", 1, 35), dm("favipiravir", 2, 35),
        dm("combination", 1, 35), dm("combination", 2, 35)), "35-day ordering")

hum = {num(r["antibody_scaling"]): r for r in rd("humoral_scan.csv")}
says("antibodies perform %.0f per cent of viral clearance"
     % num(hum[100]["antibody_share_day1"]), "antibody share at 100x")

# ---------------------------------------------------------------- 3.5 ablation
print("\n== 3.5 ribavirin non-antiviral terms ==")
abl = {r["scenario"]: r for r in rd("ablation_summary.csv")}


def am_(k):
    return 100 * num(abl[k]["mortality_prob"])


says("mortality from %.3f to %.3f per cent"
     % (am_("ribavirin_day1_default"), am_("ribavirin_day1_no_immuno")),
     "immunomodulatory ablation")
says("renal injury unchanged at %.2f" % num(abl["ribavirin_day1_default"]["K_peak"]),
     "renal injury unchanged")
says("raised mortality to %.3f per cent and renal injury to %.2f"
     % (am_("ribavirin_day1_no_endothelial"),
        num(abl["ribavirin_day1_no_endothelial"]["K_peak"])), "endothelial ablation")
says("together (%.3f per cent)" % am_("ribavirin_day1_no_nonantiviral"), "all three")
says("(%.3f per cent by default, %.3f per cent without the endothelial term)"
     % (am_("combination_day1_default"), am_("combination_day1_no_endothelial")),
     "combination ablation")

# ---------------------------------------------------------------- 3.6 adaptive
print("\n== 3.6 adaptive immunity ==")
dev = max(abs(100 * num(r[k]) / am[k] - 100) for r in ad
          for k in ("CD8_E_peak", "CD4_peak"))
m = re.search(r"within ([\d.]+) per cent of placebo across all arms", BODY)
if m and num(m.group(1)) >= dev:
    print("  ok    %-46s stated %s%% holds (actual %.4f%%)"
          % ("CD8 and CD4 invariance", m.group(1), dev))
else:
    fails.append("CD8 and CD4 invariance")
    print("  FAIL  CD8 and CD4 invariance: actual %.4f%%" % dev)


def apct(arm, day, key):
    r = [x for x in ad if x["arm"] == arm and int(num(x["treatment_day"])) == day][0]
    return 100 * num(r[key]) / am[key]


says("reaching %.0f per cent of placebo with day-1 favipiravir and %.0f per cent"
     % (apct("favipiravir", 1, "IgM_peak"), apct("combination", 1, "IgM_peak")),
     "IgM modulation")
igg = [100 * num(r["IgG_peak"]) / am["IgG_peak"] for r in ad]
says("modest (%.0f to %.0f per cent)" % (min(igg), max(igg)), "IgG modulation")

# ------------------------------------------------------------ 3.7 prophylaxis
print("\n== 3.7 post-exposure prophylaxis ==")
pep = rd("preexposure_results.csv")
says("All %d prophylaxis scenarios" % len(pep), "scenario count")
early = [num(r["V_peak"]) for r in pep
         if r["arm"] != "placebo" and int(num(r["pep_day"])) <= 1]
says("at %.1f to %.1f × 10^3^ copies ml^−1^"
     % (min(early) / 1e3, max(early) / 1e3), "day 0 to 1 peaks")


def pv(arm, day, key="V_peak"):
    return num([r for r in pep if r["arm"] == arm
                and int(num(r["pep_day"])) == day][0][key])


says("against %.2f × 10^7^ for placebo" % (pv("placebo", 0) / 1e7), "placebo peak")
d2 = min(pv(a, 2) for a in ("favipiravir", "combination"))
says("reached only %.2f × 10^4^" % (d2 / 1e4), "day 2 lowest peak")
import math
for arm, sci, tail in (("combination", 1e6, "and mortality to %.3f per cent"),
                       ("favipiravir", 1e7, "and %.2f per cent"),
                       ("ribavirin", 1e7, "and %.2f per cent")):
    says(("%.2f × 10^%d^ " + tail)
         % (pv(arm, 5) / sci, int(round(math.log10(sci))),
            100 * pv(arm, 5, "mortality_prob")), "day 5 %s" % arm)

# ------------------------------------------------------- 3.8 identifiability
print("\n== 3.8 what the calibration datum identifies ==")
ridge = {num(r["f"]): r for r in rd("identifiability_ridge.csv")}
a, b = ridge[1], ridge[100]
says("mortality from %.3f to %.3f per cent"
     % (num(a["placebo_mortality"]), num(b["placebo_mortality"])), "ridge mortality")
says("renal injury from %.2f to %.2f" % (num(a["K_peak"]), num(b["K_peak"])),
     "ridge renal")
says("pulmonary injury from %.1f to %.1f" % (num(a["L_peak"]), num(b["L_peak"])),
     "ridge pulmonary")
says("nadir from {:,.0f} to {:,.0f}".format(num(a["PLT_nadir"]), num(b["PLT_nadir"])),
     "ridge platelet nadir")
says("dialysis risk from %.2f to %.2f per cent"
     % (num(a["dialysis_risk"]), num(b["dialysis_risk"])), "ridge dialysis")
says("falls from %.2f × 10^7^ to %.2f × 10^5^ copies"
     % (num(a["V_peak"]) / 1e7, num(b["V_peak"]) / 1e5), "ridge viral peak")
says("reduction from %.1f to %.1f per cent"
     % (num(a["combo_day1_RRR"]), num(b["combo_day1_RRR"])), "ridge RRR")

prof = [r for r in rd("identifiability_profile.csv")
        if r["consistent"] == "TRUE" and r["combo_RRR"] not in ("", "NA")]
by = {}
for r in prof:
    by.setdefault(r["label"], []).append(num(r["combo_RRR"]))
says("between %.1f and %.1f per cent"
     % (min(by["Viral production rate (p)"]), max(by["Viral production rate (p)"])),
     "viral production span")
says("constant between %.1f and %.1f per cent"
     % (min(by["Infection rate (beta)"]), max(by["Infection rate (beta)"])),
     "infection rate span")
allr = [x for v in by.values() for x in v]
import math
says("benefits from %d to %d per cent" % (math.floor(min(allr)), math.ceil(max(allr))),
     "overall span")

# ------------------------------------------------------------- 3.9 sensitivity
# the joint region: the same criterion applied to the global sensitivity
# designs, in which all 32 parameters vary together
J = {r["quantity"]: num(r["value"]) for r in rd("identifiability_joint.csv")}
says("interval runs from %.2f to %.2f per cent" % (J["ci_lower_pct"], J["ci_upper_pct"]),
     "binomial interval")
says("rejects %d of the %d designs" % (J["n_rejected"], J["n_established"]),
     "designs rejected")
says("remaining %d the predicted" % J["n_consistent"], "designs retained")
says("runs from %.1f to %.1f per cent" % (J["rrr_min"], J["rrr_max"]),
     "joint span")
says("interquartile range of %.1f to %.1f and a median of %.1f"
     % (J["rrr_q1"], J["rrr_q3"], J["rrr_median"]), "joint quartiles")

print("\n== 3.9 sensitivity ==")
gsa = io.open(os.path.join(OUT, "gsa_summary.txt"), encoding="utf-8",
              errors="replace").read()
m = re.search(r"n_established .*?:\s+(\d+) \(", gsa)
n_est = int(m.group(1))
m = re.search(r"excluded \(never established\):\s+(\d+)", gsa)
n_exc = int(m.group(1))
says("established on the placebo arm in %d of %d sampled parameter sets"
     % (n_est, n_est + n_exc), "establishment count")
says("undefined for the other %d" % n_exc, "excluded count")
m = re.search(r"mortality \[ALL designs.*?median=([\d.]+)%\s+IQR=\[([\d.]+), ([\d.]+)\]", gsa)
says("median of %s per cent (interquartile range %s to %s) across all"
     % m.groups(), "placebo mortality all")
m = re.search(r"mortality \[ESTABLISHING.*?median=([\d.]+)%\s+IQR=\[([\d.]+), ([\d.]+)\]", gsa)
says("%s per cent (%s to %s) across the" % m.groups(), "placebo mortality subset")

prcc = {(r["outcome"], r["parameter"]): num(r["prcc"]) for r in rd("gsa_prcc.csv")
        if "prcc" in r and r["prcc"] not in ("", "NA")}
out_rrr = [k for k in prcc if "RRR" in k[0]]
oc = out_rrr[0][0] if out_rrr else None
if oc:
    says("viral production (−%.2f), infectivity (−%.2f) and clearance (+%.2f)"
         % (abs(prcc[(oc, "p")]), abs(prcc[(oc, "beta")]), prcc[(oc, "c")]),
         "RRR drivers")
    says("survived adjustment at −%.2f" % abs(prcc[(oc, "EC50_FAV")]),
         "favipiravir EC50")

# -------------------------------------------- supplementary table S3, by hand
# S3 is the one table in the package whose model column is typed rather than
# generated, because the other column is read out of the literature. Check the
# typed half against the same representative-patient record the text uses.
# ------------------------------------------------------------ 4.1 limitations
print("\n== 4.1 limitations ==")
says("It peaks near %d " % round(D["C_RBV_peak_ribavirin_day1"]),
     "ribavirin exposure peak")
says("peaking at %d, %d and %d rather than"
     % (round(D["P_peak"]), round(D["K_peak"]), round(D["L_peak"])),
     "arbitrary-unit peaks")

print("\n== supplementary table S3 ==")
SUP = io.open(os.path.join(HERE, "supplementary.md"), encoding="utf-8").read()
S3 = SUP.split("Table S3.")[1].split("Table S4.")[0]


def s3(fragment, label):
    if fragment in S3:
        print("  ok    %-46s %s" % (label, fragment[:76]))
    else:
        fails.append("S3 " + label)
        print("  FAIL  %-46s not in table S3: %s" % (label, fragment[:76]))


s3("{:,.0f}".format(D["PLT_nadir"]), "platelet nadir")
s3("Day %d to %d" % (int(D["t_PLT_nadir"]), int(D["t_PLT_nadir"]) + 1),
   "platelet nadir timing")
s3("Day %d" % round(D["t_V_peak"]), "viral peak timing")
s3("$%.2f \\times 10^{7}$ copies" % (D["V_peak"] / 1e7), "peak viraemia")
s3("%.2f per cent" % D["mortality_pct"], "placebo mortality")
s3("%.1f per cent" % D["dialysis_pct"], "dialysis requirement")
s3("%.1f g dl" % D["Hgb_drop_ribavirin_day1"], "haemoglobin decline")

print("\n== summary ==")
print("  %d of the checks above failed" % len(fails))
sys.exit(1 if fails else 0)
