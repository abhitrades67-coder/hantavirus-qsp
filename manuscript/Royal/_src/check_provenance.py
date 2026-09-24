# -*- coding: utf-8 -*-
"""Check that every display item in the package is produced by the archived code.

The data accessibility statement claims that every figure and table in the
article and the supplement is generated from the archived code. This checks
that claim item by item: each display item is mapped to the R script that
writes it and, for tables, to the output file it is built from. The check fails
if a source script or output file is missing, if a script is not run by
R/reproduce_manuscript.R, or if an output is older than the code that writes it.

Usage, from the repository root:
    python "manuscript/Royal/_src/check_provenance.py"
"""
import io
import os
import sys

SRC = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(SRC)
REPO = os.path.dirname(os.path.dirname(PKG))

# item -> (R script that writes it, output file it is built from or copied from)
ITEMS = [
    ("Figure 1",  "R/run_pipeline.R",                 "outputs/viral_kinetics.png"),
    ("Figure 2",  "R/run_pipeline.R",                 "outputs/treatment_window.png"),
    ("Figure 3",  "R/plot_organ_heatmap.R",           "outputs/organ_injury_heatmap.png"),
    ("Figure 4",  "R/duration_start_interaction.R",   "outputs/duration_start_interaction.png"),
    ("Figure 5",  "R/identifiability_profile.R",      "outputs/identifiability_profile.png"),
    ("Figure S1", "R/run_pipeline.R",                 "outputs/pk_profiles.png"),
    ("Figure S2", "R/plot_adaptive_heatmap.R",        "outputs/adaptive_immunity_heatmap.png"),
    ("Figure S3", "R/ablation_analysis.R",
     "manuscript/supplementary/FigureS4_ribavirin_ablation.png"),
    ("Figure S4", "R/preexposure_analysis.R",         "outputs/preexposure_viral_peak.png"),
    ("Figure S5", "R/gsa_prcc.R",                     "outputs/gsa_prcc_tornado.png"),
    ("Figure S6", "R/make_s7.R",                      "outputs/figure_S7_external_validation.png"),

    ("Table 1",   "R/run_pipeline.R",                 "outputs/clinical_endpoints.csv"),
    ("Table 2",   "R/run_pipeline.R",                 "outputs/treatment_window.csv"),
    ("Table 3",   "R/identifiability_profile.R",      "outputs/identifiability_profile.csv"),
    ("Table S1",  "R/run_pipeline.R",                 "outputs/parameter_table.csv"),
    ("Table S2",  "R/preexposure_analysis.R",         "outputs/preexposure_results.csv"),
    ("Table S4",  "R/gsa_prcc.R",                     "outputs/gsa_prcc.csv"),
    ("Table S5",  "R/mechanism_weights.R",            "outputs/mechanism_weights.csv"),
    ("Table S5",  "R/ablation_analysis.R",            "outputs/ablation_summary.csv"),
    ("Table S6",  "R/duration_start_interaction.R",   "outputs/duration_start_data.csv"),
    ("Table S7",  "R/identifiability_and_humoral.R",  "outputs/humoral_scan.csv"),
    ("Table S8",  "R/identifiability_and_humoral.R",  "outputs/identifiability_ridge.csv"),
]

# Table S3 compares the model with published observations. Its model column is
# typed, and check_numbers.py verifies it against this file.
TYPED = [("Table S3", "R/mechanism_weights.R", "outputs/representative_placebo.csv")]

# In-text values that no table carries.
INTEXT = [
    ("section 3.1", "R/mechanism_weights.R",   "outputs/representative_placebo.csv"),
    ("section 3.2", "R/mechanism_weights.R",   "outputs/mechanism_weights.csv"),
    ("section 3.3", "R/plot_organ_heatmap.R",  "outputs/organ_peak_data.csv"),
    ("section 3.4", "R/rebound_diagnostics.R", "outputs/rebound_diagnostics.csv"),
    ("section 3.6", "R/plot_adaptive_heatmap.R", "outputs/adaptive_peak_data.csv"),
    ("section 3.8", "R/identifiability_joint.R", "outputs/identifiability_joint.csv"),
    ("section 3.9", "R/sensitivity_analysis.R", "outputs/sensitivity_local.csv"),
]

fails = []


def bad(m):
    fails.append(m)
    print("  FAIL  " + m)


driver_path = os.path.join(REPO, "R", "reproduce_manuscript.R")
driver = io.open(driver_path, encoding="utf-8").read()

print("\n== every display item traces to a script and an output ==")
for label, script, out in ITEMS + TYPED + INTEXT:
    sp = os.path.join(REPO, script.replace("/", os.sep))
    op = os.path.join(REPO, out.replace("/", os.sep))
    if not os.path.exists(sp):
        bad("%s: %s is missing" % (label, script))
        continue
    if not os.path.exists(op):
        bad("%s: %s is missing" % (label, out))
        continue
    if ('"%s"' % script) not in driver:
        bad("%s: %s is not run by R/reproduce_manuscript.R" % (label, script))
        continue
    if os.path.getmtime(op) < os.path.getmtime(sp):
        bad("%s: %s is older than %s; regenerate it" % (label, out, script))
        continue
    print("  ok    %-10s %-34s <- %s" % (label, out, script))

print("\n== the package figures match their sources ==")
figdir = os.path.join(PKG, "Figures")
srcs = {l: o for l, s, o in ITEMS if l.startswith("Figure")}
for label, out in sorted(srcs.items()):
    png = os.path.join(figdir, label.replace(" ", "_") + ".png")
    op = os.path.join(REPO, out.replace("/", os.sep))
    if not os.path.exists(png):
        bad("%s: not in the package" % label)
    elif os.path.getsize(png) != os.path.getsize(op):
        bad("%s: differs from %s; rerun build.py" % (label, out))
    else:
        print("  ok    %-10s matches %s" % (label, out))

print("\n== summary ==")
print("  %d failures" % len(fails))
sys.exit(1 if fails else 0)
