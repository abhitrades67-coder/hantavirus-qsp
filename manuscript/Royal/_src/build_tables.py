"""Generate every data table in the manuscript and supplement directly from
outputs/, then substitute them into the Markdown sources.

Nothing in the submission is typed by hand twice: each table below is built
from the CSV or text file the analysis wrote, so a table cannot drift from the
code that produced it. Placeholders in the .md.in sources look like
{{TABLE_NAME}} and are replaced here.

Usage, from the repository root:
    python "manuscript/Royal/_src/build_tables.py"
"""
import csv
import io
import os
import re
from math import floor, log10

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(PKG))
OUT = os.path.join(REPO, "outputs")

NL = chr(10)
ARMS = ["placebo", "ribavirin", "favipiravir", "combination"]

# Placebo mortality of the representative patient, printed in the Table S5
# stub head. Taken from the ridge scan at f = 1, which is that patient.
MORTALITY_PLACEBO = None  # filled in main()


def rd(name):
    with io.open(os.path.join(OUT, name), encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def sci(x, sig=3):
    """3.41e7 -> '3.41 x 10^7^' in the Markdown superscript style used here."""
    x = float(x)
    if x == 0:
        return "0"
    e = int(floor(log10(abs(x))))
    m = x / 10 ** e
    return "%.*f × 10^%d^" % (sig - 1, m, e)


def pct(x, dp=1):
    return "%.*f" % (dp, 100 * float(x))


# ---------------------------------------------------------------- main tables

def table_endpoints():
    """Clinical endpoints by arm, from the virtual trial."""
    ce = {r["arm"]: r for r in rd("clinical_endpoints.csv")}
    rows = [
        ("Patients (*n*)", lambda r: r["n"]),
        ("Peak viral load (copies ml^−1^)", lambda r: sci(r["V_peak_median"])),
        ("Viral load integral (copies day ml^−1^)", lambda r: sci(r["V_AUC_median"])),
        ("Viral clearance by day 21 (%)", lambda r: pct(r["clearance_rate"])),
        ("Platelet nadir (µl^−1^)", lambda r: "{:,.0f}".format(float(r["PLT_nadir_median"]))),
        ("Peak renal injury (AU)", lambda r: "%.1f" % float(r["K_peak_median"])),
        ("Peak pulmonary injury (AU)", lambda r: "%.1f" % float(r["L_peak_median"])),
        ("Dialysis risk, mean (%)", lambda r: pct(r["dialysis_risk_mean"])),
        ("Extracorporeal support risk, mean (%)", lambda r: pct(r["ecmo_risk_mean"])),
        ("Mortality, mean (%)", lambda r: pct(r["mortality_risk_mean"])),
        ("Mortality, median (%)", lambda r: pct(r["mortality_median"])),
    ]
    out = ["| Endpoint | Placebo | Ribavirin | Favipiravir | Combination |",
           "|---|---|---|---|---|"]
    for label, f in rows:
        out.append("| %s | %s |" % (label, " | ".join(f(ce[a]) for a in ARMS)))
    return "\n".join(out)


def table_window():
    """Mortality and relative risk reduction by regimen and start day.

    Read from outputs/treatment_window.csv rather than the text summary,
    which rounds mortality to one decimal place and would turn 0.96 into 1.00.
    """
    rows = rd("treatment_window.csv")
    plac = 100 * float(rows[0]["placebo_mortality"])
    get = lambda a, d: next(r for r in rows if r["arm"] == a
                            and int(float(r["treatment_day"])) == d)
    days = sorted({int(float(r["treatment_day"])) for r in rows})
    out = ["| Start day | Ribavirin | RRR | Favipiravir | RRR | Combination | RRR |",
           "|---|---|---|---|---|---|---|"]
    for d in days:
        cells = []
        for a in ("ribavirin", "favipiravir", "combination"):
            r = get(a, d)
            cells += ["%.2f" % (100 * float(r["mortality_prob"])),
                      "%.1f" % float(r["RRR"])]
        out.append("| %d | %s |" % (d, " | ".join(cells)))
    return "\n".join(out), plac


def table_profile():
    """What the calibration datum identifies: RRR span per parameter."""
    prof = [r for r in rd("identifiability_profile.csv")
            if r["consistent"] == "TRUE" and r["combo_RRR"] not in ("", "NA")]
    by = {}
    for r in prof:
        by.setdefault((r["parameter"], r["label"]), []).append(float(r["combo_RRR"]))
    rows = []
    for (param, label), v in by.items():
        rows.append((label, param, len(v), min(v), max(v), max(v) - min(v)))
    rows.sort(key=lambda t: -t[5])
    out = ["| Parameter | Values consistent with the trial (*n*) | Lowest predicted RRR (%) | Highest predicted RRR (%) | Span (percentage points) |",
           "|---|---|---|---|---|"]
    for label, param, n, lo, hi, span in rows:
        out.append("| %s | %d | %.1f | %.1f | %.1f |" % (label, n, lo, hi, span))
    return "\n".join(out)


# -------------------------------------------------------- supplementary tables

def table_parameters():
    rows = rd("parameter_table.csv")
    out = ["| Parameter | Value | Units | Description | Source | Confidence |",
           "|---|---|---|---|---|---|"]
    for r in rows:
        val = r["value"]
        try:
            f = float(val)
            val = ("%g" % f) if 1e-3 <= abs(f) < 1e5 or f == 0 else sci(f, 3)
        except ValueError:
            pass
        out.append("| %s | %s | %s | %s | %s | %s |" % (
            r["parameter"], val, r["units"], r["description"], r["source"],
            r["confidence"].capitalize()))
    return "\n".join(out), len(rows)


def table_pep():
    rows = rd("preexposure_results.csv")
    days = sorted({int(float(r["pep_day"])) for r in rows})
    get = lambda a, d: next(r for r in rows
                            if r["arm"] == a and int(float(r["pep_day"])) == d)
    out = ["| Day | " + " | ".join("*V*~peak~, %s" % a for a in ARMS) +
           " | " + " | ".join("Mortality, %s" % a for a in ARMS) + " |",
           "|" + "---|" * (1 + 2 * len(ARMS))]
    for d in days:
        vp, mo = [], []
        for a in ARMS:
            r = get(a, d)
            v = float(r["V_peak"])
            vp.append(("{:,.0f}*".format(v)) if v < 1e4 else sci(v))
            mo.append("%.4f%%" % (100 * float(r["mortality_prob"])))
        out.append("| %d | %s | %s |" % (d, " | ".join(vp), " | ".join(mo)))
    return "\n".join(out), len(rows)


def table_gsa():
    rows = rd("gsa_prcc.csv")
    blocks, n = [], None
    for outcome in sorted({r["outcome"] for r in rows}):
        sel = [r for r in rows if r["outcome"] == outcome and r["sig"] == "*"]
        sel.sort(key=lambda r: -abs(float(r["prcc"])))
        if not sel:
            continue
        n = sel[0]["n"]
        blocks.append("*%s* (%s; %s of the sampled parameters significant after "
                      "adjustment, *n* = %s designs)\n" %
                      (outcome, sel[0]["method"], len(sel), n))
        blocks.append("| Parameter | Coefficient | Adjusted *p* |")
        blocks.append("|---|---|---|")
        for r in sel:
            blocks.append("| %s | %+.3f | %.2g |" %
                          (r["parameter"], float(r["prcc"]), float(r["p_adj"])))
        blocks.append("")
    return "\n".join(blocks), n


def table_ridge():
    rows = rd("identifiability_ridge.csv")
    out = ["| *f* | Peak viral load (copies ml^−1^) | Placebo mortality (%) | Peak renal injury | Peak pulmonary injury | Platelet nadir (µl^−1^) | Dialysis risk (%) | Day-1 combination RRR (%) |",
           "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        out.append("| %g | %s | %.3f | %.2f | %.1f | %s | %.2f | %.1f |" % (
            float(r["f"]), sci(r["V_peak"]), float(r["placebo_mortality"]),
            float(r["K_peak"]), float(r["L_peak"]),
            "{:,.0f}".format(float(r["PLT_nadir"])),
            float(r["dialysis_risk"]), float(r["combo_day1_RRR"])))
    return "\n".join(out)


def table_humoral():
    rows = rd("humoral_scan.csv")
    out = ["| Antibody scaling | Placebo mortality (%) | Day-1 mortality (%) | Day-2 mortality (%) | Day 1 minus day 2 (points) | Antibody share of clearance, day 1 (%) | Antibody share, day 2 (%) |",
           "|---|---|---|---|---|---|---|"]
    for r in rows:
        s = float(r["antibody_scaling"])
        lab = "×1 (as calibrated)" if s == 1 else "×%g" % s
        f = lambda v: ("%.2f" % v) if v < 10 else ("%.1f" % v)
        out.append("| %s | %.3f | %.3f | %.3f | %+.2f | %s | %s |" % (
            lab, float(r["placebo_mortality"]), float(r["day1_mortality"]),
            float(r["day2_mortality"]), float(r["day1_minus_day2_points"]),
            f(float(r["antibody_share_day1"])), f(float(r["antibody_share_day2"]))))
    return "\n".join(out)


def table_mechanisms():
    """Table S5: what carries the output, and what the ribavirin ablations cost.

    Rate-term and mortality shares come from R/mechanism_weights.R; the
    ablation rows from R/ablation_analysis.R. Nothing here is typed by hand.
    """
    mech = rd("mechanism_weights.csv")
    abl = {r["scenario"]: 100 * float(r["mortality_prob"])
           for r in rd("ablation_summary.csv")}

    def sig(v):
        """Two decimals at or above 1 per cent, else two significant figures.

        Written out rather than left to %g so that a trailing zero survives:
        0.090455 reads as 0.090, which is two significant figures, not 0.09.
        """
        v = float(v)
        if v >= 1:
            return "%.2f" % v
        from math import floor, log10
        return "%.*f" % (-int(floor(log10(v))) + 1, v)

    ITALIC = {"Non-specific (c)": "Non-specific ($c$)",
              "CD8+ killing": "CD8$^+$ killing"}
    out = ["| Quantity | Decomposition | Contribution |", "|---|---|---|"]
    seen = set()
    for r in mech:
        q = r["quantity"]
        head = q if q not in seen else ""
        if q == "Mortality probability" and head:
            head = "Mortality probability (%.5f)" % (
                sum(float(x["percent"]) for x in mech
                    if x["quantity"] == q) / 100 * MORTALITY_PLACEBO)
        seen.add(q)
        comp = ITALIC.get(r["component"], r["component"])
        out.append("| %s | %s | %s%% |" % (head, comp, sig(r["percent"])))

    for arm, label, rows in (
            ("ribavirin", "Ribavirin day-1 mortality",
             [("Without immunomodulatory effect", "no_immuno"),
              ("Without endothelial-protective effect", "no_endothelial"),
              ("Without organ-recovery enhancement", "no_recovery"),
              ("Without all three", "no_nonantiviral")]),
            ("combination", "Combination day-1 mortality",
             [("Without immunomodulatory effect", "no_immuno"),
              ("Without endothelial-protective effect", "no_endothelial"),
              ("Without organ-recovery enhancement", "no_recovery"),
              ("Without all three", "no_rbv_nonantiviral")])):
        base = abl["%s_day1_default" % arm]
        head = "%s (%.3f%%)" % (label, base)
        for text, key in rows:
            out.append("| %s | %s | %.3f%% |"
                       % (head, text, abl["%s_day1_%s" % (arm, key)]))
            head = ""
    return NL.join(out)


def table_duration():
    rows = rd("duration_start_data.csv")
    durs = sorted({int(float(r["duration_days"])) for r in rows})
    out = ["| Regimen | Start day | " + " | ".join("%d days" % d for d in durs) + " |",
           "|---|---|" + "---|" * len(durs)]
    for arm in ("ribavirin", "favipiravir", "combination"):
        starts = sorted({int(float(r["start_day"])) for r in rows if r["arm"] == arm})
        for i, sd in enumerate(starts):
            vals = []
            for d in durs:
                m = [r for r in rows if r["arm"] == arm
                     and int(float(r["start_day"])) == sd
                     and int(float(r["duration_days"])) == d]
                vals.append(100 * float(m[0]["mortality_prob"]) if m else float("nan"))
            best = min(range(len(vals)), key=lambda k: vals[k])
            cells = ["**%.2f**" % v if k == best else "%.2f" % v
                     for k, v in enumerate(vals)]
            out.append("| %s | %d | %s |" %
                       (arm.capitalize() if i == 0 else "", sd, " | ".join(cells)))
    placebo = 100 * float(rows[0]["placebo_mortality"])
    return "\n".join(out), placebo


def main():
    param_md, n_param = table_parameters()
    pep_md, n_pep = table_pep()
    gsa_md, gsa_n = table_gsa()
    dur_md, dur_placebo = table_duration()
    global MORTALITY_PLACEBO
    MORTALITY_PLACEBO = float([r for r in rd("identifiability_ridge.csv")
                               if float(r["f"]) == 1][0]["placebo_mortality"]) / 100
    win_md, win_placebo = table_window()

    assert n_pep == 24, "expected 24 prophylaxis scenarios, found %d" % n_pep
    assert gsa_md, "global sensitivity table is empty"

    subs = {
        "TABLE_ENDPOINTS": table_endpoints(),
        "TABLE_WINDOW": win_md,
        "TABLE_PROFILE": table_profile(),
        "TABLE_PARAMETERS": param_md,
        "TABLE_PEP": pep_md,
        "TABLE_GSA": gsa_md,
        "TABLE_RIDGE": table_ridge(),
        "TABLE_HUMORAL": table_humoral(),
        "TABLE_DURATION": dur_md,
        "TABLE_MECHANISMS": table_mechanisms(),
        "N_PARAM": str(n_param),
        "N_GSA": str(gsa_n),
        "DUR_PLACEBO": "%.2f" % dur_placebo,
        "WIN_PLACEBO": "%.2f" % win_placebo,
    }

    # Word counts cannot be known until the .docx files exist, so they are
    # deferred to finalise_counts.py, which runs after the first build.
    DEFERRED = {"TOTAL_WORDS", "TITLEPAGE_WORDS", "BODY_WORDS", "ABSTRACT_WORDS"}

    import datetime
    subs["LETTER_DATE"] = datetime.date.today().strftime("%-d %B %Y")         if os.name != "nt" else datetime.date.today().strftime("%d %B %Y").lstrip("0")

    man = io.open(os.path.join(HERE, "manuscript.md.in"), encoding="utf-8").read()
    subs["N_REFS"] = str(len(re.findall(r"^(\d+)\.\s", man.split("## References")[1]
                                        .split("## Tables")[0], re.M)))

    for src in ("manuscript.md.in", "supplementary.md.in", "title_page.md.in",
                "cover_letter.md.in", "readme.md.in"):
        spath = os.path.join(HERE, src)
        if not os.path.exists(spath):
            print("  SKIP %s (not written yet)" % src)
            continue
        text = io.open(spath, encoding="utf-8").read()
        for key, val in subs.items():
            text = text.replace("{{%s}}" % key, val)
        left = set(re.findall(r"\{\{([A-Z_]+)\}\}", text)) - DEFERRED
        assert not left, "unsubstituted placeholders in %s: %s" % (src, sorted(left))
        dst = os.path.join(HERE, src[:-3])
        io.open(dst, "w", encoding="utf-8").write(text)
        print("  OK   %-22s -> %s (%d chars)" % (src, os.path.basename(dst), len(text)))

    print("tables built from outputs/: %d parameters, %d prophylaxis scenarios, "
          "%s global-sensitivity designs" % (n_param, n_pep, gsa_n))


if __name__ == "__main__":
    main()
