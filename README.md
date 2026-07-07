# Hantavirus QSP Model — Ribavirin & Favipiravir

A mechanistic **Quantitative Systems Pharmacology (QSP)** model of Hantavirus infection,
innate immune response, endothelial dysfunction, organ-specific injury, and antiviral
treatment with **Ribavirin** and **Favipiravir**.

## Model Architecture

```
┌──────────────────────────────────────────────────────────────────────┐
│                        VIRAL DYNAMICS                                │
│  dT/dt = -β·V·T                                                     │
│  dI/dt = β·V·T - δ_nat·I - δ_NK·NK·I  (non-cytopathic)             │
│  dV/dt = p·I · IFN_antiviral · (1-E_RBV) · (1-E_FAV) · (1-ψ) - c·V │
└────────────┬───────────────────────────────────┬─────────────────────┘
             │                                   │
             ▼                                   ▼
┌──────────────────────────┐        ┌──────────────────────────┐
│  INNATE IMMUNE (6-var)   │        │  ENDOTHELIAL + PLT       │
│  NSs → IFN antagonist    │        │  dP/dt = k_PV·V          │
│  F_I → Type I IFN (α/β)  │        │        + k_PC·C_pro/Cref │
│  NK  → NK cells          │        │  dPLT = prod-loss-cons   │
│  F_II → Type II IFN (γ)  │        └──────────┬───────────────┘
│  C_pro → Pro-infl cytokines │                 │
│         + auto-amplification│                 ▼
│  C_anti → IL-10 resolution │      ┌─────────────────────────┐
└──────────┬─────────────────┘      │    ORGAN INJURY          │
           │                        │  dK/dt (Renal)           │
           │                        │  dL/dt (Lung/Cardio)     │
           ▼                        └─────────────────────────┘
    IFN_antiviral reduces viral production
           ▲                                ▲
           │                                │
┌──────────┴─────────┐      ┌───────────────┴──────────┐
│  RIBAVIRIN PK/PD   │      │  FAVIPIRAVIR PK/PD       │
│  1-comp IV         │      │  Oral absorption + gut   │
│  Emax PD           │      │  Nonlinear CL + RTP met. │
│                    │      │  Emax PD (via RTP)        │
└────────────────────┘      └──────────────────────────┘
```

### State Variables (23 compartments)

| Symbol | Description | Units |
|--------|-------------|-------|
| T | Susceptible target cells | cells |
| I | Infected cells | cells |
| V | Viral RNA (proxy for viremia) | copies/mL |
| NSs | Hantavirus NSs protein (IFN antagonist) | AU |
| F_I | Type I IFN (IFN-α/β, antiviral) | AU |
| NK | NK cell activation | AU |
| F_II | Type II IFN (IFN-γ, pro-inflammatory) | AU |
| C_pro | Pro-inflammatory cytokine burden | AU |
| C_anti | Anti-inflammatory cytokine burden (IL-10-like) | AU |
| P | Vascular permeability index | 0–1 |
| PLT | Platelet count | /μL |
| K | Renal injury index | 0–1 |
| L | Lung injury index | 0–1 |
| CD8_N | CD8+ naive/primed T cells | AU |
| CD8_E | CD8+ effector T cells | AU |
| CD4 | CD4+ helper T cells | AU |
| IgM | IgM antibodies | AU |
| IgG | Neutralizing IgG antibodies | AU |
| C_RBV | Ribavirin plasma concentration | μg/mL |
| C_FAV_gut | Favipiravir gut compartment | mg |
| C_FAV | Favipiravir plasma concentration | μg/mL |
| C_FAVI_RTP | Favipiravir active metabolite (RTP) | μg/mL |
| Hgb_drop | Hemoglobin drop (ribavirin toxicity) | g/dL |

### Immune Module (v2 — expanded)

The immune module was upgraded from a 2-variable to a 6-variable model:

- **NSs protein**: Hantavirus IFN antagonist that suppresses Type I IFN production
- **Type I IFN (F_I)**: Antiviral — feeds back to reduce viral production via E_max model
- **NK cells**: Activated by infected cells and IFN-I; kill infected cells; produce IFN-γ
- **Type II IFN (F_II)**: Pro-inflammatory — amplifies cytokine production
- **C_pro**: Pro-inflammatory cytokines with Hill-2 auto-amplification (cytokine storm) and IL-10 suppression
- **C_anti**: Anti-inflammatory cytokines (IL-10-like) — delayed resolution

Key features:
- Non-cytopathic infected cell clearance (δ_natural = 0.03/day, t½ ~23 days)
- NK-mediated immune clearance of infected cells
- Cytokine storm via positive feedback (Hill-2 auto-amplification)
- Anti-inflammatory counter-regulation via IL-10
- Unit-normalized production terms using reference values

### Antiviral PD

- **Ribavirin**: Emax model with Emax = 0.70, EC50 = 8 μg/mL, γ = 1.5
- **Favipiravir**: Emax model on RTP metabolite, EC50 = 0.5 μg/mL, γ = 1.2
- **Combination**: Bliss independence + synergy term (ψ = 0.1)

### Clinical Endpoints

- **Dialysis probability**: sigmoidal function of renal injury K
- **ECMO probability**: sigmoidal function of lung injury L
- **Mortality probability**: capped composite of dialysis risk, ECMO risk, and cytokine burden

## Project Structure

```
Hantavirus_QSP/
├── model/
│   ├── hantavirus_qsp.R     # Full 23-compartment ODE model
│   └── parameters.R          # All parameters with sources & confidence
├── R/
│   ├── pk_models.R           # Ribavirin & Favipiravir PK/dosing
│   ├── pd_models.R           # PD, combination, clinical-endpoint functions
│   ├── virtual_population.R  # Virtual patient generation (N=1000)
│   ├── simulate_trial.R      # Virtual trial simulation engine
│   ├── analysis.R            # Plotting & summary functions
│   ├── run_pipeline.R        # Core pipeline: trial + main figures/tables
│   ├── sensitivity_analysis.R / uncertainty_quantification.R / gsa_prcc.R
│   ├── ablation_analysis.R / preexposure_analysis.R / optimal_duration.R / vpc_analysis.R
│   ├── plot_*.R              # Organ/adaptive heatmaps & decluttered trajectories
│   ├── regen_all_aux.R / regen_cached_figs.R   # Re-run auxiliaries / cached figures
│   └── make_s7.R             # Supplementary Fig. S7 (external-corroboration overlay)
├── tests/testthat/           # Unit tests (parameter & model invariants)
├── outputs/                  # Generated figures and tables
├── DESCRIPTION
└── README.md
```

## Installation

### Required R packages

```r
install.packages(c("deSolve", "ggplot2", "dplyr", "tidyr", "gridExtra",
                   "foreach", "doParallel"))
```

### R version

Requires R >= 4.2.0. Tested on R 4.3+.

## Usage

### Run the full pipeline

```bash
cd /path/to/Hantavirus_QSP
Rscript R/run_pipeline.R
```

Or from within R:

```r
setwd("path/to/Hantavirus_QSP")   # set to your local checkout
source("R/run_pipeline.R")
main()
```

`run_pipeline.R` runs the core virtual trial and produces the main figures and
tables. The additional analyses (local/global sensitivity, uncertainty
quantification, mechanism ablation, post-exposure prophylaxis, optimal-duration,
and VPC) live in separate scripts; regenerate them with:

```bash
Rscript R/regen_all_aux.R        # re-runs the auxiliary analyses
Rscript R/regen_cached_figs.R    # regenerates the cached-data figures
Rscript R/make_s7.R              # regenerates Supplementary Fig. S7 (external corroboration)
```

> Note: `R/sensitivity_analysis.R` and `R/uncertainty_quantification.R` must be
> run as standalone `Rscript` invocations (they fail if sourced inside another
> script's environment).

### Regenerate the submission figures

The four main-text figures are produced entirely by code and copied into the
manuscript submission folder by `R/assemble_figures.R`:

```bash
Rscript R/run_pipeline.R           # Figure 1 (viral kinetics) + Figure 2 (treatment window)
Rscript R/plot_organ_heatmap.R     # Figure 3 (organ injury heatmap)
Rscript R/plot_adaptive_heatmap.R  # Figure 4 (adaptive immunity heatmap)
Rscript R/assemble_figures.R       # copy outputs/*.png -> Figures/Figure_1..4.png
```

Figures carry no in-figure title (the caption is supplied in the manuscript).
Figures 1–2 are seeded (`virtual_population` seed 42, `simulate_trial` seed 123),
and the heatmap scripts read git-tracked per-patient caches
(`outputs/{organ,adaptive}_peak_data.csv`), so re-running reproduces the
submitted figures exactly. Delete those caches (or run `R/regen_cached_figs.R`)
to re-simulate the heatmap inputs from scratch.

### Individual components

```r
# Load parameters
source("model/parameters.R")
pars <- get_parameters()

# Generate virtual population
source("R/virtual_population.R")
pop <- generate_virtual_population(N = 1000)

# Run single simulation
source("model/hantavirus_qsp.R")
source("R/simulate_trial.R")
sim <- simulate_patient(pars, arm = "ribavirin", t_start = 3, t_end = 21)

# Plot results
source("R/analysis.R")
plot_viral_kinetics(sim, outfile = "outputs/viral_kinetics.png")
```

## Virtual Trial Design

### Treatment Arms

1. **Placebo** — Supportive care only
2. **Ribavirin** — IV HFRS regimen (Huggins et al.):
   - Loading: 33 mg/kg IV
   - Days 1–4: 16 mg/kg q6h
   - Days 5–10: 8 mg/kg q8h
3. **Favipiravir** — Oral standard regimen:
   - Day 1: 1600 mg BID
   - Days 2–15: 600 mg BID
4. **Combination** — Both drugs at above doses

### Simulation Scenarios

- **Early treatment**: Days 1–3 from symptom onset
- **Mid treatment**: Days 4–5
- **Late treatment**: Days 6–7

### Virtual Population (N = 1000)

| Covariate | Distribution | Range |
|-----------|-------------|-------|
| Age | Normal(45, 15) | [18, 80] |
| Body weight | Normal(75, 12) kg | [40, 120] |
| eGFR | Normal(90, 25) mL/min | [15, 150] |
| Time to treatment | Uniform(1, 7) days | — |
| Initial viral load | LogNormal(log(100), 1) | [10, 10000] |
| Immune strength | LogNormal(0, 0.3) | [0.1, 5] |
| Endothelial sensitivity | LogNormal(0, 0.3) | [0.1, 5] |
| Baseline platelets | Normal(250k, 50k) | [100k, 400k] |
| Syndrome | 60% HFRS / 40% HCPS | — |

## Output Files

| File | Description |
|------|-------------|
| `outputs/viral_kinetics.png` | Viral load trajectories by arm |
| `outputs/biomarker_trajectories.png` | Platelets, cytokines, permeability |
| `outputs/organ_injury.png` | Renal (K) and Lung (L) injury |
| `outputs/pk_profiles.png` | Drug concentration profiles |
| `outputs/treatment_window.png` | Benefit by treatment start day |
| `outputs/clinical_endpoints.csv` | Summary endpoints by arm |
| `outputs/virtual_population.csv` | Virtual patient characteristics |
| `outputs/simulation_summary.txt` | Text summary of findings |
| `outputs/parameter_table.csv` | All parameters with sources |

## Parameter Sources

See `outputs/parameter_table.csv` for the complete parameter table with:
- Parameter name and value
- Units
- Description
- Literature source
- Confidence rating (high / medium / low)

Key references:
- Huggins et al. — Ribavirin in HFRS (J. Infect. Dis., 1991)
- Furuta et al. — Favipiravir (T-705), broad-spectrum viral RdRp inhibitor (Proc. Jpn. Acad. Ser. B, 2017)
- Safronetz et al. — Favipiravir against hantaviruses (Antimicrob. Agents Chemother., 2013)
- Williams et al. — Zoonotic spillover & rodent-borne RNA viruses (Viruses, 2021)
- KDIGO — AKI staging criteria

## Model Limitations

1. **Parameter uncertainty**: Many parameters (especially organ injury rates and immune
   module rates) are assumed due to limited quantitative clinical data in Hantavirus.
2. **Composite cytokines**: C_pro is a composite variable, not individual cytokines
   (IL-6, TNF-α, etc.).
3. **Reduced immune representation**: Adaptive immunity (CD4⁺/CD8⁺ T cells, IgM, IgG) is
   included but simplified — CD8⁺ cytotoxicity uses a two-pool maturation chain rather than
   full clonal/repertoire dynamics, antibodies are lumped into IgM and IgG, and NK cells
   are a single activation variable.
4. **Single-strain model**: Does not distinguish between Hantaan, Dobrava,
   Puumala, Sin Nombre, Andes, etc.
5. **Deterministic ODEs**: Stochastic effects at low viral loads are not captured.
6. **Peak viral load magnitude**: Simulated peak viremia (~10⁷ copies/mL) runs higher than
   typical clinical reports (~10⁴–10⁶), reflecting NSs-mediated near-complete Type I IFN
   suppression. Clinical endpoints derive from downstream injury and cytokine states rather
   than absolute viral titre, so the absolute value of V should not be over-interpreted.

## Version History

- **v2.1** (2026-06): Adaptive immunity, calibration, and public release (software v1.0.0):
  - Added CD4⁺/CD8⁺ T cells (two-pool maturation chain) and IgM/IgG humoral response,
    bringing the model to its current 23 compartments
  - 60-day viral clearance with logistic target-cell regeneration
  - Platelet-baseline recalibration; dialysis/ECMO thresholds decoupled from mortality
  - Calibrated to Huggins JID 1991 (~9.6% placebo mortality); archived at Zenodo
    (DOI 10.5281/zenodo.20558774)
- **v2.0** (2026-05-11): Expanded 6-variable innate immune module:
  - Added NSs protein (IFN antagonist), Type I/II IFN split, NK cells
  - Cytokine auto-amplification (Hill-2) and IL-10 anti-inflammatory resolution
  - IFN feedback on viral production; non-cytopathic infected cell clearance
  - Unit-normalized production terms; sigmoidal cytokine mortality
- **v1.0** (2026-05-08): Initial model (14 compartments):
  - Viral dynamics, 2-variable innate immune (IFN + composite cytokine),
    endothelial permeability, organ injury, Ribavirin/Favipiravir PK/PD

## License

MIT License
