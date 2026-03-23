---
name: pk-modeling-tools
description: Pharmacokinetic modeling for drug discovery — non-compartmental analysis (NCA), one- and two-compartment PK fitting, exposure metrics (AUC, Cmax, t½, clearance, Vd), dose-response simulation, and allometric scaling. Use when the user has concentration-time data and wants PK parameter estimates or PK/PD simulation outputs.
source: drugclaw
updated_at: "2026-03-23"
---

# PK Modeling Tools

Use this skill when the user has concentration-time data or PK parameters and wants to perform pharmacokinetic analysis, simulate exposure profiles, or estimate dosing.

Typical triggers:
- compute NCA parameters (AUC, Cmax, Tmax, t½, clearance) from concentration-time data
- fit one-compartment or two-compartment PK models to time-course data
- simulate plasma concentration profiles for different dosing regimens
- perform allometric scaling from preclinical species to human
- generate PK summary tables for study reports
- compare exposure metrics across formulations or species

## Environment Check

```bash
which python3 || true
python3 - <<'PY'
mods = ["numpy", "pandas", "scipy"]
for name in mods:
    try:
        __import__(name)
        print(f"{name}: ok")
    except Exception as exc:
        print(f"{name}: missing ({exc})")
PY
```

## Bundled Assets

- `templates/pk_analysis.py`

## Preferred Workflow

1. Confirm data format: concentration-time table with subject/group identifiers.
2. Confirm analysis type: NCA, compartmental fit, or simulation.
3. For NCA: compute AUC (linear-log trapezoidal), Cmax, Tmax, terminal t½, clearance, Vd.
4. For compartmental: fit the selected model, report parameter estimates and goodness-of-fit.
5. For simulation: accept PK parameters and dosing regimen, output predicted concentration-time curves.
6. Export PK parameter tables and optional plots.

## Quick Start

Non-compartmental analysis:

```bash
python3 templates/pk_analysis.py \
  --input pk_data.csv \
  --time-column time_h \
  --conc-column concentration_ng_ml \
  --group-column subject_id \
  --dose 10 \
  --dose-unit mg/kg \
  --mode nca \
  --output pk/nca_results.csv \
  --summary pk/nca_summary.json
```

One-compartment IV bolus fit:

```bash
python3 templates/pk_analysis.py \
  --input pk_data.csv \
  --time-column time_h \
  --conc-column concentration_ng_ml \
  --group-column subject_id \
  --dose 10 \
  --mode compartmental \
  --compartments 1 \
  --route iv-bolus \
  --output pk/comp1_fit.csv \
  --summary pk/comp1_summary.json
```

Dose-response simulation:

```bash
python3 templates/pk_analysis.py \
  --mode simulate \
  --compartments 1 \
  --route oral \
  --ka 1.2 \
  --ke 0.35 \
  --vd 50 \
  --bioavailability 0.6 \
  --dose 100 \
  --dose-unit mg \
  --interval 24 \
  --doses 7 \
  --output pk/simulation.csv \
  --summary pk/simulation.json
```

## Output Expectations

Good answers should mention:
- data source and number of subjects or groups analyzed
- analysis type (NCA, compartmental, simulation)
- key PK parameters with units: AUC, Cmax, Tmax, t½, CL, Vd
- goodness-of-fit metrics for compartmental models (R², RMSE)
- dosing regimen for simulations
- where output files were saved

## Working Principles

- NCA is model-independent and robust for standard PK summary; use it as the default when compartmental modeling is not specifically requested.
- Always report units for every PK parameter.
- Linear-log trapezoidal rule is standard for NCA AUC computation.
- Compartmental model selection should be justified by data (not assumed to be two-compartment).
- Allometric scaling is approximate; species-specific binding and metabolism differences are not captured.
- Simulation outputs are predictions, not experimental measurements.
- PK modeling does not replace clinical pharmacology expertise.

## Failure Modes

- Too few time points for reliable NCA (need ≥3 in terminal phase for t½)
- Concentration below LLOQ: report handling strategy (BLQ = 0, BLQ = LLOQ/2, or excluded)
- Compartmental fit not converging: try different initial estimates or reduce model complexity
- Mixed units in input data: standardize before analysis
- Missing dose information: cannot compute clearance or volume of distribution

## Related Skills

For ADMET property prediction from structure, activate `chem-tools`.
For clinical-trial design and endpoints, activate `clinical-research-tools`.
For survival analysis with time-to-event data, activate `survival-analysis-tools`.
For statistical modeling of dose-response relationships, activate `stat-modeling-tools`.
