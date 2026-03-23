---
name: pharmacovigilance-tools
description: Pharmacovigilance and safety signal analysis — adverse-event database querying (FAERS/openFDA), disproportionality analysis (PRR, ROR, EBGM), signal detection, and safety-profile comparison for post-market drug surveillance research. Use when the user wants to analyze drug safety signals, compare adverse-event profiles, or prepare pharmacovigilance briefs from public data.
source: drugclaw
updated_at: "2026-03-23"
---

# Pharmacovigilance Tools

Use this skill when the user wants to analyze drug safety data, detect adverse-event signals, or compare safety profiles across drugs from public post-market surveillance databases.

Typical triggers:
- query adverse events for a drug from FAERS/openFDA
- compute disproportionality metrics (PRR, ROR) for a drug-event pair
- compare adverse-event profiles of two drugs in the same class
- generate a safety signal summary report for a target compound
- screen for unexpected safety signals in a therapeutic class

## Environment Check

```bash
which python3 || true
python3 - <<'PY'
mods = ["requests", "pandas", "numpy"]
for name in mods:
    try:
        __import__(name)
        print(f"{name}: ok")
    except Exception as exc:
        print(f"{name}: missing ({exc})")
PY
```

## Bundled Assets

- `templates/pharmacovigilance_analyze.py`

## Preferred Workflow

1. Confirm the drug(s) of interest and the safety question.
2. Query openFDA FAERS API for adverse-event reports.
3. Aggregate event counts by preferred term (MedDRA).
4. Compute disproportionality metrics: PRR, ROR, and optionally EBGM.
5. Flag signals meeting detection thresholds (e.g., PRR ≥ 2, chi² ≥ 4, N ≥ 3).
6. Export a signal table and summary brief.
7. All outputs are research signals, not confirmed causal associations.

## Quick Start

Single drug adverse-event profile:

```bash
python3 templates/pharmacovigilance_analyze.py \
  --drug "ibuprofen" \
  --source openfda \
  --max-records 1000 \
  --output pv/ibuprofen_events.csv \
  --summary pv/ibuprofen_summary.json
```

Disproportionality analysis:

```bash
python3 templates/pharmacovigilance_analyze.py \
  --drug "ibuprofen" \
  --source openfda \
  --mode disproportionality \
  --prr-threshold 2.0 \
  --min-count 3 \
  --output pv/ibuprofen_signals.csv \
  --summary pv/ibuprofen_signals.json
```

Two-drug comparison:

```bash
python3 templates/pharmacovigilance_analyze.py \
  --drugs "ibuprofen,naproxen" \
  --source openfda \
  --mode compare \
  --output pv/nsaid_comparison.csv \
  --summary pv/nsaid_comparison.json
```

## Output Expectations

Good answers should mention:
- drug(s) queried and data source
- number of adverse-event reports retrieved
- top adverse events by frequency
- disproportionality metrics (PRR, ROR) for flagged signals
- signal detection thresholds applied
- explicit caveat that FAERS reports are spontaneous and do not prove causation
- where output files were saved

## Working Principles

- FAERS is a spontaneous reporting system; report counts reflect reporting behavior, not incidence rates.
- Disproportionality does not equal causality; it identifies statistical signals for further investigation.
- PRR and ROR require sufficient background counts to be meaningful; small-N signals are noisy.
- Drug name normalization matters: brand names, generic names, and combination products create mapping challenges.
- Always report the query date; FAERS data is updated quarterly.
- This tool supports research pharmacovigilance, not regulatory submission-grade signal detection.

## Failure Modes

- Drug name not found in FAERS: try alternative names (brand, generic, active ingredient)
- openFDA API rate limits (40 requests/minute without key): throttle and report
- Low report counts (<100 total): disproportionality metrics are unreliable
- MedDRA preferred term granularity: high-level terms may mask specific signals
- Network access unavailable: cannot query openFDA

## Related Skills

For clinical-trial data and study design, activate `clinical-research-tools`.
For drug regulatory lookups, activate `pharma-db-tools`.
For statistical modeling of safety data, activate `stat-modeling-tools`.
For medical QMS documentation, activate `medical-qms-tools`.
