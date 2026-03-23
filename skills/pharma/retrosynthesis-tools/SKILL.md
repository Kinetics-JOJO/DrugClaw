---
name: retrosynthesis-tools
description: Retrosynthetic analysis for drug discovery — disconnection strategy, synthetic-route scoring, reaction-class annotation, and building-block availability checks using RDKit reaction SMARTS and public catalogs. Use when the user asks how to synthesize a target molecule, needs a retrosynthetic tree, or wants to evaluate synthetic accessibility of hit compounds.
source: drugclaw
updated_at: "2026-03-23"
---

# Retrosynthesis Tools

Use this skill when the user wants to explore how a target molecule can be synthesized, assess the synthetic feasibility of virtual-screening hits, or plan a synthesis campaign.

Typical triggers:
- suggest synthetic routes for a target compound
- score synthetic accessibility (SA score) for a set of molecules
- generate a retrosynthetic disconnection tree for a drug candidate
- check building-block availability for proposed synthesis routes
- compare route complexity across a compound series

## Environment Check

```bash
which python3 || true
python3 - <<'PY'
mods = ["rdkit", "pandas", "numpy"]
for name in mods:
    try:
        __import__(name)
        print(f"{name}: ok")
    except Exception as exc:
        print(f"{name}: missing ({exc})")
PY
```

## Bundled Assets

- `templates/retrosynthesis_analyze.py`

## Preferred Workflow

1. Accept target molecule(s) as SMILES, CSV, or drawn structure.
2. Compute synthetic accessibility scores (SA score, SCScore if available).
3. Apply retrosynthetic disconnection rules using reaction SMARTS templates.
4. Score proposed routes by step count, known-reaction coverage, and building-block availability.
5. Export a ranked route table and per-molecule SA score summary.
6. Treat all routes as computational proposals; wet-lab feasibility requires chemist review.

## Quick Start

SA score batch computation:

```bash
python3 templates/retrosynthesis_analyze.py \
  --input hits.csv \
  --smiles-column smiles \
  --id-column compound_id \
  --mode sa-score \
  --output retro/sa_scores.csv \
  --summary retro/sa_scores.json
```

Retrosynthetic disconnection:

```bash
python3 templates/retrosynthesis_analyze.py \
  --smiles "CC(=O)Nc1ccc(O)cc1" \
  --mode disconnection \
  --max-depth 3 \
  --output retro/paracetamol_routes.json \
  --summary retro/paracetamol_summary.json
```

Building-block check:

```bash
python3 templates/retrosynthesis_analyze.py \
  --input retro/paracetamol_routes.json \
  --mode bb-check \
  --catalog building_blocks.csv \
  --catalog-smiles-column smiles \
  --output retro/bb_availability.csv \
  --summary retro/bb_availability.json
```

## Output Expectations

Good answers should mention:
- input molecule(s) and representation
- SA score distribution (mean, median, range)
- number of disconnection routes explored
- route depth and step count
- building-block coverage percentage when catalog is provided
- where output files were saved

## Working Principles

- SA scores are heuristic estimates, not synthesis predictions.
- Rule-based retrosynthesis covers common reaction classes but misses novel transformations.
- Building-block catalogs should be current; stale catalogs give false availability signals.
- Multi-step routes with all commercially available starting materials are preferred over shorter routes requiring custom intermediates.
- Always note whether the analysis is rule-based, template-based, or ML-based.

## Failure Modes

- Complex natural-product-like scaffolds: rule coverage may be insufficient
- SMILES parsing failure: report invalid molecules
- No building-block catalog provided: cannot assess commercial availability
- RDKit not installed: cannot compute SA scores or apply reaction SMARTS

## Related Skills

For ADMET and chemistry property screening, activate `chem-tools`.
For molecular generation of novel analogs, activate `molecular-generation-tools`.
For docking validation of synthesizable hits, activate `docking-tools`.
For compound database lookups, activate `pharma-db-tools`.
