---
name: molecular-generation-tools
description: De novo molecular generation and enumeration — scaffold hopping, R-group enumeration, matched molecular pair analysis, SMILES-based library design, and diversity-filtered analog generation using RDKit. Use when the user wants to generate novel molecules, expand a hit series, or design focused libraries around a lead compound.
source: drugclaw
updated_at: "2026-03-23"
---

# Molecular Generation Tools

Use this skill when the user wants to computationally generate new molecular structures for drug discovery, expand a hit series, or build combinatorial libraries.

Typical triggers:
- enumerate R-group variations around a core scaffold
- generate analogs of a hit compound via matched molecular pair transforms
- perform scaffold hopping from one chemotype to another
- build a focused combinatorial library from a set of building blocks and reaction schemes
- filter generated molecules by novelty, diversity, and drug-likeness

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

- `templates/molecular_generate.py`

## Preferred Workflow

1. Accept seed molecule(s) as SMILES, CSV, or inline structures.
2. Confirm the generation strategy: R-group enumeration, MMP transforms, scaffold hopping, or combinatorial.
3. Apply the generation rules and collect candidate structures.
4. Filter by validity (parseable SMILES), drug-likeness (Lipinski/Veber), and diversity (Tanimoto cutoff).
5. Deduplicate by canonical SMILES.
6. Export the library as CSV with computed properties.
7. Treat generated molecules as design hypotheses requiring computational and experimental validation.

## Quick Start

R-group enumeration:

```bash
python3 templates/molecular_generate.py \
  --mode rgroup \
  --core "c1ccc([*:1])cc1" \
  --rgroups '[["F", "Cl", "Br", "OC", "N(C)C"]]' \
  --output molgen/rgroup_library.csv \
  --summary molgen/rgroup_summary.json
```

Matched molecular pair transforms:

```bash
python3 templates/molecular_generate.py \
  --mode mmp \
  --input hits.csv \
  --smiles-column smiles \
  --id-column compound_id \
  --output molgen/mmp_analogs.csv \
  --summary molgen/mmp_summary.json
```

Diversity-filtered library from seeds:

```bash
python3 templates/molecular_generate.py \
  --mode diversify \
  --input seed_compounds.csv \
  --smiles-column smiles \
  --id-column compound_id \
  --diversity-cutoff 0.4 \
  --max-analogs 500 \
  --apply-lipinski \
  --output molgen/diverse_library.csv \
  --summary molgen/diverse_summary.json
```

## Output Expectations

Good answers should mention:
- seed molecule(s) and generation strategy
- number of raw candidates generated
- number passing validity, drug-likeness, and diversity filters
- key property distributions (MW, LogP, HBD, HBA, TPSA)
- novelty relative to input seeds (Tanimoto similarity range)
- where the library CSV and summary JSON were saved

## Working Principles

- Generated molecules are hypotheses, not validated drug candidates.
- Enumerate conservatively; a library of 500 well-filtered compounds is more useful than 50,000 unfiltered ones.
- Always report the Tanimoto similarity distribution between generated and seed molecules.
- Filter out known problematic substructures (PAINS, BRENK) when possible.
- R-group enumeration produces more predictable SAR series; scaffold hopping introduces more novelty but more risk.
- Canonical SMILES deduplication is mandatory before downstream scoring.

## Failure Modes

- Invalid SMILES in seed set: report and skip bad inputs
- Core SMARTS does not match any seed: verify core definition
- Too many enumerated products (combinatorial explosion): limit by diversity or Lipinski filtering
- RDKit not installed: cannot generate or filter molecules

## Related Skills

For synthetic feasibility assessment of generated molecules, activate `retrosynthesis-tools`.
For ADMET triage of generated libraries, activate `chem-tools`.
For docking generated molecules, activate `docking-tools`.
For multi-objective ranking of generated compounds, activate `multi-objective-optimization-tools`.
