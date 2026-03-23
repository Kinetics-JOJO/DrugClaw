---
name: active-learning-tools
description: Active learning for drug discovery — uncertainty-guided compound selection, query-by-committee, expected model change, and batch acquisition strategies to minimize experimental cost. Use when the user wants to intelligently select which compounds to test next from an unlabeled library based on model uncertainty rather than random or exhaustive screening.
source: drugclaw
updated_at: "2026-03-23"
---

# Active Learning Tools

Use this skill when the user wants to select the most informative compounds to test next, reducing the number of experiments needed to find hits or build reliable models.

Typical triggers:
- select the next batch of compounds for experimental testing from a large unlabeled library
- use model uncertainty to prioritize which molecules to assay
- balance exploitation (predicted-best) with exploration (most-uncertain) in compound selection
- run iterative active-learning cycles with incoming experimental results
- compare random vs uncertainty-guided selection strategies retrospectively

## Environment Check

```bash
which python3 || true
python3 - <<'PY'
mods = ["numpy", "pandas", "sklearn"]
for name in mods:
    try:
        __import__(name)
        print(f"{name}: ok")
    except Exception as exc:
        print(f"{name}: missing ({exc})")
PY
```

## Bundled Assets

- `templates/active_learning_select.py`

## Preferred Workflow

1. Confirm the labeled (already-tested) and unlabeled (candidate) datasets.
2. Confirm the objective column, feature strategy, and acquisition function.
3. Train a model on the labeled set and predict uncertainty on the unlabeled set.
4. Rank unlabeled compounds by the chosen acquisition strategy.
5. Export a batch of top-ranked candidates for the next experimental round.
6. When new results arrive, append to the labeled set and repeat.

## Quick Start

Uncertainty-based selection:

```bash
python3 templates/active_learning_select.py \
  --labeled tested_compounds.csv \
  --unlabeled candidate_library.csv \
  --smiles-column smiles \
  --id-column compound_id \
  --target-column activity \
  --acquisition uncertainty \
  --batch-size 50 \
  --feature-backend rdkit-morgan \
  --output active_learning/round1_selections.csv \
  --summary active_learning/round1_summary.json
```

Greedy exploitation (predicted-best):

```bash
python3 templates/active_learning_select.py \
  --labeled tested_compounds.csv \
  --unlabeled candidate_library.csv \
  --smiles-column smiles \
  --id-column compound_id \
  --target-column pIC50 \
  --acquisition greedy \
  --batch-size 20 \
  --output active_learning/greedy_picks.csv \
  --summary active_learning/greedy_summary.json
```

Balanced exploration-exploitation:

```bash
python3 templates/active_learning_select.py \
  --labeled tested_compounds.csv \
  --unlabeled candidate_library.csv \
  --smiles-column smiles \
  --id-column compound_id \
  --target-column activity \
  --acquisition balanced \
  --explore-fraction 0.3 \
  --batch-size 50 \
  --output active_learning/balanced_picks.csv \
  --summary active_learning/balanced_summary.json
```

## Output Expectations

Good answers should mention:
- the acquisition strategy and any exploration-exploitation balance parameters
- model type and feature representation used
- how many labeled and unlabeled compounds were processed
- the batch size and selection criteria
- summary statistics of predicted values and uncertainty estimates for the selected batch
- where the selection CSV and summary JSON were saved

## Working Principles

- Active learning reduces experimental cost by focusing on the most informative compounds.
- Uncertainty estimates from ensemble models (e.g., random forest variance) are practical defaults.
- Pure exploitation (greedy) risks missing diverse chemotypes; pure exploration (max-uncertainty) is inefficient.
- Always report the labeled set size — active learning with <20 labeled points is unreliable.
- Each round's selections should be appended to the labeled set after testing for the next iteration.

## Failure Modes

- Too few labeled compounds (<20): model uncertainty estimates are unreliable
- All predictions have similar uncertainty: feature space may be too coarse
- Labeled and unlabeled sets have disjoint chemistry: domain shift undermines selection
- Missing SMILES or features: stop and report invalid rows

## Related Skills

For single-objective Bayesian optimization over continuous parameters, activate `bayesian-optimization-tools`.
For multi-objective compound ranking, activate `multi-objective-optimization-tools`.
For chemistry featurization, activate `chem-tools`.
For bioactivity model training, activate `chem-tools` (bioactivity_predict template).
