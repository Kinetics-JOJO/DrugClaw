---
name: multi-objective-optimization-tools
description: Multi-objective optimization for drug discovery — Pareto-front construction, NSGA-II style ranking, weighted-sum scalarization, and trade-off visualization for balancing competing drug properties (potency, selectivity, ADMET, synthetic accessibility). Use when the user needs to balance multiple objectives simultaneously rather than optimizing a single metric.
source: drugclaw
updated_at: "2026-03-23"
---

# Multi-Objective Optimization Tools

Use this skill when the user has multiple conflicting objectives to balance, such as potency vs selectivity, efficacy vs toxicity, or activity vs synthetic accessibility.

Typical triggers:
- find Pareto-optimal compounds balancing potency and selectivity
- rank a library across ADMET, activity, and docking scores simultaneously
- visualize trade-offs between two or more drug property objectives
- apply NSGA-II-style non-dominated sorting to a compound table
- scalarize multiple objectives with user-specified weights for downstream ranking

## Environment Check

```bash
which python3 || true
python3 - <<'PY'
mods = ["numpy", "pandas"]
for name in mods:
    try:
        __import__(name)
        print(f"{name}: ok")
    except Exception as exc:
        print(f"{name}: missing ({exc})")
PY
```

## Bundled Assets

- `templates/multi_objective_optimize.py`

## Preferred Workflow

1. Identify the objective columns and whether each should be maximized or minimized.
2. Confirm normalization strategy (min-max, z-score, or rank-based).
3. Run non-dominated sorting to identify Pareto fronts.
4. Optionally compute a weighted-sum scalarized score for a single ranking.
5. Export the Pareto-annotated table and a summary of the non-dominated set.
6. Present trade-offs explicitly; do not collapse multi-objective results into a single score without user consent.

## Quick Start

Pareto-front extraction:

```bash
python3 templates/multi_objective_optimize.py \
  --input screening_results.csv \
  --id-column compound_id \
  --objective "potency:maximize" \
  --objective "selectivity:maximize" \
  --objective "tox_score:minimize" \
  --output moo/pareto_results.csv \
  --summary moo/pareto_summary.json
```

Weighted scalarization:

```bash
python3 templates/multi_objective_optimize.py \
  --input screening_results.csv \
  --id-column compound_id \
  --objective "potency:maximize" \
  --objective "clearance:minimize" \
  --weight potency:0.6 \
  --weight clearance:0.4 \
  --normalize minmax \
  --output moo/weighted_rank.csv \
  --summary moo/weighted_rank.json
```

## Output Expectations

Good answers should mention:
- the objectives and their optimization directions
- normalization method applied
- number of Pareto-optimal (rank-1) solutions found
- key trade-off patterns observed
- if scalarization was used, the weights and their justification
- where the annotated table and summary were written

## Working Principles

- Multi-objective problems do not have a single best answer; present the Pareto set, not just one pick.
- Scalarization weights embed subjective preferences — make them explicit.
- Normalization affects rankings — report the method used.
- Do not silently drop rows with missing objective values; report them.
- Pareto dominance analysis does not require ML; it is a deterministic ranking.

## Failure Modes

- All compounds dominated by one: verify data scale and direction
- Missing objective columns: stop and report
- Extreme outliers distorting normalization: consider rank-based normalization
- Too many objectives (>5): Pareto front becomes degenerate; suggest objective clustering

## Related Skills

For single-objective Bayesian optimization, activate `bayesian-optimization-tools`.
For ADMET and chemistry scoring inputs, activate `chem-tools`.
For statistical modeling of property relationships, activate `stat-modeling-tools`.
For virtual screening pipelines, activate `chem-tools` (virtual_screen template).
