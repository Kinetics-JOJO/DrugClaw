---
name: competitive-intelligence-tools
description: Pharmaceutical competitive-intelligence workflow — clinical-trial pipeline mapping from ClinicalTrials.gov v2 API, phase progression analysis, sponsor benchmarking, modality trend extraction, and publication-activity signals from OpenAlex. Use when the user asks about competitor pipelines, therapeutic-area clinical landscapes, or development-stage distribution for a target or indication.
source: drugclaw
updated_at: "2026-03-23"
---

# Competitive Intelligence Tools

Use this skill when the user wants a data-grounded view of who is developing what, where, and at which stage for a target or indication.

Typical triggers:
- map all active clinical trials for a therapeutic target or indication
- compare sponsor pipelines by phase, modality, and status
- identify under-explored combination strategies or patient populations
- generate a competitive-landscape brief before portfolio or investment review
- track phase progression trends over time for a mechanism class

## Environment Check

```bash
which python3 || true
python3 - <<'PY'
mods = ["requests"]
for name in mods:
    try:
        __import__(name)
        print(f"{name}: ok")
    except Exception as exc:
        print(f"{name}: missing ({exc})")
PY
```

If ClinicalTrials.gov or OpenAlex APIs are unreachable, say so before claiming results.

## Bundled Assets

- `templates/competitive_landscape.py`

## Preferred Workflow

1. Confirm scope: target gene symbol, indication term, mechanism class, modality filter, or combination thereof.
2. Pull clinical-trial records from the ClinicalTrials.gov v2 API with structured field extraction.
3. Parse protocol sections: phase, status, sponsor, collaborators, conditions, interventions, enrollment, dates.
4. Aggregate: phase distribution, sponsor ranking, modality breakdown, enrollment totals, year-over-year filing trends.
5. Optionally cross-reference OpenAlex for publication activity by institutions and years.
6. Export a trial-level table plus a landscape summary JSON. Optionally render a markdown brief.
7. All output is point-in-time; always include the retrieval date.

## Quick Start

Target-focused landscape:

```bash
python3 templates/competitive_landscape.py \
  --target "KRAS G12C" \
  --max-trials 200 \
  --output competitive/kras_g12c_trials.csv \
  --summary competitive/kras_g12c_summary.json \
  --brief competitive/kras_g12c_brief.md
```

Indication-focused landscape:

```bash
python3 templates/competitive_landscape.py \
  --indication "non-small cell lung cancer" \
  --phase PHASE3 \
  --status RECRUITING \
  --max-trials 200 \
  --output competitive/nsclc_p3_trials.csv \
  --summary competitive/nsclc_p3_summary.json
```

Combined query with publication signal:

```bash
python3 templates/competitive_landscape.py \
  --target "PD-1" \
  --indication "melanoma" \
  --include-publications \
  --max-trials 100 \
  --max-publications 50 \
  --output competitive/pd1_melanoma.csv \
  --summary competitive/pd1_melanoma.json \
  --brief competitive/pd1_melanoma.md
```

## Output Expectations

Good answers should mention:
- the query scope (target, indication, filters applied)
- total trials retrieved and their phase distribution
- top 10 sponsors ranked by trial count, with their most advanced phase
- modality breakdown when detectable (small molecule, antibody, ADC, cell therapy, gene therapy, vaccine, etc.)
- enrollment total across active trials as a market-size proxy
- year-over-year first-posted trend for the query
- notable gaps: phases, populations, geographies, or combinations with no registrations
- retrieval date and data source
- where the CSV, JSON, and optional brief were saved

## Working Principles

- A registered trial is not a successful trial. Phase counts measure industry activity, not efficacy.
- Sponsor field shows the responsible party, which may be a CRO or academic center rather than the originator company.
- Modality detection from intervention text is heuristic; validate manually for critical decisions.
- ClinicalTrials.gov covers mostly US-registered trials; EU (CTIS/EudraCT), China (ChiCTR), and other registries are not queried.
- Publication counts reflect academic interest, not commercial viability.
- Do not speculate on unpublished pipeline stages or confidential strategies.

## Failure Modes

- Overly broad query (e.g., "cancer"): thousands of trials returned — narrow by phase, condition mesh, or intervention
- API pagination limits: ClinicalTrials.gov v2 returns max 1000 per call; the template handles pagination up to `--max-trials`
- Sponsor name normalization: the same company may appear under different names (abbreviations, subsidiaries)
- Intervention text parsing: combination regimens are a single string; splitting by "+" or "/" is approximate
- Network access unavailable: cannot query any API

## Related Skills

For patent coverage around the same target or compound class, activate `patent-landscape-tools`.
For compound-level data from ChEMBL, BindingDB, or openFDA, activate `pharma-db-tools`.
For clinical-trial protocol design and reporting guidance, activate `clinical-research-tools`.
For target biology and validation background, activate `target-intelligence-tools`.
For drug repurposing hypothesis generation, activate `drug-repurposing-tools`.
