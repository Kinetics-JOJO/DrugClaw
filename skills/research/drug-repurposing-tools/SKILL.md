---
name: drug-repurposing-tools
description: Drug repurposing research workflow — systematic identification of new therapeutic indications for existing drugs via target overlap, pathway similarity, transcriptomic signature reversal, and clinical-evidence triangulation from public databases (OpenTargets, ChEMBL, ClinicalTrials.gov, CMap). Use when the user asks about repurposing candidates, off-label evidence, or repositioning strategies for known drugs or compound classes.
source: drugclaw
updated_at: "2026-03-23"
---

# Drug Repurposing Tools

Use this skill when the user wants to systematically evaluate whether an existing approved or investigational drug could work in a new indication.

Typical triggers:
- find repurposing candidates for a disease by searching for drugs that hit related targets
- assess whether a marketed drug has off-label evidence in a different therapeutic area
- compare the target profile of a drug against the target landscape of a disease
- build a repurposing hypothesis brief from public data for portfolio review
- check ClinicalTrials.gov for existing off-label or investigator-initiated trials

## Environment Check

```bash
which python3 || true
python3 - <<'PY'
mods = ["requests", "json", "csv"]
for name in mods:
    try:
        __import__(name)
        print(f"{name}: ok")
    except Exception as exc:
        print(f"{name}: missing ({exc})")
PY
```

If outbound network access is blocked, say so explicitly before claiming results were fetched.

## Bundled Assets

- `templates/drug_repurpose.py`

## Preferred Workflow

1. Start from either a **drug** (seeking new indications) or a **disease** (seeking existing drugs).
2. Resolve the drug to ChEMBL ID and its known targets; resolve the disease to OpenTargets EFO ID and its associated targets.
3. Compute **target overlap** between the drug's mechanism and the disease's genetic and functional target landscape.
4. Check for **existing clinical evidence**: ClinicalTrials.gov registrations, case reports in OpenAlex, or pharmacovigilance signal inversions.
5. Score each drug-disease pair on a composite of target-overlap, genetic-association, and clinical-evidence signals.
6. Export a ranked repurposing brief with provenance for each signal.
7. Treat every output as a research hypothesis requiring clinical validation.

## Quick Start

Disease-centric (find drugs for a disease):

```bash
python3 templates/drug_repurpose.py \
  --mode disease-to-drugs \
  --disease "Crohn's disease" \
  --max-targets 30 \
  --max-drugs 20 \
  --output repurposing/crohn_candidates.csv \
  --summary repurposing/crohn_summary.json \
  --dossier repurposing/crohn_brief.md
```

Drug-centric (find new indications for a drug):

```bash
python3 templates/drug_repurpose.py \
  --mode drug-to-diseases \
  --drug metformin \
  --max-diseases 20 \
  --output repurposing/metformin_indications.csv \
  --summary repurposing/metformin_summary.json \
  --dossier repurposing/metformin_brief.md
```

Pair check (assess a specific drug-disease hypothesis):

```bash
python3 templates/drug_repurpose.py \
  --mode pair-check \
  --drug imatinib \
  --disease "pulmonary arterial hypertension" \
  --output repurposing/imatinib_pah.csv \
  --summary repurposing/imatinib_pah.json \
  --dossier repurposing/imatinib_pah.md
```

## Output Expectations

Good answers should mention:
- the drug and/or disease queried
- how many targets were resolved for the drug and disease respectively
- target overlap count and the overlapping gene symbols
- whether any clinical trials exist for the drug-disease combination
- the composite repurposing score and which signals contributed
- explicit limitations: no animal model evidence, no pharmacovigilance signal reversal computed if data is unavailable
- where the output files were saved

## Working Principles

- Drug repurposing hypotheses are cheap to generate and expensive to validate; quality > quantity.
- Target overlap is necessary but not sufficient — the drug must reach the target at safe human exposures.
- Genetic association evidence (GWAS, OpenTargets score) adds independent support beyond pharmacology.
- Clinical-trial registration for off-label use is the strongest public signal short of published results.
- Never present a repurposing hypothesis as a validated finding.
- Provenance is essential: record every data source and identifier for each evidence line.

## Failure Modes

- Drug not in ChEMBL: limited target information; check DrugBank or PubChem instead
- Disease EFO term not found: try synonyms or parent terms
- No target overlap: does not mean no repurposing potential (mechanism may be indirect)
- ClinicalTrials.gov search too broad: refine by intervention name and condition MeSH term
- Network access unavailable: cannot run any data retrieval

## Related Skills

For raw target data across UniProt, PDB, ClinVar, gnomAD, Reactome, STRING, OpenTargets, activate `bio-db-tools`.
For compound and bioactivity lookups in ChEMBL, BindingDB, openFDA, activate `pharma-db-tools`.
For integrated single-target dossiers, activate `target-intelligence-tools`.
For knowledge-graph-based path analysis between drugs and diseases, activate `knowledge-graph-tools`.
For competitive-landscape context on who else is pursuing this repositioning, activate `competitive-intelligence-tools`.
