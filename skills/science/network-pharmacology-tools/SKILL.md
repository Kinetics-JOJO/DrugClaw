---
name: network-pharmacology-tools
description: Network pharmacology analysis — multi-target drug-pathway mapping, polypharmacology profiling, drug-target-disease network construction, and network topology metrics for understanding drug mechanisms across biological systems. Use when the user wants to analyze how a drug or compound interacts with multiple targets across interconnected pathways.
source: drugclaw
updated_at: "2026-03-23"
---

# Network Pharmacology Tools

Use this skill when the user wants to understand the network-level mechanisms of drug action, analyze polypharmacology, or map drug-target-pathway-disease interactions.

Typical triggers:
- build a drug-target-pathway-disease interaction network for a compound or drug
- identify key targets and pathways in a drug's mechanism of action
- analyze polypharmacology profile of a multi-target compound
- compare network topology of different drugs in the same therapeutic area
- find synergistic drug combination targets through network overlap

## Environment Check

```bash
which python3 || true
python3 - <<'PY'
mods = ["networkx", "pandas", "requests"]
for name in mods:
    try:
        __import__(name)
        print(f"{name}: ok")
    except Exception as exc:
        print(f"{name}: missing ({exc})")
PY
```

## Bundled Assets

- `templates/network_pharmacology.py`

## Preferred Workflow

1. Accept drug/compound identifier(s): name, SMILES, ChEMBL ID, or DrugBank ID.
2. Retrieve known and predicted targets from public databases (ChEMBL, STITCH, SwissTargetPrediction, SEA).
3. Map targets to pathways (KEGG, Reactome) and diseases (OpenTargets, DisGeNET).
4. Construct the multi-layer network: drug → targets → pathways → diseases.
5. Compute network topology metrics: degree, betweenness, closeness, clustering.
6. Identify hub targets and critical pathways.
7. Export the network and analysis results.

## Quick Start

Single drug network:

```bash
python3 templates/network_pharmacology.py \
  --drug "metformin" \
  --source chembl \
  --activity-cutoff 1000 \
  --pathway-db KEGG \
  --output netpharm/metformin_network.graphml \
  --summary netpharm/metformin_summary.json \
  --table netpharm/metformin_targets.csv
```

Multi-drug comparison:

```bash
python3 templates/network_pharmacology.py \
  --drugs metformin,pioglitazone,empagliflozin \
  --source chembl \
  --mode compare \
  --output netpharm/diabetes_drugs_network.graphml \
  --summary netpharm/diabetes_drugs_summary.json \
  --table netpharm/diabetes_drugs_targets.csv
```

Drug combination synergy analysis:

```bash
python3 templates/network_pharmacology.py \
  --drugs "drug_a,drug_b" \
  --mode combination \
  --overlap-analysis \
  --output netpharm/combo_network.graphml \
  --summary netpharm/combo_summary.json
```

## Output Expectations

Good answers should mention:
- drug(s) analyzed and data sources queried
- number of known targets retrieved and activity threshold used
- pathway mapping coverage (how many targets mapped to at least one pathway)
- network size: nodes, edges, connected components
- hub targets (highest degree or betweenness centrality)
- enriched pathways from the target set
- where network files and summaries were saved

## Working Principles

- Network pharmacology reveals potential mechanisms; it does not prove them.
- Activity cutoffs for target retrieval directly affect network construction; always report the threshold.
- Predicted targets (from structure similarity) should be flagged separately from experimentally validated ones.
- Hub targets in the network may reflect database bias (well-studied proteins) rather than biological importance.
- Network topology metrics complement but do not replace pathway enrichment analysis.
- Drug combination network overlap suggests shared biology, not guaranteed synergy.

## Failure Modes

- No targets found for the query: check identifier format, try alternative databases
- Network too sparse (<5 targets): compound may be highly selective or poorly characterized
- Network too dense (>500 targets with low cutoff): increase activity cutoff stringency
- API rate limits on ChEMBL or STITCH: batch and retry
- Predicted targets conflated with experimental: clearly separate confidence levels

## Related Skills

For pathway-level enrichment from gene lists, activate `pathway-enrichment-tools`.
For knowledge-graph construction and querying, activate `knowledge-graph-tools`.
For compound and target database lookups, activate `pharma-db-tools` and `bio-db-tools`.
For target intelligence dossiers, activate `target-intelligence-tools`.
