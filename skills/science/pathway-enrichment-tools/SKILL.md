---
name: pathway-enrichment-tools
description: Pathway and gene-set enrichment analysis — over-representation analysis (ORA), gene-set enrichment analysis (GSEA), pathway visualization, and functional annotation from public pathway databases (KEGG, Reactome, GO, MSigDB). Use when the user has a gene list or ranked gene list and wants to identify enriched biological pathways or functional categories.
source: drugclaw
updated_at: "2026-03-23"
---

# Pathway Enrichment Tools

Use this skill when the user has differential expression results, a gene hit list, or a ranked gene list and wants to understand the biological pathways and functions represented.

Typical triggers:
- run over-representation analysis on a list of differentially expressed genes
- perform GSEA on a ranked gene list from RNA-seq or proteomics
- identify enriched KEGG, Reactome, or GO pathways in a target gene set
- compare pathway enrichment across experimental conditions
- generate pathway summary tables for publication or target-selection briefs

## Environment Check

```bash
which python3 || true
python3 - <<'PY'
mods = ["scipy", "pandas", "numpy"]
for name in mods:
    try:
        __import__(name)
        print(f"{name}: ok")
    except Exception as exc:
        print(f"{name}: missing ({exc})")
PY
```

Optional enrichment libraries:

```bash
python3 - <<'PY'
optional = ["gseapy", "statsmodels"]
for name in optional:
    try:
        __import__(name)
        print(f"{name}: ok")
    except Exception as exc:
        print(f"{name}: optional, missing ({exc})")
PY
```

## Bundled Assets

- `templates/pathway_enrichment.py`

## Preferred Workflow

1. Accept a gene list (for ORA) or a ranked gene table with scores (for GSEA).
2. Confirm the organism (default: Homo sapiens) and gene ID type (symbol, Ensembl, Entrez).
3. Confirm pathway databases to query (KEGG, Reactome, GO_BP, GO_MF, GO_CC, MSigDB Hallmarks).
4. Run enrichment analysis with multiple-testing correction (Benjamini-Hochberg).
5. Export a ranked pathway table with enrichment scores, p-values, FDR, and gene overlap.
6. Optionally generate a summary of the top enriched pathways for system-prompt context.

## Quick Start

Over-representation analysis:

```bash
python3 templates/pathway_enrichment.py \
  --input deg_list.txt \
  --mode ora \
  --organism "Homo sapiens" \
  --gene-sets KEGG_2021_Human Reactome_2022 GO_Biological_Process_2023 \
  --background-size 20000 \
  --fdr-cutoff 0.05 \
  --output enrichment/ora_results.csv \
  --summary enrichment/ora_summary.json
```

GSEA from ranked list:

```bash
python3 templates/pathway_enrichment.py \
  --input ranked_genes.csv \
  --gene-column gene_symbol \
  --rank-column log2fc \
  --mode gsea \
  --organism "Homo sapiens" \
  --gene-sets MSigDB_Hallmark_2020 KEGG_2021_Human \
  --permutations 1000 \
  --output enrichment/gsea_results.csv \
  --summary enrichment/gsea_summary.json
```

## Output Expectations

Good answers should mention:
- number of input genes and how many mapped to pathway databases
- pathway databases queried
- number of significantly enriched pathways (FDR < cutoff)
- top enriched pathways with direction (up/down for GSEA)
- multiple-testing correction method
- where output files were saved

## Working Principles

- ORA requires a clear threshold for the gene list (e.g., FDR < 0.05 and |log2FC| > 1); do not mix unfiltered lists with ORA.
- GSEA uses all genes ranked by a metric; do not pre-filter for GSEA.
- Background gene set matters for ORA; use the detected transcriptome, not the full genome, when possible.
- Multiple-testing correction is mandatory; do not report uncorrected p-values as final results.
- Pathway databases have version-specific gene-set compositions; report the version used.
- Enrichment is statistical association, not proof of pathway activation.

## Failure Modes

- Gene ID mismatch (symbols vs Ensembl): low mapping rate, unreliable results
- No enriched pathways at FDR < 0.05: lower the cutoff or check input quality
- gseapy not installed: fall back to scipy-based Fisher's exact test for ORA
- Too few genes (<5) in input list: ORA is underpowered

## Related Skills

For upstream omics processing (differential expression), activate `omics-tools`.
For network-level pharmacology analysis, activate `network-pharmacology-tools`.
For knowledge-graph integration of pathway results, activate `knowledge-graph-tools`.
For target intelligence summaries, activate `target-intelligence-tools`.
