---
name: epigenomics-tools
description: Epigenomics analysis for drug discovery — DNA methylation profiling, ChIP-seq peak analysis, chromatin accessibility (ATAC-seq) summaries, histone mark enrichment, and epigenetic target identification. Use when the user has methylation arrays, ChIP-seq peaks, or ATAC-seq data and wants epigenomic characterization relevant to drug-target discovery.
source: drugclaw
updated_at: "2026-03-23"
---

# Epigenomics Tools

Use this skill when the user has epigenomic data (methylation, ChIP-seq, ATAC-seq) and wants to analyze it in the context of drug discovery — identifying epigenetic drug targets, characterizing chromatin states, or profiling methylation changes.

Typical triggers:
- summarize DNA methylation profiles from array data (450K/EPIC)
- identify differentially methylated regions or CpG sites
- analyze ChIP-seq peak files for histone mark or TF enrichment
- summarize ATAC-seq accessibility at promoters or enhancers
- identify epigenetic drug targets (HDACs, DNMTs, BET proteins, etc.) from epigenomic signatures
- integrate methylation data with gene expression for epigenetic regulation analysis

## Environment Check

```bash
which python3 || true
python3 - <<'PY'
mods = ["pandas", "numpy", "scipy"]
for name in mods:
    try:
        __import__(name)
        print(f"{name}: ok")
    except Exception as exc:
        print(f"{name}: missing ({exc})")
PY
```

Optional genomics libraries:

```bash
python3 - <<'PY'
optional = ["pysam", "pyBigWig"]
for name in optional:
    try:
        __import__(name)
        print(f"{name}: ok")
    except Exception as exc:
        print(f"{name}: optional, missing ({exc})")
PY
```

## Bundled Assets

- `templates/epigenomics_analyze.py`

## Preferred Workflow

1. Confirm data type: methylation beta-values, ChIP-seq peaks (BED/narrowPeak), or ATAC-seq peaks.
2. For methylation: compute differential methylation (DMP/DMR) between groups, annotate CpG context and gene proximity.
3. For ChIP-seq: summarize peak distribution (promoter, enhancer, intergenic), compute enrichment at gene sets.
4. For ATAC-seq: quantify accessibility at regulatory elements, compare open chromatin between conditions.
5. Identify epigenetic targets implied by the data (e.g., H3K27me3 enrichment → EZH2 as potential target).
6. Export summary tables, annotated peak/region lists, and target-implication briefs.

## Quick Start

Methylation differential analysis:

```bash
python3 templates/epigenomics_analyze.py \
  --input methylation_betas.csv \
  --mode methylation \
  --group-column condition \
  --group-a treatment \
  --group-b control \
  --delta-beta-cutoff 0.2 \
  --pvalue-cutoff 0.01 \
  --output epigenomics/dmps.csv \
  --summary epigenomics/methylation_summary.json
```

ChIP-seq peak annotation:

```bash
python3 templates/epigenomics_analyze.py \
  --input h3k27ac_peaks.narrowPeak \
  --mode chipseq \
  --annotation genes.gtf \
  --promoter-window 2000 \
  --output epigenomics/chipseq_annotated.csv \
  --summary epigenomics/chipseq_summary.json
```

ATAC-seq accessibility profile:

```bash
python3 templates/epigenomics_analyze.py \
  --input atac_peaks.bed \
  --mode atacseq \
  --annotation genes.gtf \
  --output epigenomics/atac_profile.csv \
  --summary epigenomics/atac_summary.json
```

## Output Expectations

Good answers should mention:
- data type and number of features analyzed (CpGs, peaks, regions)
- comparison groups and statistical test used
- number of significant hits (DMPs, enriched regions, accessible loci)
- genomic distribution of results (promoter, enhancer, intergenic)
- implicated epigenetic regulators or drug targets
- where output files were saved

## Working Principles

- Methylation beta-values are bounded [0,1]; M-values are better for statistical tests but less intuitive.
- Differential methylation needs biological replicates; single-sample comparisons are descriptive only.
- ChIP-seq peak calls depend heavily on the caller and parameters used; report the source pipeline.
- ATAC-seq accessibility ≠ active transcription; it reflects chromatin openness.
- Epigenetic drug-target implication (e.g., EZH2 from H3K27me3 data) is a hypothesis, not validation.
- Always report the genome assembly version and annotation source.

## Failure Modes

- Mismatched genome assembly between peaks and annotation: coordinates will not map
- Too few replicates for statistical comparison: report descriptive statistics only
- Very large peak files (>1M peaks): subsample or filter by score before analysis
- GTF annotation not provided: peak genomic context cannot be determined
- pysam or pyBigWig not installed: BED file processing still works; bigWig signal extraction does not

## Related Skills

For gene regulatory network inference, activate `grn-tools`.
For variant analysis in genomic regions, activate `variant-analysis-tools`.
For pathway enrichment of epigenetically regulated genes, activate `pathway-enrichment-tools`.
For omics data processing (RNA-seq, single-cell), activate `omics-tools`.
