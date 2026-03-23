#!/usr/bin/env python3
"""Epigenomics analysis: methylation, ChIP-seq peak annotation, and ATAC-seq profiling."""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

try:
    import numpy as np
except Exception:
    np = None

try:
    import pandas as pd
except Exception:
    pd = None


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Epigenomics analysis for drug-target discovery")
    p.add_argument("--input", required=True, help="Input data file")
    p.add_argument("--mode", choices=["methylation", "chipseq", "atacseq"], required=True)

    p.add_argument("--group-column", help="Group/condition column for differential analysis")
    p.add_argument("--group-a", help="Treatment/case group name")
    p.add_argument("--group-b", help="Control group name")
    p.add_argument("--delta-beta-cutoff", type=float, default=0.2)
    p.add_argument("--pvalue-cutoff", type=float, default=0.01)

    p.add_argument("--annotation", help="GTF annotation file for peak context")
    p.add_argument("--promoter-window", type=int, default=2000, help="Promoter window upstream of TSS")

    p.add_argument("--output", default="epigenomics_results.csv")
    p.add_argument("--summary", default="epigenomics_summary.json")
    return p.parse_args()


def require_numpy():
    if np is None:
        raise SystemExit("numpy is required")
    return np


def read_methylation(path: str) -> Any:
    if pd is None:
        raise SystemExit("pandas required for methylation analysis")
    p = Path(path)
    sep = "\t" if p.suffix.lower() in (".tsv", ".tab") else ","
    return pd.read_csv(p, sep=sep)


def differential_methylation(df, group_col: str, group_a: str, group_b: str,
                             delta_cutoff: float, pval_cutoff: float) -> tuple[Any, dict]:
    _np = require_numpy()
    from scipy.stats import mannwhitneyu

    sample_cols_a = []
    sample_cols_b = []
    if group_col in df.columns:
        raise SystemExit("Group column found in rows — expected wide-format beta matrix with samples as columns")

    meta_cols = [c for c in df.columns if not c.replace(".", "").replace("-", "").replace("_", "").replace(" ", "").isdigit()
                 and df[c].dtype == object][:3]
    sample_cols = [c for c in df.columns if c not in meta_cols]

    for c in sample_cols:
        if group_a.lower() in c.lower():
            sample_cols_a.append(c)
        elif group_b.lower() in c.lower():
            sample_cols_b.append(c)

    if not sample_cols_a or not sample_cols_b:
        sample_cols_a = sample_cols[:len(sample_cols) // 2]
        sample_cols_b = sample_cols[len(sample_cols) // 2:]

    results = []
    for idx, row in df.iterrows():
        vals_a = row[sample_cols_a].astype(float).dropna().values
        vals_b = row[sample_cols_b].astype(float).dropna().values
        if len(vals_a) < 2 or len(vals_b) < 2:
            continue
        mean_a = float(_np.mean(vals_a))
        mean_b = float(_np.mean(vals_b))
        delta_beta = mean_a - mean_b
        try:
            _, pval = mannwhitneyu(vals_a, vals_b, alternative="two-sided")
        except Exception:
            pval = 1.0

        sig = abs(delta_beta) >= delta_cutoff and pval <= pval_cutoff
        entry = {col: row[col] for col in meta_cols if col in row}
        entry.update({
            "mean_a": round(mean_a, 4),
            "mean_b": round(mean_b, 4),
            "delta_beta": round(delta_beta, 4),
            "p_value": round(pval, 6),
            "significant": sig,
        })
        results.append(entry)

    sig_count = sum(1 for r in results if r["significant"])
    summary = {
        "total_features": len(results),
        "significant": sig_count,
        "group_a_samples": len(sample_cols_a),
        "group_b_samples": len(sample_cols_b),
        "delta_beta_cutoff": delta_cutoff,
        "pvalue_cutoff": pval_cutoff,
    }
    return results, summary


def parse_bed_peaks(path: str) -> list[dict[str, Any]]:
    peaks = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            peak: dict[str, Any] = {
                "chrom": parts[0],
                "start": int(parts[1]),
                "end": int(parts[2]),
            }
            if len(parts) > 3:
                peak["name"] = parts[3]
            if len(parts) > 4:
                try:
                    peak["score"] = float(parts[4])
                except ValueError:
                    peak["score"] = 0
            if len(parts) > 6:
                try:
                    peak["signal_value"] = float(parts[6])
                except (ValueError, IndexError):
                    pass
            peaks.append(peak)
    return peaks


def parse_gtf_genes(path: str, promoter_window: int) -> list[dict[str, Any]]:
    genes = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attrs = parts[8]
            gene_name = ""
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("gene_name"):
                    gene_name = attr.split('"')[1] if '"' in attr else attr.split(" ")[-1]
                    break

            if strand == "+":
                prom_start = max(0, start - promoter_window)
                prom_end = start
            else:
                prom_start = end
                prom_end = end + promoter_window

            genes.append({
                "chrom": chrom, "start": start, "end": end, "strand": strand,
                "gene_name": gene_name,
                "promoter_start": prom_start, "promoter_end": prom_end,
            })
    return genes


def annotate_peaks(peaks: list[dict], genes: list[dict], promoter_window: int) -> list[dict]:
    gene_lookup: dict[str, list[dict]] = {}
    for g in genes:
        gene_lookup.setdefault(g["chrom"], []).append(g)

    for peak in peaks:
        chrom = peak["chrom"]
        mid = (peak["start"] + peak["end"]) // 2
        peak["genomic_context"] = "intergenic"
        peak["nearest_gene"] = ""
        best_dist = float("inf")

        for g in gene_lookup.get(chrom, []):
            if g["promoter_start"] <= mid <= g["promoter_end"]:
                peak["genomic_context"] = "promoter"
                peak["nearest_gene"] = g["gene_name"]
                best_dist = 0
                break
            elif g["start"] <= mid <= g["end"]:
                peak["genomic_context"] = "gene_body"
                peak["nearest_gene"] = g["gene_name"]
                best_dist = 0
                break
            else:
                dist = min(abs(mid - g["start"]), abs(mid - g["end"]))
                if dist < best_dist:
                    best_dist = dist
                    peak["nearest_gene"] = g["gene_name"]
    return peaks


def main() -> None:
    args = parse_args()

    if args.mode == "methylation":
        print(f"Loading methylation data: {args.input}")
        df = read_methylation(args.input)
        print(f"  {df.shape[0]} features, {df.shape[1]} columns")

        if args.group_a and args.group_b:
            results, meth_summary = differential_methylation(
                df, args.group_column or "", args.group_a, args.group_b,
                args.delta_beta_cutoff, args.pvalue_cutoff)
        else:
            _np = require_numpy()
            sample_cols = [c for c in df.columns if df[c].dtype in (float, "float64", "float32")]
            results = []
            for idx, row in df.iterrows():
                vals = row[sample_cols].astype(float).dropna().values
                entry = {c: row[c] for c in df.columns if c not in sample_cols}
                entry["mean_beta"] = round(float(_np.mean(vals)), 4) if len(vals) > 0 else None
                entry["std_beta"] = round(float(_np.std(vals)), 4) if len(vals) > 0 else None
                results.append(entry)
            meth_summary = {"total_features": len(results), "mode": "descriptive"}

        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        if results:
            keys = list(results[0].keys())
            with open(args.output, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=keys)
                w.writeheader()
                w.writerows(results)

        summary = {"mode": "methylation", **meth_summary, "output_file": args.output}

    elif args.mode in ("chipseq", "atacseq"):
        label = "ChIP-seq" if args.mode == "chipseq" else "ATAC-seq"
        print(f"Loading {label} peaks: {args.input}")
        peaks = parse_bed_peaks(args.input)
        print(f"  {len(peaks)} peaks loaded")

        if args.annotation:
            genes = parse_gtf_genes(args.annotation, args.promoter_window)
            print(f"  {len(genes)} genes from annotation")
            peaks = annotate_peaks(peaks, genes, args.promoter_window)

        context_counts = Counter(p.get("genomic_context", "unknown") for p in peaks)

        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        if peaks:
            keys = list(peaks[0].keys())
            with open(args.output, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=keys)
                w.writeheader()
                w.writerows(peaks)

        promoter_genes = sorted(set(p.get("nearest_gene", "") for p in peaks
                                    if p.get("genomic_context") == "promoter"))

        summary = {
            "mode": args.mode,
            "total_peaks": len(peaks),
            "genomic_context_distribution": dict(context_counts),
            "promoter_genes": promoter_genes[:50],
            "annotation_used": bool(args.annotation),
            "output_file": args.output,
        }
    else:
        raise SystemExit(f"Unknown mode: {args.mode}")

    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary).write_text(json.dumps(summary, indent=2, default=str))
    print(f"Summary: {args.summary}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
