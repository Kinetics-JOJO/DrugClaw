#!/usr/bin/env python3
"""Pathway and gene-set enrichment analysis: ORA and GSEA."""
from __future__ import annotations

import argparse
import csv
import json
import math
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
    p = argparse.ArgumentParser(description="Pathway enrichment analysis (ORA / GSEA)")
    p.add_argument("--input", required=True, help="Gene list (txt, one per line) or ranked gene table (CSV)")
    p.add_argument("--gene-column", help="Gene column in CSV (for GSEA)")
    p.add_argument("--rank-column", help="Rank metric column (for GSEA, e.g. log2fc)")
    p.add_argument("--mode", choices=["ora", "gsea"], default="ora")
    p.add_argument("--organism", default="Homo sapiens")
    p.add_argument("--gene-sets", nargs="+", default=["KEGG_2021_Human"],
                   help="Gene-set library names for Enrichr/gseapy")
    p.add_argument("--background-size", type=int, default=20000, help="Background gene count for ORA")
    p.add_argument("--fdr-cutoff", type=float, default=0.05)
    p.add_argument("--permutations", type=int, default=1000, help="Permutations for GSEA")
    p.add_argument("--output", default="enrichment_results.csv")
    p.add_argument("--summary", default="enrichment_summary.json")
    return p.parse_args()


def read_gene_list(path: str) -> list[str]:
    lines = Path(path).read_text().strip().split("\n")
    return [g.strip() for g in lines if g.strip() and not g.startswith("#")]


def read_ranked_genes(path: str, gene_col: str, rank_col: str) -> list[tuple[str, float]]:
    delimiter = "\t" if path.endswith((".tsv", ".tab")) else ","
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter=delimiter))
    genes = []
    for r in rows:
        g = r.get(gene_col, "").strip()
        try:
            score = float(r[rank_col])
        except (ValueError, KeyError):
            continue
        if g:
            genes.append((g, score))
    genes.sort(key=lambda x: x[1], reverse=True)
    return genes


def try_gseapy_ora(gene_list: list[str], gene_sets: list[str], organism: str) -> list[dict[str, Any]] | None:
    try:
        import gseapy as gp
    except ImportError:
        return None

    results = []
    for gs in gene_sets:
        try:
            enr = gp.enrich(gene_list=gene_list, gene_sets=gs, organism=organism,
                            outdir=None, no_plot=True, cutoff=1.0)
            if enr.results is not None and not enr.results.empty:
                for _, row in enr.results.iterrows():
                    results.append({
                        "term": row.get("Term", ""),
                        "gene_set": gs,
                        "overlap": row.get("Overlap", ""),
                        "p_value": row.get("P-value", 1.0),
                        "adjusted_p_value": row.get("Adjusted P-value", 1.0),
                        "odds_ratio": row.get("Odds Ratio", ""),
                        "combined_score": row.get("Combined Score", ""),
                        "genes": row.get("Genes", ""),
                    })
        except Exception:
            continue
    return results if results else None


def try_gseapy_gsea(ranked_genes: list[tuple[str, float]], gene_sets: list[str],
                    permutations: int) -> list[dict[str, Any]] | None:
    try:
        import gseapy as gp
        import pandas as _pd
    except ImportError:
        return None

    rnk = _pd.DataFrame(ranked_genes, columns=["gene", "score"])
    results = []
    for gs in gene_sets:
        try:
            pre_res = gp.prerank(rnk=rnk, gene_sets=gs, permutation_num=permutations,
                                 outdir=None, no_plot=True, seed=42)
            if pre_res.res2d is not None and not pre_res.res2d.empty:
                for _, row in pre_res.res2d.iterrows():
                    results.append({
                        "term": row.get("Term", ""),
                        "gene_set": gs,
                        "es": row.get("ES", ""),
                        "nes": row.get("NES", ""),
                        "p_value": row.get("NOM p-val", 1.0),
                        "fdr": row.get("FDR q-val", 1.0),
                        "size": row.get("Tag %", ""),
                        "genes": row.get("Lead_genes", ""),
                    })
        except Exception:
            continue
    return results if results else None


def fisher_ora_fallback(gene_list: list[str], background_size: int) -> list[dict[str, Any]]:
    """Minimal ORA using scipy Fisher exact test with a stub gene-set."""
    try:
        from scipy.stats import fisher_exact
    except ImportError:
        return [{"note": "scipy not available; install gseapy or scipy for enrichment analysis"}]

    placeholder_sets = {
        "Apoptosis": {"TP53", "BCL2", "BAX", "CASP3", "CASP9", "CYCS", "APAF1", "BID", "XIAP", "MCL1"},
        "Cell_Cycle": {"CDK1", "CDK2", "CCND1", "CCNE1", "RB1", "E2F1", "TP53", "CDKN1A", "CDKN2A", "MYC"},
        "PI3K_AKT_Signaling": {"PIK3CA", "AKT1", "MTOR", "PTEN", "TSC1", "TSC2", "GSK3B", "FOXO3", "BAD", "MDM2"},
        "MAPK_Signaling": {"KRAS", "BRAF", "MAP2K1", "MAPK1", "MAPK3", "RAF1", "HRAS", "NRAS", "ERK1", "ERK2"},
        "DNA_Repair": {"BRCA1", "BRCA2", "ATM", "ATR", "CHEK1", "CHEK2", "RAD51", "PARP1", "MLH1", "MSH2"},
    }

    gene_set_upper = set(g.upper() for g in gene_list)
    n_list = len(gene_list)
    results = []

    for name, genes in placeholder_sets.items():
        genes_upper = set(g.upper() for g in genes)
        overlap = gene_set_upper & genes_upper
        a = len(overlap)
        b = len(genes_upper) - a
        c = n_list - a
        d = background_size - a - b - c
        if d < 0:
            d = 0
        _, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
        results.append({
            "term": name,
            "gene_set": "stub_pathways",
            "overlap": f"{a}/{len(genes_upper)}",
            "p_value": round(pval, 6),
            "adjusted_p_value": round(min(pval * len(placeholder_sets), 1.0), 6),
            "genes": ",".join(sorted(overlap)),
            "note": "stub gene sets — install gseapy for real pathway databases",
        })
    results.sort(key=lambda x: x["p_value"])
    return results


def main() -> None:
    args = parse_args()

    if args.mode == "ora":
        gene_list = read_gene_list(args.input)
        print(f"ORA: {len(gene_list)} genes, gene sets: {args.gene_sets}")
        results = try_gseapy_ora(gene_list, args.gene_sets, args.organism)
        if results is None:
            print("gseapy not available, falling back to Fisher exact test with stub sets")
            results = fisher_ora_fallback(gene_list, args.background_size)
    else:
        if not args.gene_column or not args.rank_column:
            raise SystemExit("--gene-column and --rank-column required for GSEA mode")
        ranked = read_ranked_genes(args.input, args.gene_column, args.rank_column)
        print(f"GSEA: {len(ranked)} ranked genes, gene sets: {args.gene_sets}")
        results = try_gseapy_gsea(ranked, args.gene_sets, args.permutations)
        if results is None:
            raise SystemExit("gseapy required for GSEA mode. Install with: pip install gseapy")

    significant = [r for r in results if float(r.get("adjusted_p_value", r.get("fdr", 1.0))) <= args.fdr_cutoff]
    print(f"Total results: {len(results)}, significant (FDR ≤ {args.fdr_cutoff}): {len(significant)}")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    if results:
        keys = list(results[0].keys())
        with open(args.output, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(results)

    summary = {
        "mode": args.mode,
        "organism": args.organism,
        "gene_sets": args.gene_sets,
        "total_results": len(results),
        "significant_count": len(significant),
        "fdr_cutoff": args.fdr_cutoff,
        "top_terms": [{"term": r.get("term", ""), "p": r.get("p_value", "")} for r in results[:10]],
        "output_file": args.output,
    }

    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary).write_text(json.dumps(summary, indent=2, default=str))
    print(f"Summary: {args.summary}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
