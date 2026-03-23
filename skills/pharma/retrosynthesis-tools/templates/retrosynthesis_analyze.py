#!/usr/bin/env python3
"""Retrosynthetic analysis: SA scores, disconnection trees, and building-block checks."""
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


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Retrosynthetic analysis and synthetic accessibility scoring")
    p.add_argument("--input", help="CSV/TSV of molecules (sa-score and bb-check modes)")
    p.add_argument("--smiles", nargs="*", help="Inline SMILES (sa-score and disconnection modes)")
    p.add_argument("--smiles-column", default="smiles")
    p.add_argument("--id-column", default="compound_id")
    p.add_argument("--mode", choices=["sa-score", "disconnection", "bb-check"], default="sa-score")
    p.add_argument("--max-depth", type=int, default=3, help="Max disconnection depth")
    p.add_argument("--catalog", help="Building-block catalog CSV for bb-check mode")
    p.add_argument("--catalog-smiles-column", default="smiles")
    p.add_argument("--output", default="retrosynthesis_results.csv")
    p.add_argument("--summary", default="retrosynthesis_summary.json")
    return p.parse_args()


def require_rdkit():
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors
        return Chem, Descriptors, rdMolDescriptors
    except ImportError:
        raise SystemExit("RDKit is required for retrosynthetic analysis")


def compute_sa_score(mol) -> float:
    """Heuristic synthetic accessibility score (1=easy, 10=hard).

    Approximation inspired by Ertl & Schuffenhauer (J. Cheminform. 2009)
    using ring complexity, stereo centers, heavy atom count, and sp3 fraction.
    """
    Chem, Descriptors, rdMolDescriptors = require_rdkit()
    if mol is None:
        return float("nan")

    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    num_heavy = mol.GetNumHeavyAtoms()
    num_stereo = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
    num_rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
    fraction_sp3 = rdMolDescriptors.CalcFractionCSP3(mol)

    ring_penalty = min(num_rings * 0.5, 3.0)
    size_term = math.log10(max(num_heavy, 1)) * 1.5
    stereo_penalty = num_stereo * 0.4
    flexibility = min(num_rot * 0.1, 1.5)
    sp3_bonus = (1.0 - fraction_sp3) * 0.5

    score = 1.0 + ring_penalty + size_term + stereo_penalty + flexibility + sp3_bonus
    return min(max(round(score, 2), 1.0), 10.0)


def sa_score_batch(smiles_list: list[str], id_list: list[str]) -> list[dict[str, Any]]:
    Chem, _, _ = require_rdkit()
    results = []
    for smi, cid in zip(smiles_list, id_list):
        mol = Chem.MolFromSmiles(smi)
        sa = compute_sa_score(mol)
        results.append({
            "id": cid,
            "smiles": smi,
            "sa_score": sa,
            "valid": mol is not None,
        })
    return results


def disconnection_analysis(smiles: str, max_depth: int) -> dict[str, Any]:
    Chem, _, _ = require_rdkit()
    from rdkit.Chem import BRICS

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"smiles": smiles, "error": "invalid SMILES", "fragments": []}

    try:
        fragments = list(BRICS.BRICSDecompose(mol, minFragmentSize=3))
    except Exception as e:
        return {"smiles": smiles, "error": str(e), "fragments": []}

    sa = compute_sa_score(mol)
    return {
        "smiles": smiles,
        "sa_score": sa,
        "num_fragments": len(fragments),
        "fragments": sorted(fragments),
        "max_depth_requested": max_depth,
    }


def bb_check(routes_path: str, catalog_path: str, catalog_smi_col: str) -> list[dict[str, Any]]:
    Chem, _, _ = require_rdkit()

    catalog_smiles = set()
    with open(catalog_path, newline="") as f:
        for row in csv.DictReader(f):
            smi = row.get(catalog_smi_col, "").strip()
            if smi:
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    catalog_smiles.add(Chem.MolToSmiles(mol))

    routes = json.loads(Path(routes_path).read_text())
    if isinstance(routes, dict):
        routes = [routes]

    results = []
    for route in routes:
        frags = route.get("fragments", [])
        available = []
        missing = []
        for frag in frags:
            clean = frag.replace("[1*]", "[H]").replace("[2*]", "[H]").replace("[3*]", "[H]")
            mol = Chem.MolFromSmiles(clean)
            canon = Chem.MolToSmiles(mol) if mol else clean
            if canon in catalog_smiles:
                available.append(frag)
            else:
                missing.append(frag)
        results.append({
            "smiles": route.get("smiles", ""),
            "total_fragments": len(frags),
            "available": len(available),
            "missing": len(missing),
            "coverage_pct": round(len(available) / max(len(frags), 1) * 100, 1),
            "missing_fragments": missing,
        })
    return results


def read_smiles(args) -> tuple[list[str], list[str]]:
    if args.smiles:
        return args.smiles, [f"mol_{i}" for i in range(len(args.smiles))]
    if args.input:
        p = Path(args.input)
        delimiter = "\t" if p.suffix.lower() in (".tsv", ".tab") else ","
        with open(p, newline="") as f:
            rows = list(csv.DictReader(f, delimiter=delimiter))
        smiles = [r.get(args.smiles_column, "") for r in rows]
        ids = [r.get(args.id_column, f"mol_{i}") for i, r in enumerate(rows)]
        return smiles, ids
    raise SystemExit("Provide --input or --smiles")


def main() -> None:
    args = parse_args()

    if args.mode == "sa-score":
        smiles_list, id_list = read_smiles(args)
        print(f"Computing SA scores for {len(smiles_list)} molecules...")
        results = sa_score_batch(smiles_list, id_list)
        valid = [r for r in results if r["valid"]]
        scores = [r["sa_score"] for r in valid if not math.isnan(r["sa_score"])]

        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["id", "smiles", "sa_score", "valid"])
            w.writeheader()
            w.writerows(results)

        summary = {
            "mode": "sa-score",
            "total": len(results),
            "valid": len(valid),
            "invalid": len(results) - len(valid),
            "sa_score_mean": round(sum(scores) / max(len(scores), 1), 2),
            "sa_score_min": round(min(scores), 2) if scores else None,
            "sa_score_max": round(max(scores), 2) if scores else None,
            "output_file": args.output,
        }

    elif args.mode == "disconnection":
        smiles_list, _ = read_smiles(args)
        print(f"Running BRICS disconnection for {len(smiles_list)} molecules...")
        all_routes = [disconnection_analysis(smi, args.max_depth) for smi in smiles_list]

        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        if args.output.endswith(".json"):
            Path(args.output).write_text(json.dumps(all_routes, indent=2))
        else:
            flat = []
            for r in all_routes:
                flat.append({
                    "smiles": r["smiles"],
                    "sa_score": r.get("sa_score", ""),
                    "num_fragments": r.get("num_fragments", 0),
                    "fragments": "; ".join(r.get("fragments", [])),
                })
            with open(args.output, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(flat[0].keys()))
                w.writeheader()
                w.writerows(flat)

        summary = {
            "mode": "disconnection",
            "total_molecules": len(all_routes),
            "routes": all_routes,
            "output_file": args.output,
        }

    elif args.mode == "bb-check":
        if not args.input or not args.catalog:
            raise SystemExit("--input (routes JSON) and --catalog (CSV) required for bb-check")
        results = bb_check(args.input, args.catalog, args.catalog_smiles_column)

        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w", newline="") as f:
            keys = ["smiles", "total_fragments", "available", "missing", "coverage_pct", "missing_fragments"]
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for r in results:
                row = dict(r)
                row["missing_fragments"] = "; ".join(row["missing_fragments"])
                w.writerow(row)

        summary = {
            "mode": "bb-check",
            "total_routes": len(results),
            "results": results,
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
