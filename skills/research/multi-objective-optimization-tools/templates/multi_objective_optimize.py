#!/usr/bin/env python3
"""Multi-objective optimization: Pareto-front extraction and weighted scalarization."""
from __future__ import annotations

import argparse
import csv
import json
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
    p = argparse.ArgumentParser(description="Multi-objective Pareto and scalarization ranking")
    p.add_argument("--input", required=True, help="CSV/TSV with objective columns")
    p.add_argument("--id-column", help="Row identifier column")
    p.add_argument("--objective", action="append", required=True,
                   help="Objective spec: column_name:maximize or column_name:minimize. Repeat for each objective.")
    p.add_argument("--weight", action="append", default=[],
                   help="Scalarization weight: column_name:weight. Repeat for each. If omitted, equal weights.")
    p.add_argument("--normalize", choices=["minmax", "zscore", "rank"], default="minmax")
    p.add_argument("--output", default="moo_results.csv")
    p.add_argument("--summary", default="moo_summary.json")
    return p.parse_args()


def parse_objectives(specs: list[str]) -> list[tuple[str, str]]:
    objectives = []
    for s in specs:
        parts = s.rsplit(":", 1)
        if len(parts) != 2 or parts[1] not in ("maximize", "minimize"):
            raise SystemExit(f"Invalid objective spec: {s}. Use column_name:maximize or column_name:minimize")
        objectives.append((parts[0], parts[1]))
    return objectives


def parse_weights(specs: list[str], objectives: list[tuple[str, str]]) -> dict[str, float]:
    if not specs:
        return {name: 1.0 / len(objectives) for name, _ in objectives}
    weights = {}
    for s in specs:
        parts = s.rsplit(":", 1)
        if len(parts) != 2:
            raise SystemExit(f"Invalid weight spec: {s}")
        weights[parts[0]] = float(parts[1])
    total = sum(weights.values())
    if total > 0:
        weights = {k: v / total for k, v in weights.items()}
    for name, _ in objectives:
        if name not in weights:
            weights[name] = 0.0
    return weights


def read_data(path: str) -> list[dict[str, Any]]:
    p = Path(path)
    delimiter = "\t" if p.suffix.lower() in (".tsv", ".tab") else ","
    with open(p, newline="") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        return list(reader)


def normalize_column(values: list[float], method: str) -> list[float]:
    if np is None:
        raise SystemExit("numpy required for normalization")
    arr = np.array(values, dtype=float)
    if method == "minmax":
        lo, hi = np.nanmin(arr), np.nanmax(arr)
        if hi - lo < 1e-12:
            return [0.5] * len(values)
        return ((arr - lo) / (hi - lo)).tolist()
    elif method == "zscore":
        m, s = np.nanmean(arr), np.nanstd(arr)
        if s < 1e-12:
            return [0.0] * len(values)
        return ((arr - m) / s).tolist()
    elif method == "rank":
        from scipy.stats import rankdata
        ranks = rankdata(arr, method="average", nan_policy="omit")
        lo, hi = ranks.min(), ranks.max()
        if hi - lo < 1e-12:
            return [0.5] * len(values)
        return ((ranks - lo) / (hi - lo)).tolist()
    return values


def dominates(a: list[float], b: list[float]) -> bool:
    dominated = False
    for ai, bi in zip(a, b):
        if ai < bi:
            return False
        if ai > bi:
            dominated = True
    return dominated


def pareto_rank(rows: list[dict[str, Any]], objectives: list[tuple[str, str]],
                normalize: str) -> list[dict[str, Any]]:
    if np is None:
        raise SystemExit("numpy is required")

    obj_values: dict[str, list[float]] = {}
    for col, _ in objectives:
        vals = []
        for r in rows:
            try:
                vals.append(float(r[col]))
            except (ValueError, KeyError):
                vals.append(float("nan"))
        obj_values[col] = vals

    norm_values: dict[str, list[float]] = {}
    for col, direction in objectives:
        normed = normalize_column(obj_values[col], normalize)
        if direction == "minimize":
            normed = [1.0 - v for v in normed]
        norm_values[col] = normed

    n = len(rows)
    vectors = []
    for i in range(n):
        vectors.append([norm_values[col][i] for col, _ in objectives])

    ranks = [0] * n
    assigned = [False] * n
    current_rank = 1
    remaining = n

    while remaining > 0:
        front = []
        for i in range(n):
            if assigned[i]:
                continue
            is_dominated = False
            for j in range(n):
                if assigned[j] or i == j:
                    continue
                if dominates(vectors[j], vectors[i]):
                    is_dominated = True
                    break
            if not is_dominated:
                front.append(i)
        for idx in front:
            ranks[idx] = current_rank
            assigned[idx] = True
            remaining -= 1
        current_rank += 1

    for i, r in enumerate(rows):
        r["pareto_rank"] = ranks[i]
        for col, _ in objectives:
            r[f"{col}_normalized"] = round(norm_values[col][i], 6)

    rows.sort(key=lambda x: x["pareto_rank"])
    return rows


def scalarize(rows: list[dict[str, Any]], objectives: list[tuple[str, str]],
              weights: dict[str, float]) -> list[dict[str, Any]]:
    for r in rows:
        score = 0.0
        for col, _ in objectives:
            nkey = f"{col}_normalized"
            if nkey in r:
                try:
                    score += weights.get(col, 0.0) * float(r[nkey])
                except (ValueError, TypeError):
                    pass
        r["scalarized_score"] = round(score, 6)
    rows.sort(key=lambda x: x.get("scalarized_score", 0.0), reverse=True)
    return rows


def main() -> None:
    args = parse_args()
    objectives = parse_objectives(args.objective)
    weights = parse_weights(args.weight, objectives)

    print(f"Reading {args.input}")
    rows = read_data(args.input)
    print(f"  {len(rows)} rows, {len(objectives)} objectives")

    rows = pareto_rank(rows, objectives, args.normalize)
    rows = scalarize(rows, objectives, weights)

    rank1 = [r for r in rows if r.get("pareto_rank") == 1]
    print(f"  Pareto rank-1 solutions: {len(rank1)}")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    if rows:
        keys = list(rows[0].keys())
        with open(args.output, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)

    id_col = args.id_column
    summary = {
        "objectives": [{"column": c, "direction": d} for c, d in objectives],
        "normalization": args.normalize,
        "weights": weights,
        "total_rows": len(rows),
        "pareto_rank1_count": len(rank1),
        "pareto_rank1_ids": [r.get(id_col, i) for i, r in enumerate(rank1)] if id_col else list(range(len(rank1))),
        "output_file": args.output,
    }
    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary).write_text(json.dumps(summary, indent=2, default=str))
    print(f"Summary: {args.summary}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
