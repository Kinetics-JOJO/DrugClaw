#!/usr/bin/env python3
"""Active learning compound selection: uncertainty, greedy, and balanced strategies."""
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
    p = argparse.ArgumentParser(description="Select compounds for testing via active learning")
    p.add_argument("--labeled", required=True, help="CSV of already-tested compounds with labels")
    p.add_argument("--unlabeled", required=True, help="CSV of candidate compounds to select from")
    p.add_argument("--smiles-column", default="smiles")
    p.add_argument("--id-column", default="compound_id")
    p.add_argument("--target-column", required=True, help="Label column in the labeled set")
    p.add_argument("--acquisition", choices=["uncertainty", "greedy", "balanced"], default="uncertainty")
    p.add_argument("--explore-fraction", type=float, default=0.5,
                   help="Fraction of batch for exploration in balanced mode")
    p.add_argument("--batch-size", type=int, default=50)
    p.add_argument("--feature-backend", choices=["rdkit-morgan", "rdkit-maccs", "ecfp"], default="rdkit-morgan")
    p.add_argument("--n-estimators", type=int, default=100, help="Number of trees in the ensemble")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--output", default="active_learning_selections.csv")
    p.add_argument("--summary", default="active_learning_summary.json")
    return p.parse_args()


def require_deps():
    try:
        import numpy as _np
        from sklearn.ensemble import RandomForestRegressor
        return _np, RandomForestRegressor
    except ImportError as e:
        raise SystemExit(f"numpy and scikit-learn are required: {e}")


def featurize_smiles(smiles_list: list[str], backend: str) -> Any:
    _np, _ = require_deps()
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        raise SystemExit("RDKit required for molecular featurization")

    fps = []
    valid_mask = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            fps.append(None)
            valid_mask.append(False)
            continue
        if backend in ("rdkit-morgan", "ecfp"):
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        else:
            from rdkit.Chem import MACCSkeys
            fp = MACCSkeys.GenMACCSKeys(mol)
        fps.append(list(fp))
        valid_mask.append(True)

    valid_fps = [f for f in fps if f is not None]
    if not valid_fps:
        raise SystemExit("No valid SMILES could be featurized")
    X = _np.array(valid_fps, dtype=float)
    return X, valid_mask


def read_csv(path: str) -> list[dict[str, Any]]:
    p = Path(path)
    delimiter = "\t" if p.suffix.lower() in (".tsv", ".tab") else ","
    with open(p, newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def main() -> None:
    args = parse_args()
    _np, RandomForestRegressor = require_deps()
    _np.random.seed(args.seed)

    print(f"Loading labeled set: {args.labeled}")
    labeled_rows = read_csv(args.labeled)
    print(f"Loading unlabeled set: {args.unlabeled}")
    unlabeled_rows = read_csv(args.unlabeled)

    labeled_smiles = [r.get(args.smiles_column, "") for r in labeled_rows]
    labeled_labels = []
    for r in labeled_rows:
        try:
            labeled_labels.append(float(r[args.target_column]))
        except (ValueError, KeyError):
            labeled_labels.append(float("nan"))

    unlabeled_smiles = [r.get(args.smiles_column, "") for r in unlabeled_rows]

    print(f"Featurizing {len(labeled_smiles)} labeled + {len(unlabeled_smiles)} unlabeled compounds...")
    X_labeled, mask_labeled = featurize_smiles(labeled_smiles, args.feature_backend)
    X_unlabeled, mask_unlabeled = featurize_smiles(unlabeled_smiles, args.feature_backend)

    y_labeled = _np.array([labeled_labels[i] for i in range(len(mask_labeled)) if mask_labeled[i]])
    valid_idx = _np.isfinite(y_labeled)
    X_train = X_labeled[valid_idx]
    y_train = y_labeled[valid_idx]
    print(f"Training on {len(y_train)} valid labeled compounds")

    if len(y_train) < 5:
        print("WARNING: Very few labeled compounds — uncertainty estimates will be unreliable")

    model = RandomForestRegressor(n_estimators=args.n_estimators, random_state=args.seed, n_jobs=-1)
    model.fit(X_train, y_train)

    tree_preds = _np.array([tree.predict(X_unlabeled) for tree in model.estimators_])
    mean_pred = tree_preds.mean(axis=0)
    std_pred = tree_preds.std(axis=0)

    if args.acquisition == "uncertainty":
        scores = std_pred
    elif args.acquisition == "greedy":
        scores = mean_pred
    elif args.acquisition == "balanced":
        n_explore = int(args.batch_size * args.explore_fraction)
        n_exploit = args.batch_size - n_explore
        explore_idx = _np.argsort(-std_pred)[:n_explore]
        exploit_idx = _np.argsort(-mean_pred)
        exploit_idx = _np.array([i for i in exploit_idx if i not in set(explore_idx)])[:n_exploit]
        combined_idx = _np.concatenate([explore_idx, exploit_idx])
        scores = _np.zeros(len(mean_pred))
        for rank, idx in enumerate(combined_idx):
            scores[idx] = len(combined_idx) - rank
    else:
        scores = std_pred

    if args.acquisition != "balanced":
        selected_idx = _np.argsort(-scores)[:args.batch_size]
    else:
        selected_idx = _np.argsort(-scores)[:args.batch_size]

    results = []
    unlabeled_valid_indices = [i for i, m in enumerate(mask_unlabeled) if m]
    for rank, sel in enumerate(selected_idx):
        if sel >= len(unlabeled_valid_indices):
            continue
        orig_idx = unlabeled_valid_indices[sel]
        row = dict(unlabeled_rows[orig_idx])
        row["al_rank"] = rank + 1
        row["predicted_mean"] = round(float(mean_pred[sel]), 6)
        row["predicted_std"] = round(float(std_pred[sel]), 6)
        row["acquisition_score"] = round(float(scores[sel]), 6)
        results.append(row)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    if results:
        keys = list(results[0].keys())
        with open(args.output, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(results)

    summary = {
        "labeled_count": len(y_train),
        "unlabeled_count": len(unlabeled_smiles),
        "unlabeled_valid": sum(mask_unlabeled),
        "acquisition": args.acquisition,
        "batch_size": args.batch_size,
        "selected_count": len(results),
        "feature_backend": args.feature_backend,
        "n_estimators": args.n_estimators,
        "predicted_mean_range": [round(float(mean_pred.min()), 4), round(float(mean_pred.max()), 4)],
        "predicted_std_range": [round(float(std_pred.min()), 4), round(float(std_pred.max()), 4)],
        "output_file": args.output,
    }
    if args.acquisition == "balanced":
        summary["explore_fraction"] = args.explore_fraction

    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary).write_text(json.dumps(summary, indent=2, default=str))
    print(f"Selected {len(results)} compounds")
    print(f"Summary: {args.summary}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
