#!/usr/bin/env python3
"""De novo molecular generation: R-group enumeration, MMP transforms, and diversity filtering."""
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


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate and filter molecular libraries")
    p.add_argument("--mode", choices=["rgroup", "mmp", "diversify"], required=True)
    p.add_argument("--input", help="CSV/TSV of seed molecules (mmp/diversify modes)")
    p.add_argument("--smiles-column", default="smiles")
    p.add_argument("--id-column", default="compound_id")

    p.add_argument("--core", help="SMARTS core for R-group enumeration")
    p.add_argument("--rgroups", help="JSON array of R-group lists, e.g. [[\"F\",\"Cl\"],[\"OC\",\"NC\"]]")

    p.add_argument("--diversity-cutoff", type=float, default=0.4, help="Tanimoto cutoff for diversity filtering")
    p.add_argument("--max-analogs", type=int, default=500)
    p.add_argument("--apply-lipinski", action="store_true", help="Filter by Lipinski rule of five")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--output", default="generated_library.csv")
    p.add_argument("--summary", default="generated_summary.json")
    return p.parse_args()


def require_rdkit():
    try:
        from rdkit import Chem, DataStructs
        from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
        return Chem, AllChem, Descriptors, DataStructs, rdMolDescriptors
    except ImportError:
        raise SystemExit("RDKit is required for molecular generation")


def lipinski_filter(mol) -> bool:
    _, _, Descriptors, _, _ = require_rdkit()
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    return mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10


def compute_properties(mol) -> dict[str, Any]:
    _, _, Descriptors, _, rdMolDescriptors = require_rdkit()
    return {
        "mw": round(Descriptors.MolWt(mol), 2),
        "logp": round(Descriptors.MolLogP(mol), 2),
        "hbd": Descriptors.NumHDonors(mol),
        "hba": Descriptors.NumHAcceptors(mol),
        "tpsa": round(Descriptors.TPSA(mol), 2),
        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
    }


def rgroup_enumerate(core_smarts: str, rgroup_lists: list[list[str]]) -> list[str]:
    Chem, _, _, _, _ = require_rdkit()

    from itertools import product
    combos = list(product(*rgroup_lists))
    results = []
    for combo in combos:
        smi = core_smarts
        for i, rg in enumerate(combo):
            smi = smi.replace(f"[*:{i + 1}]", rg)
        mol = Chem.MolFromSmiles(smi)
        if mol:
            results.append(Chem.MolToSmiles(mol))
    return list(set(results))


def mmp_transforms(smiles_list: list[str]) -> list[str]:
    Chem, _, _, _, _ = require_rdkit()
    from rdkit.Chem import BRICS

    generated = set()
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        try:
            fragments = list(BRICS.BRICSDecompose(mol, minFragmentSize=3))
            recombined = list(BRICS.BRICSBuild([Chem.MolFromSmiles(f) for f in fragments if Chem.MolFromSmiles(f)]))
            for new_mol in recombined[:50]:
                try:
                    new_smi = Chem.MolToSmiles(new_mol)
                    if new_smi and new_smi not in smiles_list:
                        generated.add(new_smi)
                except Exception:
                    continue
        except Exception:
            continue
    return list(generated)


def diversity_filter(smiles_list: list[str], cutoff: float, max_keep: int) -> list[str]:
    Chem, AllChem, _, DataStructs, _ = require_rdkit()

    mols_fps = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            mols_fps.append((smi, fp))

    if not mols_fps:
        return []

    selected = [mols_fps[0]]
    for smi, fp in mols_fps[1:]:
        if len(selected) >= max_keep:
            break
        max_sim = max(DataStructs.TanimotoSimilarity(fp, sel_fp) for _, sel_fp in selected)
        if max_sim <= (1.0 - cutoff):
            selected.append((smi, fp))

    return [smi for smi, _ in selected]


def read_smiles(path: str, smi_col: str, id_col: str) -> tuple[list[str], list[str]]:
    p = Path(path)
    delimiter = "\t" if p.suffix.lower() in (".tsv", ".tab") else ","
    with open(p, newline="") as f:
        rows = list(csv.DictReader(f, delimiter=delimiter))
    return ([r.get(smi_col, "") for r in rows],
            [r.get(id_col, f"mol_{i}") for i, r in enumerate(rows)])


def main() -> None:
    args = parse_args()
    Chem, AllChem, Descriptors, DataStructs, _ = require_rdkit()

    generated_smiles: list[str] = []
    seed_count = 0

    if args.mode == "rgroup":
        if not args.core or not args.rgroups:
            raise SystemExit("--core and --rgroups required for rgroup mode")
        rgroups = json.loads(args.rgroups)
        print(f"Enumerating R-groups on core: {args.core}")
        generated_smiles = rgroup_enumerate(args.core, rgroups)
        seed_count = 1

    elif args.mode == "mmp":
        if not args.input:
            raise SystemExit("--input required for mmp mode")
        smiles_list, _ = read_smiles(args.input, args.smiles_column, args.id_column)
        seed_count = len(smiles_list)
        print(f"Generating MMP analogs from {seed_count} seeds...")
        generated_smiles = mmp_transforms(smiles_list)

    elif args.mode == "diversify":
        if not args.input:
            raise SystemExit("--input required for diversify mode")
        smiles_list, _ = read_smiles(args.input, args.smiles_column, args.id_column)
        seed_count = len(smiles_list)
        print(f"Generating diverse analogs from {seed_count} seeds...")
        mmp_results = mmp_transforms(smiles_list)
        all_candidates = list(set(smiles_list + mmp_results))
        generated_smiles = diversity_filter(all_candidates, args.diversity_cutoff, args.max_analogs)

    raw_count = len(generated_smiles)
    print(f"Raw candidates: {raw_count}")

    if args.apply_lipinski:
        filtered = []
        for smi in generated_smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol and lipinski_filter(mol):
                filtered.append(smi)
        generated_smiles = filtered
        print(f"After Lipinski filter: {len(generated_smiles)}")

    results = []
    for i, smi in enumerate(generated_smiles):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        props = compute_properties(mol)
        results.append({
            "id": f"gen_{i:05d}",
            "smiles": smi,
            **props,
        })

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    if results:
        keys = list(results[0].keys())
        with open(args.output, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(results)

    summary = {
        "mode": args.mode,
        "seed_count": seed_count,
        "raw_generated": raw_count,
        "after_filters": len(results),
        "lipinski_applied": args.apply_lipinski,
        "output_file": args.output,
    }
    if results:
        mws = [r["mw"] for r in results]
        logps = [r["logp"] for r in results]
        summary["mw_range"] = [min(mws), max(mws)]
        summary["logp_range"] = [min(logps), max(logps)]

    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary).write_text(json.dumps(summary, indent=2, default=str))
    print(f"Generated {len(results)} compounds")
    print(f"Summary: {args.summary}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
