#!/usr/bin/env python3
"""Network pharmacology: drug-target-pathway-disease network construction and analysis."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

try:
    import networkx as nx
except Exception:
    nx = None

try:
    import requests
except Exception:
    requests = None


CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
STITCH_BASE = "http://stitch.embl.de/api"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Network pharmacology analysis")
    p.add_argument("--drug", help="Single drug name or identifier")
    p.add_argument("--drugs", help="Comma-separated drug names for comparison/combination")
    p.add_argument("--source", choices=["chembl", "stitch", "local"], default="chembl")
    p.add_argument("--mode", choices=["single", "compare", "combination"], default="single")
    p.add_argument("--input", help="Local target table CSV for source=local")
    p.add_argument("--drug-column", default="drug")
    p.add_argument("--target-column", default="target")
    p.add_argument("--activity-column", help="Activity value column for filtering")
    p.add_argument("--activity-cutoff", type=float, default=1000, help="Activity cutoff in nM")
    p.add_argument("--pathway-db", choices=["KEGG", "Reactome", "none"], default="none")
    p.add_argument("--overlap-analysis", action="store_true")
    p.add_argument("--output", default="network_pharmacology.graphml")
    p.add_argument("--table", default="network_targets.csv")
    p.add_argument("--summary", default="network_summary.json")
    return p.parse_args()


def search_chembl_targets(drug_name: str, cutoff_nm: float) -> list[dict[str, Any]]:
    if requests is None:
        raise SystemExit("requests required for ChEMBL queries")

    search_url = f"{CHEMBL_BASE}/molecule/search.json"
    resp = requests.get(search_url, params={"q": drug_name, "limit": 5}, timeout=20)
    if resp.status_code != 200:
        return []

    molecules = resp.json().get("molecules", [])
    if not molecules:
        return []

    chembl_id = molecules[0].get("molecule_chembl_id", "")
    if not chembl_id:
        return []

    activity_url = f"{CHEMBL_BASE}/activity.json"
    params = {
        "molecule_chembl_id": chembl_id,
        "pchembl_value__isnull": "false",
        "limit": 100,
    }
    resp = requests.get(activity_url, params=params, timeout=30)
    if resp.status_code != 200:
        return []

    activities = resp.json().get("activities", [])
    targets = []
    seen = set()
    for act in activities:
        target_id = act.get("target_chembl_id", "")
        pchembl = act.get("pchembl_value")
        if not target_id or target_id in seen:
            continue
        if pchembl:
            try:
                nm_val = 10 ** (9 - float(pchembl))
                if nm_val > cutoff_nm:
                    continue
            except (ValueError, TypeError):
                pass
        seen.add(target_id)
        targets.append({
            "drug": drug_name,
            "drug_chembl_id": chembl_id,
            "target_chembl_id": target_id,
            "target_name": act.get("target_pref_name", ""),
            "target_organism": act.get("target_organism", ""),
            "pchembl_value": pchembl,
            "activity_type": act.get("standard_type", ""),
            "source": "chembl",
        })
    return targets


def read_local_targets(path: str, drug_col: str, target_col: str,
                       activity_col: str | None, cutoff: float) -> list[dict[str, Any]]:
    p = Path(path)
    delimiter = "\t" if p.suffix.lower() in (".tsv", ".tab") else ","
    with open(p, newline="") as f:
        rows = list(csv.DictReader(f, delimiter=delimiter))

    targets = []
    for r in rows:
        drug = r.get(drug_col, "")
        target = r.get(target_col, "")
        if not drug or not target:
            continue
        if activity_col and activity_col in r:
            try:
                val = float(r[activity_col])
                if val > cutoff:
                    continue
            except ValueError:
                pass
        targets.append({
            "drug": drug,
            "target_name": target,
            "source": "local",
            **{k: v for k, v in r.items() if k not in (drug_col, target_col)},
        })
    return targets


def build_network(targets: list[dict[str, Any]]) -> Any:
    if nx is None:
        raise SystemExit("networkx required")
    G = nx.Graph()
    for t in targets:
        drug = t.get("drug", "")
        target = t.get("target_name", "") or t.get("target_chembl_id", "")
        if drug and target:
            G.add_node(drug, entity_type="drug")
            G.add_node(target, entity_type="target")
            G.add_edge(drug, target, relation="targets",
                       pchembl=t.get("pchembl_value", ""),
                       source=t.get("source", ""))
    return G


def topology_metrics(G: Any) -> dict[str, Any]:
    degree = dict(G.degree())
    try:
        betweenness = nx.betweenness_centrality(G)
    except Exception:
        betweenness = {n: 0 for n in G.nodes()}

    target_nodes = [n for n, d in G.nodes(data=True) if d.get("entity_type") == "target"]
    hub_targets = sorted(target_nodes, key=lambda n: degree.get(n, 0), reverse=True)[:20]

    return {
        "nodes": G.number_of_nodes(),
        "edges": G.number_of_edges(),
        "connected_components": nx.number_connected_components(G),
        "hub_targets": [{"target": n, "degree": degree[n],
                         "betweenness": round(betweenness.get(n, 0), 4)} for n in hub_targets],
    }


def overlap_analysis(targets: list[dict[str, Any]]) -> dict[str, Any]:
    drug_targets: dict[str, set[str]] = {}
    for t in targets:
        drug = t.get("drug", "")
        target = t.get("target_name", "") or t.get("target_chembl_id", "")
        if drug and target:
            drug_targets.setdefault(drug, set()).add(target)

    drugs = list(drug_targets.keys())
    overlap = {}
    for i in range(len(drugs)):
        for j in range(i + 1, len(drugs)):
            shared = drug_targets[drugs[i]] & drug_targets[drugs[j]]
            overlap[f"{drugs[i]}__vs__{drugs[j]}"] = {
                "shared_count": len(shared),
                "shared_targets": sorted(shared),
                "unique_to_first": len(drug_targets[drugs[i]] - drug_targets[drugs[j]]),
                "unique_to_second": len(drug_targets[drugs[j]] - drug_targets[drugs[i]]),
            }
    return {"drug_target_counts": {d: len(t) for d, t in drug_targets.items()}, "overlaps": overlap}


def main() -> None:
    args = parse_args()

    drug_list = []
    if args.drug:
        drug_list = [args.drug]
    elif args.drugs:
        drug_list = [d.strip() for d in args.drugs.split(",")]

    if not drug_list and not args.input:
        raise SystemExit("Provide --drug, --drugs, or --input")

    all_targets: list[dict[str, Any]] = []
    if args.source == "local" and args.input:
        all_targets = read_local_targets(args.input, args.drug_column, args.target_column,
                                         args.activity_column, args.activity_cutoff)
    else:
        for drug in drug_list:
            print(f"Querying {args.source} targets for: {drug}")
            if args.source == "chembl":
                targets = search_chembl_targets(drug, args.activity_cutoff)
            else:
                targets = search_chembl_targets(drug, args.activity_cutoff)
            print(f"  Found {len(targets)} targets")
            all_targets.extend(targets)

    if not all_targets:
        print("No targets found.")

    G = build_network(all_targets)
    metrics = topology_metrics(G)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    if nx is not None:
        nx.write_graphml(G, args.output)
    print(f"Network: {metrics['nodes']} nodes, {metrics['edges']} edges")

    Path(args.table).parent.mkdir(parents=True, exist_ok=True)
    if all_targets:
        keys = list(all_targets[0].keys())
        with open(args.table, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(all_targets)

    summary: dict[str, Any] = {
        "drugs": drug_list,
        "source": args.source,
        "activity_cutoff_nm": args.activity_cutoff,
        "total_targets": len(all_targets),
        "network": metrics,
        "output_file": args.output,
        "table_file": args.table,
    }

    if args.overlap_analysis or args.mode in ("compare", "combination"):
        ol = overlap_analysis(all_targets)
        summary["overlap"] = ol

    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary).write_text(json.dumps(summary, indent=2, default=str))
    print(f"Summary: {args.summary}")


if __name__ == "__main__":
    main()
