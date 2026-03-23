#!/usr/bin/env python3
"""Systematic drug repurposing: target overlap, genetic association, and clinical evidence triangulation."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

try:
    import requests
except Exception:
    requests = None

OPENTARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"
CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
CTGOV_API = "https://clinicaltrials.gov/api/v2/studies"

DISEASE_TARGETS_QUERY = """
query DiseaseTargets($efoId: String!, $size: Int!) {
  disease(efoId: $efoId) {
    id
    name
    associatedTargets(page: {index: 0, size: $size}) {
      count
      rows { target { id approvedSymbol } score }
    }
  }
}
"""

DISEASE_SEARCH_QUERY = """
query SearchDisease($q: String!, $size: Int!) {
  search(queryString: $q, entityNames: ["disease"], page: {index: 0, size: $size}) {
    hits { id name description entity }
  }
}
"""

TARGET_SEARCH_QUERY = """
query SearchTarget($q: String!, $size: Int!) {
  search(queryString: $q, entityNames: ["target"], page: {index: 0, size: $size}) {
    hits { id name description entity }
  }
}
"""

TARGET_DISEASES_QUERY = """
query TargetDiseases($ensemblId: String!, $size: Int!) {
  target(ensemblId: $ensemblId) {
    id approvedSymbol approvedName
    associatedDiseases(page: {index: 0, size: $size}) {
      count
      rows { disease { id name } score }
    }
  }
}
"""


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Drug repurposing via target overlap and clinical evidence")
    p.add_argument("--mode", choices=["disease-to-drugs", "drug-to-diseases", "pair-check"], required=True)
    p.add_argument("--drug", help="Drug name (INN or brand) or ChEMBL ID")
    p.add_argument("--disease", help="Disease name or EFO ID")
    p.add_argument("--max-targets", type=int, default=30)
    p.add_argument("--max-drugs", type=int, default=20)
    p.add_argument("--max-diseases", type=int, default=20)
    p.add_argument("--timeout", type=int, default=30)
    p.add_argument("--output", default="repurposing_results.csv")
    p.add_argument("--summary", default="repurposing_summary.json")
    p.add_argument("--dossier", help="Optional markdown dossier output")
    return p.parse_args()


def require_requests():
    if requests is None:
        raise SystemExit("requests library required for network queries")
    return requests


def ot_graphql(query: str, variables: dict[str, Any], timeout: int) -> Any:
    req = require_requests()
    resp = req.post(OPENTARGETS_URL, json={"query": query, "variables": variables},
                    headers={"Content-Type": "application/json"}, timeout=timeout)
    resp.raise_for_status()
    return resp.json()


def resolve_disease(name: str, timeout: int) -> dict[str, Any]:
    if name.startswith("EFO_") or name.startswith("MONDO_") or name.startswith("Orphanet_"):
        data = ot_graphql(DISEASE_TARGETS_QUERY, {"efoId": name, "size": 1}, timeout)
        disease = data.get("data", {}).get("disease") or {}
        return {"efo_id": disease.get("id", name), "name": disease.get("name", name)}

    data = ot_graphql(DISEASE_SEARCH_QUERY, {"q": name, "size": 5}, timeout)
    hits = data.get("data", {}).get("search", {}).get("hits", [])
    if not hits:
        return {"efo_id": "", "name": name, "error": "no OpenTargets disease match"}
    best = hits[0]
    return {"efo_id": best["id"], "name": best.get("name", name)}


def fetch_disease_targets(efo_id: str, limit: int, timeout: int) -> list[dict[str, Any]]:
    if not efo_id:
        return []
    data = ot_graphql(DISEASE_TARGETS_QUERY, {"efoId": efo_id, "size": limit}, timeout)
    rows_raw = data.get("data", {}).get("disease", {}).get("associatedTargets", {}).get("rows", [])
    return [{"ensembl_id": r["target"]["id"], "symbol": r["target"]["approvedSymbol"],
             "ot_score": r.get("score", 0)} for r in rows_raw]


def resolve_drug_targets(drug_name: str, timeout: int) -> tuple[str, list[dict[str, Any]]]:
    req = require_requests()
    resp = req.get(f"{CHEMBL_BASE}/molecule/search.json",
                   params={"q": drug_name, "limit": 3}, timeout=timeout)
    if resp.status_code != 200:
        return "", []
    mols = resp.json().get("molecules", [])
    if not mols:
        return "", []
    chembl_id = mols[0].get("molecule_chembl_id", "")

    resp = req.get(f"{CHEMBL_BASE}/mechanism.json",
                   params={"molecule_chembl_id": chembl_id, "limit": 50}, timeout=timeout)
    if resp.status_code != 200:
        return chembl_id, []
    mechanisms = resp.json().get("mechanisms", [])
    targets = []
    seen = set()
    for mech in mechanisms:
        tid = mech.get("target_chembl_id", "")
        if tid in seen:
            continue
        seen.add(tid)
        targets.append({
            "target_chembl_id": tid,
            "target_name": mech.get("target_pref_name", ""),
            "mechanism": mech.get("mechanism_of_action", ""),
            "action_type": mech.get("action_type", ""),
        })

    for t in targets:
        tid = t["target_chembl_id"]
        resp2 = req.get(f"{CHEMBL_BASE}/target/{tid}.json", timeout=timeout)
        if resp2.status_code == 200:
            components = resp2.json().get("target_components") or []
            for comp in components:
                for xref in comp.get("target_component_xrefs") or []:
                    if xref.get("xref_src_db") == "EnsemblGene":
                        t["ensembl_id"] = xref.get("xref_id", "")
                        break
        if "ensembl_id" not in t:
            name = t["target_name"]
            ot_data = ot_graphql(TARGET_SEARCH_QUERY, {"q": name, "size": 3}, timeout)
            ot_hits = ot_data.get("data", {}).get("search", {}).get("hits", [])
            if ot_hits:
                t["ensembl_id"] = ot_hits[0].get("id", "")
                t["symbol"] = ot_hits[0].get("name", name)
            else:
                t["ensembl_id"] = ""
                t["symbol"] = name

    return chembl_id, targets


def search_clinical_trials(drug: str, disease: str, timeout: int) -> list[dict[str, Any]]:
    req = require_requests()
    params = {"query.term": f"{drug} {disease}", "pageSize": 20, "format": "json"}
    try:
        resp = req.get(CTGOV_API, params=params, timeout=timeout)
        if resp.status_code != 200:
            return []
        data = resp.json()
    except Exception:
        return []

    trials = []
    for study in data.get("studies", []):
        proto = study.get("protocolSection", {})
        ident = proto.get("identificationModule", {})
        status_mod = proto.get("statusModule", {})
        design = proto.get("designModule", {})
        phases = design.get("phases") or []
        trials.append({
            "nct_id": ident.get("nctId", ""),
            "title": ident.get("briefTitle", ""),
            "status": status_mod.get("overallStatus", ""),
            "phase": ", ".join(phases) if phases else "N/A",
        })
    return trials


def compute_overlap(drug_targets: list[dict], disease_targets: list[dict]) -> list[dict[str, Any]]:
    drug_ensembl = {}
    for t in drug_targets:
        eid = t.get("ensembl_id", "")
        if eid:
            drug_ensembl[eid] = t

    overlaps = []
    for dt in disease_targets:
        eid = dt.get("ensembl_id", "")
        if eid and eid in drug_ensembl:
            drt = drug_ensembl[eid]
            overlaps.append({
                "ensembl_id": eid,
                "symbol": dt.get("symbol", drt.get("symbol", "")),
                "mechanism": drt.get("mechanism", ""),
                "action_type": drt.get("action_type", ""),
                "disease_ot_score": dt.get("ot_score", 0),
            })
    overlaps.sort(key=lambda x: x.get("disease_ot_score", 0), reverse=True)
    return overlaps


def score_pair(overlaps: list[dict], disease_targets: list[dict],
               trials: list[dict]) -> dict[str, Any]:
    n_disease = max(len(disease_targets), 1)
    overlap_fraction = len(overlaps) / n_disease
    avg_ot = sum(o.get("disease_ot_score", 0) for o in overlaps) / max(len(overlaps), 1)
    trial_signal = min(len(trials) / 5, 1.0)
    composite = round(0.4 * overlap_fraction + 0.3 * avg_ot + 0.3 * trial_signal, 4)
    return {
        "overlap_count": len(overlaps),
        "overlap_fraction": round(overlap_fraction, 4),
        "avg_ot_score": round(avg_ot, 4),
        "trial_count": len(trials),
        "composite_score": composite,
    }


def render_dossier(mode: str, drug: str, disease: str, drug_chembl: str,
                   disease_info: dict, drug_targets: list, disease_targets: list,
                   overlaps: list, trials: list, score: dict) -> str:
    lines = [f"# Drug Repurposing Brief", ""]
    lines.append(f"**Mode**: {mode}")
    if drug:
        lines.append(f"**Drug**: {drug} (`{drug_chembl}`)")
    if disease:
        lines.append(f"**Disease**: {disease_info.get('name', disease)} (`{disease_info.get('efo_id', '')}`)")
    lines.append("")

    lines.append("## Target Analysis")
    lines.append(f"- Drug targets resolved: {len(drug_targets)}")
    lines.append(f"- Disease-associated targets: {len(disease_targets)}")
    lines.append(f"- **Overlapping targets: {len(overlaps)}**")
    lines.append("")
    if overlaps:
        lines.append("### Overlapping Targets")
        lines.append("| Symbol | Mechanism | Action | OT Score |")
        lines.append("|--------|-----------|--------|----------|")
        for o in overlaps:
            lines.append(f"| {o.get('symbol', '')} | {o.get('mechanism', '')} | "
                         f"{o.get('action_type', '')} | {o.get('disease_ot_score', '')} |")
        lines.append("")

    lines.append("## Clinical Evidence")
    lines.append(f"- ClinicalTrials.gov hits: {len(trials)}")
    if trials:
        for t in trials[:10]:
            lines.append(f"  - [{t.get('nct_id', '')}] {t.get('title', '')} — "
                         f"Phase: {t.get('phase', 'N/A')}, Status: {t.get('status', '')}")
    else:
        lines.append("  - No registered trials found for this drug-disease combination")
    lines.append("")

    lines.append("## Composite Score")
    lines.append(f"- Target overlap fraction: {score.get('overlap_fraction', 0)}")
    lines.append(f"- Average OT association score: {score.get('avg_ot_score', 0)}")
    lines.append(f"- Clinical trial signal: {score.get('trial_count', 0)} trials")
    lines.append(f"- **Composite**: {score.get('composite_score', 0)}")
    lines.append("")
    lines.append("## Caveats")
    lines.append("- This is a research hypothesis, not clinical evidence.")
    lines.append("- Target overlap does not account for tissue specificity or drug exposure.")
    lines.append("- Clinical trials may be unrelated to the mechanism of interest.")
    lines.append("- Genetic association scores are probabilistic, not deterministic.")
    return "\n".join(lines) + "\n"


def disease_to_drugs(args) -> dict[str, Any]:
    print(f"Resolving disease: {args.disease}")
    disease_info = resolve_disease(args.disease, args.timeout)
    efo_id = disease_info.get("efo_id", "")
    if not efo_id:
        raise SystemExit(f"Could not resolve disease: {args.disease}")
    print(f"  EFO ID: {efo_id} ({disease_info.get('name', '')})")

    disease_targets = fetch_disease_targets(efo_id, args.max_targets, args.timeout)
    print(f"  Disease targets: {len(disease_targets)}")

    req = require_requests()
    drug_rows = []
    disease_ensembl = {dt["ensembl_id"] for dt in disease_targets if dt.get("ensembl_id")}

    known_drugs_q = f"""
    query KnownDrugs($efoId: String!, $size: Int!) {{
      disease(efoId: $efoId) {{
        knownDrugs(size: $size) {{
          count
          rows {{ drug {{ id name }} mechanismOfAction phase status
                  target {{ id approvedSymbol }} }}
        }}
      }}
    }}"""
    ot = ot_graphql(known_drugs_q, {"efoId": efo_id, "size": 100}, args.timeout)
    known_rows = ot.get("data", {}).get("disease", {}).get("knownDrugs", {}).get("rows", [])

    seen_drugs: dict[str, dict] = {}
    for kr in known_rows:
        drug_obj = kr.get("drug") or {}
        did = drug_obj.get("id", "")
        if not did or did in seen_drugs:
            continue
        seen_drugs[did] = {
            "drug_chembl_id": did,
            "drug_name": drug_obj.get("name", ""),
            "mechanism": kr.get("mechanismOfAction", ""),
            "phase": kr.get("phase", ""),
            "target_id": (kr.get("target") or {}).get("id", ""),
            "target_symbol": (kr.get("target") or {}).get("approvedSymbol", ""),
        }

    for did, info in list(seen_drugs.items())[:args.max_drugs]:
        target_eid = info.get("target_id", "")
        in_disease = target_eid in disease_ensembl
        ot_score = next((dt["ot_score"] for dt in disease_targets if dt["ensembl_id"] == target_eid), 0)
        trials = search_clinical_trials(info["drug_name"], args.disease, args.timeout)
        drug_rows.append({
            **info,
            "target_in_disease_targets": in_disease,
            "disease_ot_score": ot_score,
            "trial_count": len(trials),
        })

    drug_rows.sort(key=lambda x: (x.get("disease_ot_score", 0), x.get("trial_count", 0)), reverse=True)
    return {
        "mode": "disease-to-drugs",
        "disease": disease_info,
        "disease_targets": disease_targets,
        "drug_rows": drug_rows,
    }


def drug_to_diseases(args) -> dict[str, Any]:
    print(f"Resolving drug targets: {args.drug}")
    chembl_id, drug_targets = resolve_drug_targets(args.drug, args.timeout)
    print(f"  ChEMBL ID: {chembl_id}, targets: {len(drug_targets)}")

    disease_rows = []
    for t in drug_targets:
        eid = t.get("ensembl_id", "")
        if not eid:
            continue
        ot = ot_graphql(TARGET_DISEASES_QUERY, {"ensemblId": eid, "size": args.max_diseases}, args.timeout)
        target_data = ot.get("data", {}).get("target", {}) or {}
        for row in target_data.get("associatedDiseases", {}).get("rows", []):
            disease = row.get("disease") or {}
            disease_rows.append({
                "disease_id": disease.get("id", ""),
                "disease_name": disease.get("name", ""),
                "via_target": t.get("symbol", eid),
                "mechanism": t.get("mechanism", ""),
                "ot_score": row.get("score", 0),
            })

    seen = {}
    for dr in disease_rows:
        did = dr["disease_id"]
        if did not in seen or dr["ot_score"] > seen[did]["ot_score"]:
            seen[did] = dr
    unique_diseases = sorted(seen.values(), key=lambda x: x.get("ot_score", 0), reverse=True)

    for d in unique_diseases[:args.max_diseases]:
        trials = search_clinical_trials(args.drug, d["disease_name"], args.timeout)
        d["trial_count"] = len(trials)

    return {
        "mode": "drug-to-diseases",
        "drug": args.drug,
        "drug_chembl_id": chembl_id,
        "drug_targets": drug_targets,
        "disease_rows": unique_diseases[:args.max_diseases],
    }


def pair_check(args) -> dict[str, Any]:
    print(f"Pair check: {args.drug} → {args.disease}")
    chembl_id, drug_targets = resolve_drug_targets(args.drug, args.timeout)
    disease_info = resolve_disease(args.disease, args.timeout)
    efo_id = disease_info.get("efo_id", "")
    disease_targets = fetch_disease_targets(efo_id, args.max_targets, args.timeout) if efo_id else []
    overlaps = compute_overlap(drug_targets, disease_targets)
    trials = search_clinical_trials(args.drug, args.disease, args.timeout)
    sc = score_pair(overlaps, disease_targets, trials)

    return {
        "mode": "pair-check",
        "drug": args.drug,
        "drug_chembl_id": chembl_id,
        "disease": disease_info,
        "drug_targets": drug_targets,
        "disease_targets": disease_targets,
        "overlaps": overlaps,
        "trials": trials,
        "score": sc,
    }


def write_csv_rows(path: str, rows: list[dict]) -> None:
    if not rows:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        Path(path).write_text("no_results\n")
        return
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    keys = list(rows[0].keys())
    for r in rows[1:]:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    args = parse_args()

    if args.mode == "disease-to-drugs":
        if not args.disease:
            raise SystemExit("--disease required for disease-to-drugs mode")
        result = disease_to_drugs(args)
        write_csv_rows(args.output, result["drug_rows"])
        summary = {
            "mode": result["mode"],
            "disease": result["disease"],
            "disease_target_count": len(result["disease_targets"]),
            "drug_candidates": len(result["drug_rows"]),
            "output_file": args.output,
        }

    elif args.mode == "drug-to-diseases":
        if not args.drug:
            raise SystemExit("--drug required for drug-to-diseases mode")
        result = drug_to_diseases(args)
        write_csv_rows(args.output, result["disease_rows"])
        summary = {
            "mode": result["mode"],
            "drug": result["drug"],
            "drug_chembl_id": result["drug_chembl_id"],
            "drug_target_count": len(result["drug_targets"]),
            "disease_candidates": len(result["disease_rows"]),
            "output_file": args.output,
        }

    elif args.mode == "pair-check":
        if not args.drug or not args.disease:
            raise SystemExit("--drug and --disease required for pair-check")
        result = pair_check(args)
        write_csv_rows(args.output, result["overlaps"] or [{"note": "no target overlap found"}])
        summary = {
            "mode": result["mode"],
            "drug": result["drug"],
            "drug_chembl_id": result["drug_chembl_id"],
            "disease": result["disease"],
            "drug_target_count": len(result["drug_targets"]),
            "disease_target_count": len(result["disease_targets"]),
            "score": result["score"],
            "trial_count": len(result["trials"]),
            "output_file": args.output,
        }
        if args.dossier:
            md = render_dossier(
                result["mode"], args.drug, args.disease, result["drug_chembl_id"],
                result["disease"], result["drug_targets"], result["disease_targets"],
                result["overlaps"], result["trials"], result["score"])
            Path(args.dossier).parent.mkdir(parents=True, exist_ok=True)
            Path(args.dossier).write_text(md, encoding="utf-8")
            print(f"Dossier: {args.dossier}")
    else:
        raise SystemExit(f"Unknown mode: {args.mode}")

    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary).write_text(json.dumps(summary, indent=2, default=str, ensure_ascii=False))
    print(f"Summary: {args.summary}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
