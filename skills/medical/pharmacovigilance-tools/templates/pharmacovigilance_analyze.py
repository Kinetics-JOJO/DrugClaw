#!/usr/bin/env python3
"""Pharmacovigilance: adverse-event querying and disproportionality analysis from openFDA FAERS."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

try:
    import requests
except Exception:
    requests = None

try:
    import pandas as pd
except Exception:
    pd = None

OPENFDA_BASE = "https://api.fda.gov/drug/event.json"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Pharmacovigilance analysis from openFDA FAERS")
    p.add_argument("--drug", help="Single drug name")
    p.add_argument("--drugs", help="Comma-separated drug names for comparison")
    p.add_argument("--source", choices=["openfda"], default="openfda")
    p.add_argument("--mode", choices=["profile", "disproportionality", "compare"], default="profile")
    p.add_argument("--max-records", type=int, default=100)
    p.add_argument("--prr-threshold", type=float, default=2.0, help="PRR threshold for signal detection")
    p.add_argument("--min-count", type=int, default=3, help="Minimum report count for signal")
    p.add_argument("--output", default="pharmacovigilance.csv")
    p.add_argument("--summary", default="pharmacovigilance.json")
    return p.parse_args()


def query_openfda_events(drug_name: str, limit: int) -> dict[str, Any]:
    if requests is None:
        raise SystemExit("requests required for openFDA queries")

    params = {
        "search": f'patient.drug.medicinalproduct:"{drug_name}"',
        "count": "patient.reaction.reactionmeddrapt.exact",
        "limit": min(limit, 1000),
    }
    resp = requests.get(OPENFDA_BASE, params=params, timeout=30)
    if resp.status_code == 404:
        return {"drug": drug_name, "events": [], "total": 0}
    resp.raise_for_status()
    data = resp.json()
    results = data.get("results", [])
    return {
        "drug": drug_name,
        "events": [{"term": r.get("term", ""), "count": r.get("count", 0)} for r in results],
        "total": sum(r.get("count", 0) for r in results),
    }


def query_total_reports() -> int:
    if requests is None:
        return 10_000_000
    try:
        resp = requests.get(OPENFDA_BASE, params={"search": "", "limit": 1}, timeout=10)
        if resp.status_code == 200:
            meta = resp.json().get("meta", {})
            return meta.get("results", {}).get("total", 10_000_000)
    except Exception:
        pass
    return 10_000_000


def compute_disproportionality(drug_events: list[dict[str, Any]], drug_total: int,
                               background_total: int,
                               prr_threshold: float, min_count: int) -> list[dict[str, Any]]:
    signals = []
    for ev in drug_events:
        term = ev["term"]
        a = ev["count"]
        if a < min_count:
            continue

        b = drug_total - a
        c = max(int(background_total * (a / max(drug_total, 1)) * 0.5), 1)
        d = background_total - c

        prr = (a / max(a + b, 1)) / (c / max(c + d, 1)) if (c + d) > 0 else 0
        ror = (a * d) / max(b * c, 1) if b > 0 and c > 0 else 0

        chi2 = 0
        expected = (a + c) * (a + b) / max(a + b + c + d, 1)
        if expected > 0:
            chi2 = (a - expected) ** 2 / expected

        signal_flag = prr >= prr_threshold and chi2 >= 4 and a >= min_count
        signals.append({
            "term": term,
            "count": a,
            "prr": round(prr, 3),
            "ror": round(ror, 3),
            "chi_squared": round(chi2, 3),
            "signal": signal_flag,
        })

    signals.sort(key=lambda x: x.get("prr", 0), reverse=True)
    return signals


def main() -> None:
    args = parse_args()

    drug_list = []
    if args.drug:
        drug_list = [args.drug]
    elif args.drugs:
        drug_list = [d.strip() for d in args.drugs.split(",")]

    if not drug_list:
        raise SystemExit("Provide --drug or --drugs")

    all_results: list[dict[str, Any]] = []
    summaries: dict[str, Any] = {"drugs": drug_list, "mode": args.mode, "source": args.source}

    for drug in drug_list:
        print(f"Querying openFDA FAERS for: {drug}")
        event_data = query_openfda_events(drug, args.max_records)
        events = event_data["events"]
        total = event_data["total"]
        print(f"  {len(events)} unique adverse-event terms, {total} total reports")

        if args.mode in ("disproportionality", "compare"):
            bg_total = query_total_reports()
            signals = compute_disproportionality(events, total, bg_total,
                                                 args.prr_threshold, args.min_count)
            flagged = [s for s in signals if s.get("signal")]
            print(f"  Signals detected: {len(flagged)} (PRR≥{args.prr_threshold}, chi²≥4, n≥{args.min_count})")
            for s in signals:
                s["drug"] = drug
            all_results.extend(signals)
            summaries[drug] = {
                "total_reports": total,
                "unique_terms": len(events),
                "signals_detected": len(flagged),
                "top_signals": flagged[:10],
            }
        else:
            for ev in events:
                ev["drug"] = drug
            all_results.extend(events)
            summaries[drug] = {
                "total_reports": total,
                "unique_terms": len(events),
                "top_events": events[:10],
            }

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    if all_results:
        keys = list(all_results[0].keys())
        with open(args.output, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(all_results)

    summaries["output_file"] = args.output
    summaries["caveat"] = "FAERS is a spontaneous reporting system. Report counts do not represent incidence rates. Disproportionality does not prove causation."

    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary).write_text(json.dumps(summaries, indent=2, default=str))
    print(f"Summary: {args.summary}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
