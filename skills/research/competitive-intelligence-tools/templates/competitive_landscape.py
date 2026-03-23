#!/usr/bin/env python3
"""Pharmaceutical competitive landscape from ClinicalTrials.gov v2 and OpenAlex."""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import requests
except Exception:
    requests = None

CTGOV_API = "https://clinicaltrials.gov/api/v2/studies"
OPENALEX_API = "https://api.openalex.org/works"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Competitive landscape from clinical-trial and publication data")
    p.add_argument("--target", help="Target gene symbol or mechanism keyword")
    p.add_argument("--indication", help="Disease or indication term")
    p.add_argument("--phase", help="Phase filter: EARLY_PHASE1, PHASE1, PHASE2, PHASE3, PHASE4")
    p.add_argument("--status", help="Status filter: RECRUITING, ACTIVE_NOT_RECRUITING, COMPLETED, etc.")
    p.add_argument("--max-trials", type=int, default=200)
    p.add_argument("--include-publications", action="store_true")
    p.add_argument("--max-publications", type=int, default=50)
    p.add_argument("--timeout", type=int, default=30)
    p.add_argument("--output", default="competitive_trials.csv")
    p.add_argument("--summary", default="competitive_summary.json")
    p.add_argument("--brief", help="Optional markdown brief output")
    return p.parse_args()


def require_requests():
    if requests is None:
        raise SystemExit("requests library required")
    return requests


def fetch_trials(query_parts: list[str], phase: str | None, status: str | None,
                 max_results: int, timeout: int) -> list[dict[str, Any]]:
    req = require_requests()
    query = " ".join(query_parts)
    all_trials: list[dict[str, Any]] = []
    page_token = None
    remaining = max_results

    while remaining > 0:
        params: dict[str, Any] = {
            "query.term": query,
            "pageSize": min(remaining, 100),
            "format": "json",
        }
        if phase:
            params["filter.phase"] = phase
        if status:
            params["filter.overallStatus"] = status
        if page_token:
            params["pageToken"] = page_token

        resp = req.get(CTGOV_API, params=params, timeout=timeout)
        if resp.status_code != 200:
            break
        data = resp.json()

        for study in data.get("studies", []):
            proto = study.get("protocolSection", {})
            ident = proto.get("identificationModule", {})
            status_mod = proto.get("statusModule", {})
            design = proto.get("designModule", {})
            sponsor_mod = proto.get("sponsorCollaboratorsModule", {})
            lead = sponsor_mod.get("leadSponsor", {})
            desc = proto.get("descriptionModule", {})
            enroll_mod = design.get("enrollmentInfo", {})
            conds = proto.get("conditionsModule", {})
            arms_mod = proto.get("armsInterventionsModule", {})

            phases = design.get("phases") or []
            conditions = conds.get("conditions") or []
            interventions = []
            for arm in (arms_mod.get("interventions") or []):
                interventions.append(arm.get("name", ""))

            first_posted = status_mod.get("studyFirstPostDateStruct", {}).get("date", "")
            start_date = status_mod.get("startDateStruct", {}).get("date", "")

            all_trials.append({
                "nct_id": ident.get("nctId", ""),
                "title": ident.get("briefTitle", ""),
                "status": status_mod.get("overallStatus", ""),
                "phase": ", ".join(phases) if phases else "N/A",
                "sponsor": lead.get("name", ""),
                "sponsor_class": lead.get("class", ""),
                "enrollment": enroll_mod.get("count", ""),
                "conditions": "; ".join(conditions[:5]),
                "interventions": "; ".join(interventions[:5]),
                "start_date": start_date,
                "first_posted": first_posted,
                "modality_guess": guess_modality("; ".join(interventions)),
            })

        page_token = data.get("nextPageToken")
        remaining -= len(data.get("studies", []))
        if not page_token or not data.get("studies"):
            break

    return all_trials


MODALITY_PATTERNS = [
    (r'\b(mab|umab|izumab|ximab|zumab)\b', "antibody"),
    (r'\bADC\b|antibody.drug.conjugate', "ADC"),
    (r'\bCAR.?T\b|chimeric antigen receptor', "cell therapy"),
    (r'\bgene.therapy\b|AAV|lentivir', "gene therapy"),
    (r'\bsiRNA\b|\bASO\b|antisense|oligonucleotide', "oligonucleotide"),
    (r'\bvaccine\b|mRNA.vaccine', "vaccine"),
    (r'\bdegrader\b|PROTAC|molecular.glue', "degrader"),
    (r'\bbispecific\b', "bispecific"),
    (r'\bradio.?(?:therapy|pharma|conjugate|ligand)\b|\bRLT\b', "radiopharmaceutical"),
]


def guess_modality(intervention_text: str) -> str:
    text = intervention_text.lower()
    for pattern, label in MODALITY_PATTERNS:
        if re.search(pattern, text, re.IGNORECASE):
            return label
    return "small molecule / other"


def fetch_publications(query: str, max_results: int, timeout: int) -> list[dict[str, Any]]:
    req = require_requests()
    params = {"search": query, "per_page": min(max_results, 100), "sort": "publication_date:desc"}
    try:
        resp = req.get(OPENALEX_API, params=params, timeout=timeout)
        if resp.status_code != 200:
            return []
        data = resp.json()
    except Exception:
        return []

    pubs = []
    for work in data.get("results", []):
        authorships = work.get("authorships") or []
        institutions = set()
        for a in authorships:
            for inst in a.get("institutions") or []:
                name = inst.get("display_name", "")
                if name:
                    institutions.add(name)
        pubs.append({
            "title": work.get("title", ""),
            "publication_date": work.get("publication_date", ""),
            "cited_by_count": work.get("cited_by_count", 0),
            "institutions": "; ".join(sorted(institutions)[:5]),
            "doi": work.get("doi", ""),
        })
    return pubs


def analyze_trials(trials: list[dict]) -> dict[str, Any]:
    phase_counts = Counter(t["phase"] for t in trials)
    sponsor_counts = Counter(t["sponsor"] for t in trials if t["sponsor"])
    status_counts = Counter(t["status"] for t in trials)
    modality_counts = Counter(t.get("modality_guess", "") for t in trials)

    sponsor_phases: dict[str, str] = {}
    phase_rank = {"PHASE4": 5, "PHASE3": 4, "PHASE2": 3, "PHASE1": 2, "EARLY_PHASE1": 1}
    for t in trials:
        sp = t["sponsor"]
        ph = t["phase"].replace(", ", "/")
        best_rank = max((phase_rank.get(p.strip(), 0) for p in ph.split("/")), default=0)
        if sp not in sponsor_phases or best_rank > phase_rank.get(sponsor_phases[sp], 0):
            sponsor_phases[sp] = ph

    top_sponsors = []
    for sp, count in sponsor_counts.most_common(15):
        top_sponsors.append({"sponsor": sp, "trial_count": count, "most_advanced_phase": sponsor_phases.get(sp, "")})

    year_counts: dict[str, int] = {}
    for t in trials:
        date = t.get("first_posted") or t.get("start_date") or ""
        if len(date) >= 4:
            yr = date[:4]
            year_counts[yr] = year_counts.get(yr, 0) + 1

    total_enrollment = 0
    for t in trials:
        try:
            total_enrollment += int(t.get("enrollment", 0))
        except (ValueError, TypeError):
            pass

    return {
        "total_trials": len(trials),
        "phase_distribution": dict(sorted(phase_counts.items())),
        "status_distribution": dict(status_counts),
        "modality_distribution": dict(modality_counts),
        "top_sponsors": top_sponsors,
        "filing_trend_by_year": dict(sorted(year_counts.items())),
        "total_enrollment": total_enrollment,
    }


def render_brief(query_desc: str, analysis: dict, trials: list[dict],
                 pubs: list[dict] | None) -> str:
    lines = ["# Competitive Landscape Brief", ""]
    lines.append(f"**Query**: {query_desc}")
    lines.append(f"**Retrieved**: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}")
    lines.append(f"**Total trials**: {analysis['total_trials']}")
    lines.append(f"**Total enrollment**: {analysis['total_enrollment']:,}")
    lines.append("")

    lines.append("## Phase Distribution")
    for ph, cnt in sorted(analysis["phase_distribution"].items()):
        lines.append(f"- {ph}: {cnt}")
    lines.append("")

    lines.append("## Modality Breakdown")
    for mod, cnt in sorted(analysis["modality_distribution"].items(), key=lambda x: -x[1]):
        lines.append(f"- {mod}: {cnt}")
    lines.append("")

    lines.append("## Top Sponsors")
    lines.append("| Sponsor | Trials | Most Advanced |")
    lines.append("|---------|--------|---------------|")
    for sp in analysis["top_sponsors"][:15]:
        lines.append(f"| {sp['sponsor']} | {sp['trial_count']} | {sp['most_advanced_phase']} |")
    lines.append("")

    lines.append("## Filing Trend")
    for yr, cnt in sorted(analysis["filing_trend_by_year"].items()):
        lines.append(f"- {yr}: {cnt} trials")
    lines.append("")

    if pubs:
        lines.append("## Publication Activity")
        lines.append(f"Top {len(pubs)} recent publications:")
        for p in pubs[:10]:
            lines.append(f"- [{p.get('publication_date', '')}] {p.get('title', '')} "
                         f"(cited: {p.get('cited_by_count', 0)})")
        lines.append("")

    lines.append("## Caveats")
    lines.append("- ClinicalTrials.gov coverage is US-centric; EU/China registries not queried.")
    lines.append("- Modality classification is heuristic from intervention text.")
    lines.append("- Sponsor may be CRO or academic center, not the originator company.")
    lines.append("- Trial registration does not indicate clinical success.")
    return "\n".join(lines) + "\n"


def main() -> None:
    args = parse_args()
    query_parts = []
    if args.target:
        query_parts.append(args.target)
    if args.indication:
        query_parts.append(args.indication)
    if not query_parts:
        raise SystemExit("Provide at least --target or --indication")

    query_desc = " + ".join(query_parts)
    print(f"Searching ClinicalTrials.gov: {query_desc}")
    trials = fetch_trials(query_parts, args.phase, args.status, args.max_trials, args.timeout)
    print(f"  Retrieved {len(trials)} trials")
    analysis = analyze_trials(trials)

    pubs = None
    if args.include_publications:
        print(f"Searching OpenAlex for publication signals...")
        pubs = fetch_publications(" ".join(query_parts), args.max_publications, args.timeout)
        print(f"  Retrieved {len(pubs)} publications")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    if trials:
        keys = list(trials[0].keys())
        with open(args.output, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(trials)
    else:
        Path(args.output).write_text("no_results\n")

    summary: dict[str, Any] = {
        "query": query_desc,
        "retrieval_date": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "filters": {"phase": args.phase, "status": args.status},
        **analysis,
        "output_file": args.output,
    }
    if pubs:
        summary["publication_count"] = len(pubs)
        summary["top_publications"] = pubs[:5]

    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary).write_text(json.dumps(summary, indent=2, default=str, ensure_ascii=False))
    print(f"Summary: {args.summary}")

    if args.brief:
        md = render_brief(query_desc, analysis, trials, pubs)
        Path(args.brief).parent.mkdir(parents=True, exist_ok=True)
        Path(args.brief).write_text(md, encoding="utf-8")
        print(f"Brief: {args.brief}")

    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
