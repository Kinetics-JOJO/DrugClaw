#!/usr/bin/env python3
"""Pharmacokinetic analysis: NCA, compartmental fitting, and dose-response simulation."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="PK analysis: NCA, compartmental fit, simulation")
    p.add_argument("--input", help="Concentration-time CSV (nca and compartmental modes)")
    p.add_argument("--time-column", default="time_h")
    p.add_argument("--conc-column", default="concentration_ng_ml")
    p.add_argument("--group-column", help="Subject or group column")
    p.add_argument("--dose", type=float, help="Dose amount")
    p.add_argument("--dose-unit", default="mg/kg")
    p.add_argument("--mode", choices=["nca", "compartmental", "simulate"], default="nca")

    p.add_argument("--compartments", type=int, default=1, choices=[1, 2])
    p.add_argument("--route", choices=["iv-bolus", "oral"], default="iv-bolus")

    p.add_argument("--ka", type=float, help="Absorption rate constant (simulate, oral)")
    p.add_argument("--ke", type=float, help="Elimination rate constant (simulate)")
    p.add_argument("--vd", type=float, help="Volume of distribution (simulate)")
    p.add_argument("--bioavailability", type=float, default=1.0)
    p.add_argument("--interval", type=float, default=24, help="Dosing interval in hours")
    p.add_argument("--doses", type=int, default=1, help="Number of doses")
    p.add_argument("--sim-duration", type=float, help="Total simulation hours (default: doses * interval * 1.5)")
    p.add_argument("--sim-step", type=float, default=0.1, help="Time step in hours")

    p.add_argument("--output", default="pk_results.csv")
    p.add_argument("--summary", default="pk_summary.json")
    return p.parse_args()


def require_scipy():
    try:
        import numpy as _np
        from scipy import optimize as _opt
        return _np, _opt
    except ImportError as e:
        raise SystemExit(f"numpy and scipy required: {e}")


def read_pk_data(path: str, time_col: str, conc_col: str,
                 group_col: str | None) -> dict[str, list[tuple[float, float]]]:
    p = Path(path)
    delimiter = "\t" if p.suffix.lower() in (".tsv", ".tab") else ","
    with open(p, newline="") as f:
        rows = list(csv.DictReader(f, delimiter=delimiter))

    groups: dict[str, list[tuple[float, float]]] = {}
    for r in rows:
        g = r.get(group_col, "all") if group_col else "all"
        try:
            t = float(r[time_col])
            c = float(r[conc_col])
        except (ValueError, KeyError):
            continue
        groups.setdefault(g, []).append((t, c))

    for g in groups:
        groups[g].sort()
    return groups


def nca_analysis(times: list[float], concs: list[float], dose: float | None) -> dict[str, Any]:
    np, _ = require_scipy()

    if len(times) < 3:
        return {"error": "need at least 3 time points"}

    cmax = max(concs)
    tmax = times[concs.index(cmax)]

    auc = 0.0
    for i in range(len(times) - 1):
        dt = times[i + 1] - times[i]
        c1, c2 = concs[i], concs[i + 1]
        if c1 > 0 and c2 > 0 and c2 < c1:
            auc += dt * (c1 - c2) / math.log(c1 / c2)
        else:
            auc += dt * (c1 + c2) / 2.0

    terminal_points = [(times[i], concs[i]) for i in range(len(times)) if concs[i] > 0]
    half_life = None
    ke = None
    if len(terminal_points) >= 3:
        t_term = [p[0] for p in terminal_points[-3:]]
        c_term = [math.log(p[1]) for p in terminal_points[-3:]]
        n = len(t_term)
        sum_t = sum(t_term)
        sum_c = sum(c_term)
        sum_tc = sum(t_term[i] * c_term[i] for i in range(n))
        sum_t2 = sum(t ** 2 for t in t_term)
        denom = n * sum_t2 - sum_t ** 2
        if abs(denom) > 1e-12:
            slope = (n * sum_tc - sum_t * sum_c) / denom
            if slope < 0:
                ke = -slope
                half_life = math.log(2) / ke

    result: dict[str, Any] = {
        "cmax": round(cmax, 4),
        "tmax": round(tmax, 4),
        "auc_0_last": round(auc, 4),
    }
    if half_life is not None:
        result["half_life_h"] = round(half_life, 4)
        result["ke_h"] = round(ke, 6)
        if ke and concs[-1] > 0:
            auc_inf = auc + concs[-1] / ke
            result["auc_0_inf"] = round(auc_inf, 4)
    if dose is not None and auc > 0:
        result["clearance"] = round(dose / auc, 6)
        if half_life and ke:
            result["vd"] = round(dose / (auc * ke), 4)
    return result


def simulate_pk(args) -> list[dict[str, Any]]:
    np, _ = require_scipy()

    duration = args.sim_duration or args.doses * args.interval * 1.5
    times = np.arange(0, duration, args.sim_step)
    conc = np.zeros_like(times)

    for d in range(args.doses):
        t_dose = d * args.interval
        mask = times >= t_dose
        t_rel = times[mask] - t_dose
        effective_dose = args.dose * args.bioavailability

        if args.route == "iv-bolus":
            c = (effective_dose / args.vd) * np.exp(-args.ke * t_rel)
        else:
            if args.ka is None:
                raise SystemExit("--ka required for oral route simulation")
            if abs(args.ka - args.ke) < 1e-12:
                c = (effective_dose * args.ka / args.vd) * t_rel * np.exp(-args.ke * t_rel)
            else:
                c = (effective_dose * args.ka / (args.vd * (args.ka - args.ke))) * (
                    np.exp(-args.ke * t_rel) - np.exp(-args.ka * t_rel))
        conc[mask] += c

    return [{"time_h": round(float(t), 4), "concentration": round(float(c), 6)}
            for t, c in zip(times, conc)]


def main() -> None:
    args = parse_args()

    if args.mode == "nca":
        if not args.input:
            raise SystemExit("--input required for NCA")
        groups = read_pk_data(args.input, args.time_column, args.conc_column, args.group_column)
        print(f"NCA analysis for {len(groups)} group(s)")

        all_results = []
        for gname, data in groups.items():
            times = [d[0] for d in data]
            concs = [d[1] for d in data]
            nca = nca_analysis(times, concs, args.dose)
            nca["group"] = gname
            nca["n_timepoints"] = len(times)
            if args.dose is not None:
                nca["dose"] = args.dose
                nca["dose_unit"] = args.dose_unit
            all_results.append(nca)

        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        if all_results:
            keys = list(all_results[0].keys())
            for r in all_results[1:]:
                for k in r:
                    if k not in keys:
                        keys.append(k)
            with open(args.output, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
                w.writeheader()
                w.writerows(all_results)

        summary = {"mode": "nca", "groups": len(groups), "results": all_results, "output_file": args.output}

    elif args.mode == "compartmental":
        if not args.input:
            raise SystemExit("--input required for compartmental fit")
        np, opt = require_scipy()
        groups = read_pk_data(args.input, args.time_column, args.conc_column, args.group_column)
        print(f"Compartmental ({args.compartments}-cmpt, {args.route}) fit for {len(groups)} group(s)")

        all_results = []
        for gname, data in groups.items():
            times = np.array([d[0] for d in data])
            concs = np.array([d[1] for d in data])

            if args.compartments == 1 and args.route == "iv-bolus":
                def model(t, c0, ke):
                    return c0 * np.exp(-ke * t)
                try:
                    popt, _ = opt.curve_fit(model, times, concs, p0=[concs[0], 0.1], maxfev=5000)
                    pred = model(times, *popt)
                    ss_res = np.sum((concs - pred) ** 2)
                    ss_tot = np.sum((concs - np.mean(concs)) ** 2)
                    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
                    rmse = float(np.sqrt(ss_res / len(concs)))
                    result = {
                        "group": gname, "model": "1cmpt-iv",
                        "C0": round(float(popt[0]), 4), "ke": round(float(popt[1]), 6),
                        "half_life_h": round(float(math.log(2) / popt[1]), 4) if popt[1] > 0 else None,
                        "r_squared": round(r2, 4), "rmse": round(rmse, 4),
                    }
                except Exception as e:
                    result = {"group": gname, "model": "1cmpt-iv", "error": str(e)}
            else:
                result = {"group": gname, "model": f"{args.compartments}cmpt-{args.route}",
                          "note": "Extended compartmental models require additional parameter constraints"}
            all_results.append(result)

        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        if all_results:
            keys = list(all_results[0].keys())
            with open(args.output, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
                w.writeheader()
                w.writerows(all_results)

        summary = {"mode": "compartmental", "compartments": args.compartments, "route": args.route,
                    "groups": len(groups), "results": all_results, "output_file": args.output}

    elif args.mode == "simulate":
        if not args.dose or not args.ke or not args.vd:
            raise SystemExit("--dose, --ke, --vd required for simulation")
        print(f"Simulating {args.compartments}-cmpt {args.route} PK: {args.doses} dose(s) of {args.dose} {args.dose_unit}")
        curve = simulate_pk(args)

        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["time_h", "concentration"])
            w.writeheader()
            w.writerows(curve)

        cmax = max(r["concentration"] for r in curve)
        tmax = next(r["time_h"] for r in curve if r["concentration"] == cmax)
        summary = {
            "mode": "simulate", "compartments": args.compartments, "route": args.route,
            "dose": args.dose, "dose_unit": args.dose_unit, "doses": args.doses,
            "interval_h": args.interval, "ke": args.ke, "vd": args.vd,
            "bioavailability": args.bioavailability,
            "cmax_simulated": round(cmax, 4), "tmax_simulated": round(tmax, 4),
            "duration_h": round(curve[-1]["time_h"], 2),
            "output_file": args.output,
        }
        if args.ka:
            summary["ka"] = args.ka
    else:
        raise SystemExit(f"Unknown mode: {args.mode}")

    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary).write_text(json.dumps(summary, indent=2, default=str))
    print(f"Summary: {args.summary}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
