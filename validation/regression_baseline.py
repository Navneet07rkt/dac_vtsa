"""Regression check of a run's KPIs against a stored baseline.

    python -m validation.regression_baseline results/<run>/kpis.json [--baseline validation/baseline/<file>.json]

* Baselines are KPI JSON files written by run_dac_simulation.py (``kpis.json``)
  from a run you trust; after any refactor, re-run and compare.
* ``parse_legacy_log`` reads the old ``vtsa_main_run.log`` format ("key = value"
  lines).  NOTE: that log came from the earlier dry-air / RK4 scripts (12.8 %
  CO2 feed, 150 C wall, gas concentration that does not drop during a 10x
  pressure reduction), so it is NOT a physically valid reference for this
  model -- use it only to see what changed, never as a pass/fail target.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from typing import Dict

DEFAULT_TOL = {"recovery_pct": 0.02, "purity_dry_pct": 0.005, "working_capacity_mol_kg": 0.02,
               "productivity_kg_m3_day": 0.02, "specific_thermal_kJ_per_kg": 0.02,
               "specific_electric_kJ_per_kg": 0.03, "specific_total_kJ_per_kg": 0.02}


def compare(run: Dict[str, float], base: Dict[str, float], tol: Dict[str, float] = None):
    tol = tol or DEFAULT_TOL
    rows, ok = [], True
    for k, t in tol.items():
        if k not in run or k not in base:
            continue
        a, b = float(run[k]), float(base[k])
        rel = abs(a - b) / max(abs(b), 1e-30)
        passed = rel <= t
        ok &= passed
        rows.append({"kpi": k, "run": a, "baseline": b, "rel_diff": rel, "tolerance": t, "pass": passed})
    return ok, rows


def parse_legacy_log(path: str) -> Dict[str, float]:
    out = {}
    pat = re.compile(r"^\s*([A-Za-z0-9_ ()\-/,.%₂]+?)\s*=\s*([-+0-9.eE]+)")
    with open(path, errors="ignore") as fh:
        for line in fh:
            m = pat.match(line)
            if m:
                try:
                    out[m.group(1).strip()] = float(m.group(2))
                except ValueError:
                    pass
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("kpis_json")
    ap.add_argument("--baseline", default=os.path.join(os.path.dirname(__file__), "baseline", "default_kpis.json"))
    a = ap.parse_args(argv)
    run = json.load(open(a.kpis_json))
    base = json.load(open(a.baseline))
    ok, rows = compare(run, base)
    for r in rows:
        print(f"{'PASS' if r['pass'] else 'FAIL'}  {r['kpi']:<30s} run={r['run']:.6g}  base={r['baseline']:.6g}  "
              f"rel={r['rel_diff']:.2e} (tol {r['tolerance']:.0e})")
    print("REGRESSION", "PASSED" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, root)
    sys.exit(main())
