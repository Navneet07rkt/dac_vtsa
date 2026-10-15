#!/usr/bin/env python3
"""Compare the Young-2021 1D-equivalent case with the published optimum (ESI Table S10).

    python run_dac_simulation.py --config config/cases/young2021_1d_equivalent.yaml
    python -m validation.young2021_compare [results/young2021_1d_equivalent]

Published values: Young et al. (2021) Energy Environ. Sci. 14, 5377, ESI Table S10, WADST
column, minimum-work-equivalent operating point (P_vac 0.2 bar, T_heat 99.73 C, 98.95 %
end-of-column saturation, t_heat 856 s, t_des 23033 s).  Writes young2021_comparison.csv
next to the run's kpis.json.

What agreement to expect: the case reproduces Young's 1D assumptions (one radial cell,
lumped jacket U, no wall heat capacity, their dispersion correlation) but NOT their exact
numerics or every closure (e.g. our pressure-driven blowdown/pressurisation, finite pump
and valve speeds, conservative energy balance).  Differences of a few percent are expected;
tens of percent would indicate a problem.
"""
from __future__ import annotations

import csv
import json
import os
import sys

PUBLISHED = [
    # key in kpis.json, label, published value, factor (ours -> published units), unit
    ("purity_dry_pct", "CO2 purity (dry)", 97.45, 1.0, "%"),
    ("working_capacity_mol_kg", "working capacity", 1.00, 1.0, "mol/kg"),
    ("specific_thermal_kJ_per_kg", "specific heat duty", 9.19, 1e-3, "MJ/kg CO2"),
    ("specific_electric_kJ_per_kg", "mechanical work (vacuum + blower)", 0.80, 1e-3, "MJ/kg CO2"),
    ("work_equivalent_kJ_per_kg", "work equivalent", 2.37, 1e-3, "MJ/kg CO2"),
    ("productivity_kg_m3_day", "productivity", 64.40, 1.0, "kg/m3/day"),
]


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    run_dir = argv[0] if argv else os.path.join(root, "results", "young2021_1d_equivalent")
    k = json.load(open(os.path.join(run_dir, "kpis.json")))
    rows = []
    print(f"Young et al. 2021 (ESI Table S10, WADST optimum) vs this model  [{run_dir}]")
    print(f"  cycles run: {k.get('cycles_run')}, CSS reached: {k.get('css_reached')}")
    print(f"  {'quantity':38s} {'this model':>12s} {'Young 2021':>12s} {'difference':>11s}")
    for key, label, pub, fac, unit in PUBLISHED:
        ours = float(k[key]) * fac
        diff = 100.0 * (ours - pub) / pub
        rows.append({"quantity": label, "unit": unit, "this_model": ours, "young2021": pub, "difference_pct": diff})
        print(f"  {label + ' [' + unit + ']':38s} {ours:12.3f} {pub:12.3f} {diff:+10.1f}%")
    print(f"  {'CO2 recovery [%] (not tabulated)':38s} {float(k['recovery_pct']):12.2f}")
    out = os.path.join(run_dir, "young2021_comparison.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("  ->", out)
    return rows


if __name__ == "__main__":
    main()
