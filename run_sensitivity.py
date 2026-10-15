#!/usr/bin/env python3
"""Local sensitivity analysis of the TVSA KPIs (see analysis/sensitivity.py).

    python run_sensitivity.py                                  # parameters/KPIs from defaults.yaml
    python run_sensitivity.py --config my.yaml --set sensitivity.processes=4

Writes results/<run name>_sensitivity/: sensitivity.csv (normalised sensitivities
d ln KPI / d ln p), sensitivity_summary.txt, tornado_<kpi>.png
"""
from __future__ import annotations

import argparse
import os
import sys
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from analysis.sensitivity import run_sensitivity  # noqa: E402
from config.parameters import load_config  # noqa: E402
from simio.export import write_rows  # noqa: E402
from vtsa_cycle import Logger  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--set", action="append", default=[])
    ap.add_argument("--name", default=None)
    a = ap.parse_args(argv)
    cfg = load_config(a.config, a.set)
    name = a.name or f"{cfg.run.name}_sensitivity"
    out = os.path.join(ROOT, cfg.run.results_dir, name)
    os.makedirs(out, exist_ok=True)
    log = Logger(os.path.join(out, "sensitivity.log"))
    t0 = time.time()
    rows = run_sensitivity(cfg, out, log)
    write_rows(os.path.join(out, "sensitivity.csv"), rows)
    lines = ["normalised sensitivities  S = d ln(KPI) / d ln(parameter)   (|S| > 1: KPI moves more than the parameter)",
             ""]
    kpis = list(dict.fromkeys(r["kpi"] for r in rows))
    params = list(dict.fromkeys(r["parameter"] for r in rows))
    head = f"{'parameter':<34s}" + "".join(f"{k[:22]:>24s}" for k in kpis)
    lines.append(head)
    for p in params:
        vals = {r["kpi"]: r["normalized_sensitivity"] for r in rows if r["parameter"] == p}
        lines.append(f"{p:<34s}" + "".join(f"{vals[k]:>24.3f}" for k in kpis))
    txt = "\n".join(lines)
    with open(os.path.join(out, "sensitivity_summary.txt"), "w") as fh:
        fh.write(txt + "\n")
    log(txt)
    try:
        from postprocessing.plotting import plot_tornado
        for k in kpis:
            plot_tornado(rows, k, os.path.join(out, f"tornado_{k}.png"))
    except Exception as exc:   # plotting must never lose the numbers
        log(f"(plotting skipped: {exc})")
    log(f"sensitivity analysis finished in {time.time() - t0:.0f} s -> {out}")


if __name__ == "__main__":
    main()
