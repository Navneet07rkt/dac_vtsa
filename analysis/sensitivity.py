"""Local (one-at-a-time) sensitivity analysis of the cycle KPIs.

For every parameter p_i (dotted config key) the cycle is run to (near) cyclic
steady state at p_i (1 + d) and p_i (1 - d) -- or p_i +/- delta_abs for
temperatures / humidities -- and the normalised (logarithmic) sensitivity

        S_ij = d ln K_j / d ln p_i  ~=  (K_j+ - K_j-) / (2 K_j0) * p_i0 / (p_i+ - p_i0)

is reported for every KPI K_j (recovery, purity, productivity, kJ/kg ...).
|S| = 1 means a 10 % change of the parameter moves the KPI by ~10 %.

Cost control
------------
* runs use a reduced grid (sensitivity.grid, default axisymmetric 6x1x10) and
  at most sensitivity.cycles cycles;
* every perturbed run is warm-started from the base case's cyclic steady state
  (much faster convergence);
* runs are distributed over sensitivity.processes worker processes.

Interpretation caveat: local derivatives around ONE operating point; strongly
non-linear responses (e.g. T_heat near the regeneration threshold) need a scan,
not a derivative.
"""
from __future__ import annotations

import copy
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor
from typing import Dict, List, Tuple

import numpy as np

ABS_DELTA = {"operation.T_heat": 5.0, "operation.T_cool": 5.0, "operation.T_adsorption_jacket": 2.0,
             "feed.T": 2.0, "feed.RH": 0.05, "ambient.RH": 0.05}


def _as_dict(node) -> Dict:
    if isinstance(node, dict):
        return dict(node)
    if hasattr(node, "to_dict"):
        return node.to_dict()
    return dict(node._data)


def _parse_param(entry) -> Tuple[str, Dict]:
    if isinstance(entry, str):
        return entry, {}
    if isinstance(entry, dict):
        key = entry.get("key")
        return key, {k: v for k, v in entry.items() if k != "key"}
    raise ValueError(entry)


def _run_case(args):
    """Worker: run one configuration to CSS and return its final KPIs."""
    import sys
    root, overrides, run_dir, label = args
    if root not in sys.path:
        sys.path.insert(0, root)
    from config.parameters import load_config
    from vtsa_cycle import Logger, run_to_css
    cfg = load_config(None, None, overrides)
    t0 = time.time()
    log = Logger(os.path.join(run_dir, "log.txt"), echo=False)
    res = run_to_css(cfg, run_dir, log)
    k = dict(res.cycles[-1].kpis)
    k["_cycles"] = len(res.cycles)
    k["_css"] = res.converged
    k["_wall_s"] = time.time() - t0
    final_csv = os.path.join(run_dir, f"cycle_{len(res.cycles):03d}",
                             f"state_after_{cfg.cycle[-1]['name']}.csv")
    return label, k, final_csv


def run_sensitivity(cfg, out_dir: str, log=print) -> List[Dict]:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sc = cfg.sensitivity
    base_over = {"run.max_cycles": int(sc.cycles), "run.min_cycles": min(2, int(sc.cycles)),
                 "run.make_plots": False, "run.snapshots_per_step": 3}
    for k, v in _as_dict(sc.grid).items():
        base_over[f"grid.{k}"] = v
    raw = cfg.to_dict()
    raw.pop("derived", None)

    def flat_get(key):
        node = raw
        for p in key.split("."):
            node = node[p]
        return node

    # ---------------- base case
    os.makedirs(out_dir, exist_ok=True)
    log(f"sensitivity: base case on grid {_as_dict(sc.grid)} (<= {sc.cycles} cycles)")
    base_dir = os.path.join(out_dir, "base")
    label, k0, base_state = _run_case((root, _merge(raw, base_over), base_dir, "base"))
    log(f"  base: recovery {k0['recovery_pct']:.2f} %  thermal {k0['specific_thermal_kJ_per_kg']:.0f} kJ/kg  "
        f"({k0['_cycles']} cycles, {k0['_wall_s']:.0f} s)")
    warm = dict(base_over)
    warm.update({"run.initial_state": "csv", "run.initial_state_csv": base_state,
                 "run.max_cycles": max(2, int(sc.cycles) // 2 + 1)})
    # ---------------- perturbed cases
    jobs, meta = [], []
    for entry in sc.parameters:
        key, opts = _parse_param(entry)
        p0 = float(flat_get(key))
        if "delta_abs" in opts or key in ABS_DELTA:
            d = float(opts.get("delta_abs", ABS_DELTA.get(key)))
            vals = (p0 + d, p0 - d)
        else:
            d = float(opts.get("relative_step", sc.relative_step))
            vals = (p0 * (1 + d), p0 * (1 - d))
        for tag, v in zip(("plus", "minus"), vals):
            over = dict(warm)
            over[key] = v
            rd = os.path.join(out_dir, f"{key.replace('.', '_')}_{tag}")
            jobs.append((root, _merge(raw, over), rd, f"{key}|{tag}"))
        meta.append((key, p0, vals))
    results = {}
    with ProcessPoolExecutor(max_workers=int(sc.processes)) as ex:
        for label, k, _ in ex.map(_run_case, jobs):
            results[label] = k
            log(f"  done {label:<40s} recovery {k['recovery_pct']:.2f} %  ({k['_wall_s']:.0f} s)")
    # ---------------- sensitivities
    rows = []
    for key, p0, (pp, pm) in meta:
        kp, km = results[f"{key}|plus"], results[f"{key}|minus"]
        for kpi in sc.kpis:
            K0, Kp, Km = k0[kpi], kp[kpi], km[kpi]
            dlnp = (pp - pm) / (2.0 * p0) if p0 != 0 else np.nan
            S = (Kp - Km) / (2.0 * K0) / dlnp if (K0 not in (0, None) and dlnp) else np.nan
            rows.append({"parameter": key, "p_base": p0, "p_plus": pp, "p_minus": pm, "kpi": kpi,
                         "K_base": K0, "K_plus": Kp, "K_minus": Km,
                         "normalized_sensitivity": S,
                         "change_pct_plus": 100 * (Kp - K0) / K0 if K0 else np.nan,
                         "change_pct_minus": 100 * (Km - K0) / K0 if K0 else np.nan})
    return rows


def _merge(raw: Dict, over: Dict) -> Dict:
    """Return a flat override dict = every leaf of raw (so the worker sees the same base) + over."""
    flat = {}

    def walk(node, prefix):
        for k, v in node.items():
            key = f"{prefix}.{k}" if prefix else k
            if isinstance(v, dict) and key not in ("sensitivity.grid",):
                walk(v, key)
            else:
                flat[key] = v
    walk({k: v for k, v in raw.items() if k != "cycle"}, "")
    flat["cycle"] = raw["cycle"]
    flat.update(over)
    return flat
