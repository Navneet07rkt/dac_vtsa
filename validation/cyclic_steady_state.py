"""Cyclic steady state (CSS) detection.

Two criteria, both must hold (and run.min_cycles must have been reached):

1. state change   max over fields of  ||X_N - X_{N-1}|| / scale   < run.css_tolerance
   (volume-weighted RMS; q scaled by its RMS value, T by the temperature swing
    T_heat - T_feed, gas concentrations by the RMS total concentration);
2. cycle CO2 balance  |fed - product - vented| / fed  < run.css_mb_tolerance
   (at CSS nothing accumulates over a cycle; Young et al. 2021 used 0.5 %).
"""
from __future__ import annotations

from typing import Dict

import numpy as np


def _rms(a, w):
    return float(np.sqrt(np.sum(w * a * a) / np.sum(w)))


def state_change(prev, now, grid, cfg) -> Dict[str, float]:
    w = grid.volume_full()
    out = {}
    for n, s in enumerate(now.adsorbed):
        sc = max(_rms(now.q[n], w), 1e-6)
        out[f"q_{s}"] = _rms(now.q[n] - prev.q[n], w) / sc
    dT = max(cfg.operation.T_heat - cfg.feed.T, 1.0)
    out["T"] = _rms(now.T - prev.T, w) / dT
    Csc = max(_rms(now.c.sum(axis=0), w), 1e-9)
    for k, s in enumerate(now.species):
        out[f"c_{s}"] = _rms(now.c[k] - prev.c[k], w) / Csc
    out["max"] = max(out.values())
    return out


def cycle_co2_balance(kpis: Dict[str, float]) -> float:
    fed = kpis["CO2_fed_mol"]
    return abs(fed - kpis["CO2_product_mol"] - kpis["CO2_vented_mol"]) / fed if fed > 0 else np.inf


def is_css(change: Dict[str, float], kpis: Dict[str, float], cfg, cycle: int) -> bool:
    return (cycle >= cfg.run.min_cycles and change["max"] < cfg.run.css_tolerance
            and cycle_co2_balance(kpis) < cfg.run.css_mb_tolerance)
