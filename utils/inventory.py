"""Shared inventory bookkeeping (one implementation used by every step).

compute_inventory(state)      -> moles of every species in gas and adsorbed phase [mol]
compute_co2_inventory(state)  -> total CO2 held by the bed (gas + adsorbed) [mol]
"""
from __future__ import annotations

from typing import Dict

import numpy as np


def compute_inventory(state, grid, cfg) -> Dict[str, float]:
    eps_t = cfg.sorbent.total_voidage
    rho_b = cfg.derived.bulk_density
    out = {}
    for k, s in enumerate(state.species):
        gas = float(np.sum(eps_t * state.c[k] * grid.V))
        ads = 0.0
        if s in state.adsorbed:
            ads = float(np.sum(rho_b * state.q[state.adsorbed.index(s)] * grid.V))
        out[f"gas_{s}"] = gas
        out[f"ads_{s}"] = ads
        out[f"tot_{s}"] = gas + ads
    return out


def compute_co2_inventory(state, grid, cfg) -> float:
    return compute_inventory(state, grid, cfg)["tot_CO2"]


def bed_averages(state, grid, cfg) -> Dict[str, float]:
    V = grid.volume_full()
    Vt = float(V.sum())
    out = {"T_mean_K": float(np.sum(state.T * V) / Vt),
           "T_min_K": float(state.T.min()), "T_max_K": float(state.T.max()),
           "P_mean_Pa": float(np.sum(state.P * V) / Vt)}
    for n, s in enumerate(state.adsorbed):
        out[f"q_{s}_mean"] = float(np.sum(state.q[n] * V) / Vt)
    return out
