"""Explicit stability limits (CFL / diffusion numbers) for the RK4 path.

For an explicit method the admissible step is bounded by the fastest
characteristic time of the semi-discrete system:

    convection        dt <= C  dx / v_int
    dispersion        dt <= 0.5 dx^2 / D
    conduction        dt <= 0.5 dx^2 rho_c / lam
    PRESSURE          dt <= 0.5 dx^2 / D_P,   D_P = K P / (eps_t mu),  K = eps^3 d_p^2 / (150 (1-eps)^2)
    LDF               dt <= 1 / k_LDF

The pressure-diffusion limit is the killer: for the default grid it is
~1e-7 s, so an explicit 8000 s adsorption step needs ~1e11 RK4 steps.  That
is why the production path is implicit (core/solver.py).  Use
:func:`explicit_dt` to see the numbers for your grid.
"""
from __future__ import annotations

from typing import Dict

import numpy as np

from properties import R_GAS


def explicit_limits(model, y) -> Dict[str, float]:
    g = model.grid
    cfg = model.cfg
    d = model.rhs(0.0, y, diag=True)
    T, P = d["T"], d["P"]
    mu = 1.8e-5
    eps = model.eps_b
    dp = cfg.sorbent.particle_diameter
    dx_min = min(g.dz, float(np.min(g.dr)), float(np.min(g.r_c * g.dth)) if g.Nth > 1 else np.inf)
    v = float(np.max(np.abs(d["uz"]))) / eps if d["uz"].size else 0.0
    v = max(v, cfg.feed.velocity / eps)
    K = eps ** 3 * dp ** 2 / (150.0 * (1.0 - eps) ** 2)
    D_P = K * float(np.max(P)) / (model.eps_t * mu)
    rhoc = model.rho_b * cfg.sorbent.cp_solid
    lam = 0.3
    Dm = 1e-4
    out = {
        "convection": 0.9 * dx_min / max(v, 1e-12),
        "dispersion": 0.5 * dx_min ** 2 / Dm,
        "conduction": 0.5 * dx_min ** 2 * rhoc / lam,
        "pressure_diffusion": 0.5 * dx_min ** 2 / D_P,
        "ldf": 1.0 / max(cfg.kinetics.k_ldf_co2, cfg.kinetics.k_ldf_h2o),
    }
    out["min"] = min(out.values())
    return out


def explicit_dt(model, y, safety: float = 0.5) -> float:
    return safety * explicit_limits(model, y)["min"]


def check_dt(model, y, dt: float, raise_on_violation: bool = False) -> bool:
    lim = explicit_limits(model, y)
    ok = dt <= lim["min"]
    if not ok:
        msg = (f"explicit dt = {dt:.3e} s exceeds the stability limit {lim['min']:.3e} s "
               f"(limiting process: {min((v, k) for k, v in lim.items() if k != 'min')[1]})")
        if raise_on_violation:
            raise ValueError(msg)
        print("WARNING:", msg)
    return ok
