"""Explicit RK4 stepper (quick standalone checks / non-stiff sub-steps).

Kept from the original scripts for transparency and for the verification test
that the explicit and implicit paths agree on a short horizon.  It is NOT
usable for full cycle steps: the pressure field makes the semi-discrete system
extremely stiff (see numerics/cfl.py) -- the integrator refuses to start if
the required number of steps exceeds `max_steps`.

Accumulators are integrated with the same RK4 weights, so linear invariants
(species inventories + cumulative boundary flows) are preserved exactly.
"""
from __future__ import annotations

import numpy as np

from numerics.cfl import explicit_dt


def rk4_step(model, t, y, dt):
    k1, a1 = model.rhs(t, y)
    k2, a2 = model.rhs(t + 0.5 * dt, y + 0.5 * dt * k1)
    k3, a3 = model.rhs(t + 0.5 * dt, y + 0.5 * dt * k2)
    k4, a4 = model.rhs(t + dt, y + dt * k3)
    return (y + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4),
            dt / 6.0 * (a1 + 2 * a2 + 2 * a3 + a4))


def integrate_rk4(model, y0, t_end, dt=None, max_steps=200000, safety=0.5):
    y = y0.copy()
    acc = np.zeros(model.n_acc)
    if dt is None:
        dt = explicit_dt(model, y, safety)
    n = int(np.ceil(t_end / dt))
    if n > max_steps:
        raise RuntimeError(f"RK4 would need {n:.3g} steps (dt = {dt:.3e} s) for {t_end} s; "
                           "use solver.method = trbdf2 (implicit) for full steps")
    dt = t_end / n
    t = 0.0
    for _ in range(n):
        y, da = rk4_step(model, t, y, dt)
        acc += da
        t += dt
    return y, acc, t, n
