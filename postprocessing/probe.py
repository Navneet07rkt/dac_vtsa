"""Evaluate any field at any (r, theta, z, t) of a finished run.

The solver stores cell averages on the (Nr, Ntheta, Nz) grid at the snapshot
times of every step.  `FieldProbe` interpolates them:

* space: linear in r between ring centres (theta-averaged value on the axis
  r = 0, constant extrapolation to the wall r = R), periodic-linear in theta,
  linear in z (constant to the end faces);
* time: linear between the two snapshots that bracket the requested time.

Fields: T, P, c_<gas>, y_<gas>, q_<ads>, and derived p_<gas> (partial pressure),
RH (relative humidity), q_star_CO2 / q_star_H2O (local equilibrium loading,
needs the run's config_used.yaml), Tw (wall, function of theta and z).

Python
------
>>> from postprocessing.probe import FieldProbe
>>> pr = FieldProbe("results/default_3d")
>>> pr.value("q_CO2", r=0.02, theta_deg=90, z=0.005, t_cycle=3600)
>>> pr.series("T", r=0.0, theta_deg=0, z=0.005)          # (t_cycle, values) over the final cycle

CLI
---
python -m postprocessing.probe results/default_3d --field q_CO2 --r 0.02 --theta 90 --z 0.005 --t-cycle 3600
python -m postprocessing.probe results/default_3d --field T --r 0 --z 0.005 --series
"""
from __future__ import annotations

import argparse
import os
import sys
from typing import Optional

import numpy as np

R_GAS = 8.314462618


class FieldProbe:
    def __init__(self, run_dir: str, npz_name: str = "fields_final_cycle.npz"):
        path = os.path.join(run_dir, npz_name)
        if not os.path.exists(path):
            raise FileNotFoundError(f"{path} not found -- run run_dac_simulation.py first")
        d = np.load(path, allow_pickle=False)
        self.run_dir = run_dir
        self.r, self.theta, self.z = d["r"], d["theta"], d["z"]
        self.R, self.L = float(d["R"]), float(d["L"])
        self.step = d["step"]
        self.t_cycle = d["t_cycle"]
        self.t_step = d["t_step"]
        self.T, self.P, self.c, self.q, self.Tw = d["T"], d["P"], d["c"], d["q"], d["Tw"]
        self.species = [str(s) for s in d["species"]]
        self.adsorbed = [str(s) for s in d["adsorbed"]]
        order = np.argsort(self.t_cycle, kind="stable")
        self._order = order
        self._eq = None

    # ------------------------------------------------------------------ raw fields per snapshot
    def _snapshot_field(self, name: str, k: int) -> np.ndarray:
        if name == "T":
            return self.T[k]
        if name == "P":
            return self.P[k]
        if name.startswith("c_"):
            return self.c[k][self.species.index(name[2:])]
        if name.startswith("y_"):
            C = self.c[k].sum(axis=0)
            return self.c[k][self.species.index(name[2:])] / np.maximum(C, 1e-30)
        if name.startswith("p_"):
            return self.c[k][self.species.index(name[2:])] * R_GAS * self.T[k]
        if name.startswith("q_star_"):
            return self._qstar(k)[self.adsorbed.index(name[7:])]
        if name.startswith("q_"):
            return self.q[k][self.adsorbed.index(name[2:])]
        if name == "RH":
            from properties.ideal_gas_backend import arden_buck
            return self.c[k][self.species.index("H2O")] * R_GAS * self.T[k] / arden_buck(self.T[k])
        raise KeyError(f"unknown field '{name}'")

    def _qstar(self, k):
        if self._eq is None:
            root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            if root not in sys.path:
                sys.path.insert(0, root)
            from config.parameters import load_config
            from physics.isotherms import EquilibriumModel
            from properties import get_backend
            cfg = load_config(os.path.join(self.run_dir, "config_used.yaml"))
            self._eq = (EquilibriumModel(cfg), get_backend(cfg))
        eq, props = self._eq
        T = self.T[k]
        p = self.c[k] * R_GAS * T
        ic = self.species.index("CO2")
        if "H2O" in self.species:
            iw = self.species.index("H2O")
            return eq.equilibrium(p[ic], p[iw], T, self.q[k][1], props.p_sat(T))
        return eq.equilibrium(p[ic], None, T, None, props.p_sat(T))

    # ------------------------------------------------------------------ interpolation
    def _z_theta(self, F: np.ndarray, th, z):
        """Interpolate the last two axes (theta, z) of F; returns F[..., ] at (th, z)."""
        zi = np.clip(np.interp(z, self.z, np.arange(len(self.z))), 0, len(self.z) - 1)
        k0 = int(np.floor(zi)); k1 = min(k0 + 1, len(self.z) - 1); wz = zi - k0
        Fz = (1 - wz) * F[..., k0] + wz * F[..., k1]                       # (..., Nth)
        nth = len(self.theta)
        if nth == 1:
            return Fz[..., 0], Fz
        dth = 2 * np.pi / nth
        x = ((th - self.theta[0]) % (2 * np.pi)) / dth
        j0 = int(np.floor(x)) % nth; j1 = (j0 + 1) % nth; wt = x - np.floor(x)
        return (1 - wt) * Fz[..., j0] + wt * Fz[..., j1], Fz

    def _space(self, F: np.ndarray, r, th, z) -> float:
        Ft, Fz = self._z_theta(F, th, z)                                    # (Nr,), (Nr, Nth)
        # r: axis value = theta-average of the innermost ring
        axis = Fz[0].mean()
        rr = np.concatenate([[0.0], self.r])
        vals = np.concatenate([[axis], Ft])
        return float(np.interp(min(max(r, 0.0), self.R), rr, vals))

    def _time_index(self, t_cycle: Optional[float], step: Optional[str], t_step: Optional[float]):
        if step is not None:
            idx = np.where(self.step == step)[0]
            if len(idx) == 0:
                raise KeyError(f"no snapshots for step '{step}' (available: {sorted(set(self.step))})")
            ts = self.t_step[idx]
            t = ts[-1] if t_step is None else t_step
            o = idx[np.argsort(ts)]
            tt = np.sort(ts)
        else:
            o = self._order
            tt = self.t_cycle[o]
            t = tt[-1] if t_cycle is None else t_cycle
        if t <= tt[0]:
            return o[0], o[0], 0.0
        if t >= tt[-1]:
            return o[-1], o[-1], 0.0
        i = int(np.searchsorted(tt, t)) - 1
        w = (t - tt[i]) / max(tt[i + 1] - tt[i], 1e-30)
        return o[i], o[i + 1], w

    def value(self, field: str, r: float, theta_deg: float = 0.0, z: float = 0.0,
              t_cycle: Optional[float] = None, step: Optional[str] = None, t_step: Optional[float] = None):
        th = np.radians(theta_deg)
        if field == "Tw":                       # wall temperature: function of (theta, z) only
            a, b, w = self._time_index(t_cycle, step, t_step)
            fa = float(self._z_theta(self.Tw[a], th, z)[0])
            fb = float(self._z_theta(self.Tw[b], th, z)[0])
            return (1 - w) * fa + w * fb
        a, b, w = self._time_index(t_cycle, step, t_step)
        fa = self._space(self._snapshot_field(field, a), r, th, z)
        if b == a:
            return fa
        fb = self._space(self._snapshot_field(field, b), r, th, z)
        return (1 - w) * fa + w * fb

    def series(self, field: str, r: float, theta_deg: float = 0.0, z: float = 0.0):
        o = self._order
        th = np.radians(theta_deg)
        vals = np.array([self._space(self._snapshot_field(field, k), r, th, z) for k in o])
        return self.t_cycle[o], vals, self.step[o]


def main(argv=None):
    ap = argparse.ArgumentParser(description="probe a field at (r, theta, z, t) of a finished run")
    ap.add_argument("run_dir")
    ap.add_argument("--field", default="q_CO2")
    ap.add_argument("--r", type=float, default=0.0, help="radius [m]")
    ap.add_argument("--theta", type=float, default=0.0, help="angle [deg]")
    ap.add_argument("--z", type=float, default=0.0, help="height from the top [m]")
    ap.add_argument("--t-cycle", type=float, default=None, help="time in the cycle [s]")
    ap.add_argument("--step", default=None)
    ap.add_argument("--t-step", type=float, default=None, help="time within --step [s]")
    ap.add_argument("--series", action="store_true", help="print the whole cycle history at the point")
    a = ap.parse_args(argv)
    pr = FieldProbe(a.run_dir)
    if a.series:
        t, v, s = pr.series(a.field, a.r, a.theta, a.z)
        print(f"t_cycle_s,step,{a.field}")
        for ti, si, vi in zip(t, s, v):
            print(f"{ti:.3f},{si},{vi:.8g}")
    else:
        v = pr.value(a.field, a.r, a.theta, a.z, a.t_cycle, a.step, a.t_step)
        print(f"{a.field}(r={a.r}, theta={a.theta} deg, z={a.z}, t={a.t_cycle if a.step is None else (a.step, a.t_step)}) = {v:.8g}")


if __name__ == "__main__":
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, root)
    main()
