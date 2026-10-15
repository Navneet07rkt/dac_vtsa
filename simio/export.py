"""Raw data export (CSV), independent of any plotting.

* :func:`write_rows`        list of dicts -> CSV (outlet time series, KPI tables, balances)
* :func:`fields_long_rows`  snapshots -> long-format rows  (cycle, step, t, r, theta, z, q, T, P, ...)
  so that  q(r, theta, z, t),  T(r, theta, z, t),  P(r, theta, z, t)  are available for every cell
* :func:`write_npz`         compact binary copy of the same snapshots (fast reload for probing)
"""
from __future__ import annotations

import csv
import os
from typing import Dict, Iterable, List, Sequence

import numpy as np


def write_rows(path: str, rows: Sequence[Dict], columns: Sequence[str] = None, append: bool = False) -> str:
    if not rows:
        return path
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    cols = list(columns) if columns else list(rows[0].keys())
    new = not (append and os.path.exists(path))
    with open(path, "a" if append else "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        if new:
            w.writeheader()
        for r in rows:
            w.writerow({k: _fmt(r.get(k, "")) for k in cols})
    return path


def _fmt(v):
    if isinstance(v, (float, np.floating)):
        return f"{float(v):.10g}"
    if isinstance(v, (np.integer,)):
        return int(v)
    return v


FIELD_COLUMNS_BASE = ["cycle", "step", "t_step_s", "t_cycle_s", "t_total_s", "i_r", "i_theta", "i_z",
                      "r_m", "theta_deg", "z_m", "T_K", "P_Pa"]


def field_columns(species, adsorbed):
    return (FIELD_COLUMNS_BASE + [f"c_{s}_mol_m3" for s in species] + [f"y_{s}" for s in species]
            + [f"q_{s}_mol_kg" for s in adsorbed])


def write_fields_csv(path: str, grid, snaps: List[Dict], species, adsorbed, append: bool = False) -> str:
    """snaps: list of dicts {cycle, step, t_step, t_cycle, t_total, state}  -> long CSV."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    cols = field_columns(species, adsorbed)
    r, th, z = grid.cell_centres()
    ii, jj, kk = np.meshgrid(np.arange(grid.Nr), np.arange(grid.Nth), np.arange(grid.Nz), indexing="ij")
    geo = [ii.ravel(), jj.ravel(), kk.ravel(), r.ravel(), np.degrees(th).ravel(), z.ravel()]
    new = not (append and os.path.exists(path))
    with open(path, "a" if append else "w") as fh:
        if new:
            fh.write(",".join(cols) + "\n")
        for s in snaps:
            st = s["state"]
            n = grid.ncell
            y = st.y
            data = [np.full(n, s["t_step"]), np.full(n, s["t_cycle"]), np.full(n, s["t_total"])] + geo
            data += [st.T.ravel(), st.P.ravel()]
            data += [st.c[k].ravel() for k in range(len(species))]
            data += [y[k].ravel() for k in range(len(species))]
            data += [st.q[k].ravel() for k in range(len(adsorbed))]
            arr = np.column_stack(data)
            prefix = f"{s['cycle']},{s['step']},"
            lines = []
            for row in arr:
                lines.append(prefix + "%.9g,%.9g,%.9g,%d,%d,%d," % tuple(row[:6])
                             + ",".join("%.9g" % v for v in row[6:]))
            fh.write("\n".join(lines) + "\n")
    return path


def write_npz(path: str, grid, snaps: List[Dict]) -> str:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    if not snaps:
        return path
    st0 = snaps[0]["state"]
    np.savez_compressed(
        path,
        r=grid.r_c, theta=grid.th_c, z=grid.z_c, R=grid.R, L=grid.L,
        step=np.array([s["step"] for s in snaps]), cycle=np.array([s["cycle"] for s in snaps]),
        t_step=np.array([s["t_step"] for s in snaps]), t_cycle=np.array([s["t_cycle"] for s in snaps]),
        t_total=np.array([s["t_total"] for s in snaps]),
        T=np.array([s["state"].T for s in snaps]), P=np.array([s["state"].P for s in snaps]),
        c=np.array([s["state"].c for s in snaps]), q=np.array([s["state"].q for s in snaps]),
        Tw=np.array([s["state"].Tw for s in snaps]),
        species=np.array(st0.species), adsorbed=np.array(st0.adsorbed),
    )
    return path
