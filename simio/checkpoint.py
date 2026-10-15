"""Save / load a full CycleState as CSV (the step-to-step hand-off files).

Two files per state:
    <name>.csv        one row per bed cell:  i_r, i_theta, i_z, r, theta_deg, z,
                      T, P, c_<gas>..., y_<gas>..., q_<ads>...
    <name>_wall.csv   one row per wall cell: i_theta, i_z, theta_deg, z, T_wall
plus a one-line header comment with t, cycle, step.

Floats are written with 17 significant digits, so a save/load round trip is
lossless (bit-exact): a step started from the CSV produces exactly the same
result as one handed the in-memory state.
"""
from __future__ import annotations

import csv
import os
from typing import Sequence

import numpy as np

from core.state import CycleState
from properties import R_GAS

FMT = "%.17g"


def save_state_csv(state: CycleState, grid, path: str) -> str:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    r, th, z = grid.cell_centres()
    P = state.P
    y = state.y
    cols = ["i_r", "i_theta", "i_z", "r_m", "theta_deg", "z_m", "T_K", "P_Pa"]
    cols += [f"c_{s}_mol_m3" for s in state.species]
    cols += [f"y_{s}" for s in state.species]
    cols += [f"q_{s}_mol_kg" for s in state.adsorbed]
    ii, jj, kk = np.meshgrid(np.arange(grid.Nr), np.arange(grid.Nth), np.arange(grid.Nz), indexing="ij")
    data = [ii.ravel(), jj.ravel(), kk.ravel(), r.ravel(), np.degrees(th).ravel(), z.ravel(),
            state.T.ravel(), P.ravel()]
    data += [state.c[k].ravel() for k in range(len(state.species))]
    data += [y[k].ravel() for k in range(len(state.species))]
    data += [state.q[k].ravel() for k in range(len(state.adsorbed))]
    arr = np.column_stack(data)
    header = (f"# t_s={float(state.t)!r} cycle={int(state.cycle)} step={state.step} "
              f"grid={grid.Nr}x{grid.Nth}x{grid.Nz}\n" + ",".join(cols))
    fmt = ["%d", "%d", "%d"] + [FMT] * (arr.shape[1] - 3)
    np.savetxt(path, arr, delimiter=",", header=header, comments="", fmt=fmt)
    wpath = path[:-4] + "_wall.csv"
    jw, kw = np.meshgrid(np.arange(grid.Nth), np.arange(grid.Nz), indexing="ij")
    warr = np.column_stack([jw.ravel(), kw.ravel(), np.degrees(grid.th_c)[jw.ravel()],
                            grid.z_c[kw.ravel()], state.Tw.ravel()])
    np.savetxt(wpath, warr, delimiter=",", header="i_theta,i_z,theta_deg,z_m,T_wall_K", comments="",
               fmt=["%d", "%d", FMT, FMT, FMT])
    return path


def load_state_csv(path: str, grid, species: Sequence[str], adsorbed: Sequence[str]) -> CycleState:
    with open(path) as fh:
        first = fh.readline().strip()
        header = fh.readline().strip().split(",")
    meta = {}
    for tok in first.lstrip("# ").split():
        if "=" in tok:
            k, v = tok.split("=", 1)
            meta[k] = v
    arr = np.loadtxt(path, delimiter=",", skiprows=2, ndmin=2)
    col = {h: i for i, h in enumerate(header)}
    expect = grid.Nr * grid.Nth * grid.Nz
    if arr.shape[0] != expect:
        raise ValueError(f"{path}: has {arr.shape[0]} cells but the grid has {expect} "
                         f"({grid.Nr}x{grid.Nth}x{grid.Nz}); use the same grid as the run that wrote it")
    i = arr[:, col["i_r"]].astype(int)
    j = arr[:, col["i_theta"]].astype(int)
    k = arr[:, col["i_z"]].astype(int)
    shape = grid.shape
    T = np.empty(shape)
    T[i, j, k] = arr[:, col["T_K"]]
    c = np.empty((len(species),) + shape)
    for n, s in enumerate(species):
        name = f"c_{s}_mol_m3"
        if name not in col:
            raise ValueError(f"{path}: missing column {name} (was the run dry/humid?)")
        c[n][i, j, k] = arr[:, col[name]]
    q = np.empty((len(adsorbed),) + shape)
    for n, s in enumerate(adsorbed):
        q[n][i, j, k] = arr[:, col[f"q_{s}_mol_kg"]]
    wpath = path[:-4] + "_wall.csv"
    Tw = np.empty(grid.wall_shape)
    if os.path.exists(wpath):
        w = np.loadtxt(wpath, delimiter=",", skiprows=1, ndmin=2)
        Tw[w[:, 0].astype(int), w[:, 1].astype(int)] = w[:, 4]
    else:
        Tw[:] = T[-1]
    st = CycleState(c=c, q=q, T=T, Tw=Tw, t=float(meta.get("t_s", 0.0)),
                    cycle=int(meta.get("cycle", 1)), step=meta.get("step", "loaded"),
                    species=tuple(species), adsorbed=tuple(adsorbed))
    return st
