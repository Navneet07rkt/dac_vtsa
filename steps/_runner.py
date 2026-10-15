"""Shared machinery for every cycle step (physics is NOT duplicated per step).

run_step(state_in, ...) -> (state_out, StepResult)

1. builds the semi-discrete model for the step's valve states / jacket (core/residual.py)
2. integrates it with the configured solver (TR-BDF2 Newton default, scipy BDF, RK4)
3. records the outlet time series at every accepted step and field snapshots
4. returns the end state plus everything needed for balances and KPIs
"""
from __future__ import annotations

import argparse
import os
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

from core.residual import TVSAModel
from core.solver import TRBDF2
from properties import R_GAS
from utils.inventory import bed_averages


@dataclass
class StepResult:
    name: str
    type: str
    cycle: int
    t_start_cycle: float
    duration_planned: float
    duration: float
    stopped_early: bool
    acc: Dict[str, float]
    inv_start: Dict[str, float]
    inv_end: Dict[str, float]
    outlet: List[Dict] = field(default_factory=list)
    snapshots: List[Dict] = field(default_factory=list)
    solver: Dict = field(default_factory=dict)
    wall_time: float = 0.0
    T_jacket: float = 298.15
    jacket_on: bool = True
    product_port: str = "V3"
    valves: Dict = field(default_factory=dict)
    ads_start: Dict[str, float] = field(default_factory=dict)
    ads_end: Dict[str, float] = field(default_factory=dict)


def _end_summary(model, d, end, c_feed_co2):
    b = d["bnd"][end]
    Fout = -b["F_in"]                                # (ng, Nr, Nth) leaving the bed
    tot = Fout.sum(axis=0)
    Ftot = float(tot.sum())
    k_lay = 0 if end == "top" else -1
    T_lay = d["T"][..., k_lay]
    A = model.grid.Az * np.ones_like(T_lay)
    P_face = float(np.sum(b["P_face"] * A) / np.sum(A))
    out_faces = tot > 0.0                            # faces through which gas LEAVES the bed
    if "F_outpart" in b:                             # flux-split valve: the leaving stream itself
        Fout = b["F_outpart"]
        tot = Fout.sum(axis=0)
        out_faces = tot > 0.0
    F_leaving = float(tot[out_faces].sum()) if np.any(out_faces) else 0.0
    if F_leaving > 1e-14:
        # composition of the gas actually leaving (outflow faces only: inflow faces carry external
        # gas and must not enter the average -- that is what made C_out/C_feed blow up or go
        # negative while a vent switched from back-flow to out-flow)
        ys = [float(Fout[k][out_faces].sum()) / F_leaving for k in range(model.ng)]
        T_g = float(np.sum(T_lay[out_faces] * tot[out_faces]) / F_leaving)
    else:   # no out-flow: report the composition of the gas sitting at that end
        yl = d["y"][..., k_lay]
        ys = [float(np.sum(yl[k] * A) / np.sum(A)) for k in range(model.ng)]
        T_g = float(np.sum(T_lay * A) / np.sum(A))
    out = {f"{end}_valve": model.conds[end].valve or "closed",
           f"{end}_valve_state": model.conds[end].state.value,
           f"{end}_flow_out_mol_s": Ftot,
           f"{end}_P_Pa": P_face, f"{end}_T_K": T_g}
    for k, s in enumerate(model.species):
        out[f"{end}_y_{s}"] = ys[k]
    c_co2 = ys[model.i_co2] * P_face / (R_GAS * T_g)
    out[f"{end}_c_CO2_mol_m3"] = c_co2
    out[f"{end}_C_over_Cfeed"] = c_co2 / c_feed_co2 if c_feed_co2 > 0 else np.nan
    return out


def run_step(state_in, cfg, grid, props, step: Dict, cycle_idx: int = 1, t_cycle0: float = 0.0,
             T_jacket_prev: Optional[float] = None, logger=None, snapshots: Optional[int] = None):
    log = logger or (lambda *a, **k: None)
    model = TVSAModel(cfg, grid, props, step, state_in, T_jacket_prev)
    y0 = model.pack(state_in)
    dur = float(step["duration"])
    n_snap = snapshots if snapshots is not None else cfg.run.snapshots_per_step
    t_out = np.linspace(0.0, dur, max(2, n_snap))
    from valves.valve_logic import humid_air_composition
    yf = humid_air_composition(model.species, cfg.feed.P, cfg.feed.T,
                               cfg.feed.RH if model.water else 0.0, cfg.feed.y_co2_dry,
                               float(props.p_sat(cfg.feed.T)))
    c_feed_co2 = yf[model.i_co2] * cfg.feed.P / (R_GAS * cfg.feed.T)
    outlet: List[Dict] = []
    every = cfg.solver.outlet_record_every
    counter = [0]

    def record(t, y, acc):
        d = model.rhs(t, y, diag=True)
        row = {"cycle": cycle_idx, "step": step["name"], "t_step_s": t, "t_cycle_s": t_cycle0 + t,
               "t_total_s": state_in.t + t}
        for end in ("top", "bottom"):
            row.update(_end_summary(model, d, end, c_feed_co2))
        V = grid.volume_full()
        Vt = float(V.sum())
        row["bed_T_mean_K"] = float(np.sum(d["T"] * V) / Vt)
        row["bed_P_mean_Pa"] = float(np.sum(d["P"] * V) / Vt)
        c, q, e, Tw = model.unpack(y)
        for n, s in enumerate(model.adsorbed):
            row[f"bed_q_{s}_mean"] = float(np.sum(q[n] * V) / Vt)
        ai = model.acc_index
        for v in ("V1", "V2", "V3"):
            row[f"cum_CO2_out_{v}_mol"] = acc[ai[f"n_out_{v}_CO2"]]
        row["cum_Q_jacket_J"] = acc[ai["Q_jac_pos"]] + acc[ai["Q_jac_neg"]]
        outlet.append(row)

    last_print = [time.time()]

    def on_step(t, y, acc):
        counter[0] += 1
        if counter[0] % every == 0:
            record(t, y, acc)
        if time.time() - last_print[0] > 30.0:
            last_print[0] = time.time()
            log(f"      ... {step['name']}: t = {t:.0f}/{dur:.0f} s ({100 * t / dur:.0f} %), {counter[0]} steps")

    record(0.0, y0, np.zeros(model.n_acc))
    stop = step.get("stop_on")
    stop_fn = stop_val = None
    direction = 1
    if stop:
        stop_fn = model.stop_metric
        stop_val = stop["value"]
        direction = -1 if stop["criterion"] == "min_temperature" else 1
    wall0 = time.time()
    method = cfg.solver.method
    if method == "trbdf2":
        integ = TRBDF2(model, cfg, h_max=cfg.solver.h_max_fraction * dur)
        res = integ.integrate(y0, dur, t_out, on_step=on_step, stop_fn=stop_fn, stop_value=stop_val,
                              stop_direction=direction)
        y_end, acc_end, t_end = res.y, res.acc, res.t
        snaps_raw = res.snapshots
        solver_stats = {"method": "trbdf2", "steps": res.n_steps, "rejected": res.n_reject,
                        "newton_fail": res.n_newton_fail, "krylov": res.n_krylov,
                        "jacobians": res.n_jac, "lu": res.n_lu, "nfev": res.nfev,
                        "groups": integ.jac.ngroups, **integ.stats}
        stopped = res.stopped_early
    elif method == "scipy_bdf":
        y_end, acc_end, t_end, snaps_raw, solver_stats = _run_scipy(model, y0, dur, t_out, on_step)
        stopped = False
    else:
        y_end, acc_end, t_end, snaps_raw, solver_stats = _run_rk4(model, y0, dur, t_out, on_step)
        stopped = False
    if outlet[-1]["t_step_s"] < t_end:
        record(t_end, y_end, acc_end)
    wall = time.time() - wall0
    state_out = model.to_state(y_end, state_in, state_in.t + t_end, step["name"])
    state_out.cycle = cycle_idx
    snaps = []
    for (ts, ys, _a) in snaps_raw:
        st = model.to_state(ys, state_in, state_in.t + ts, step["name"])
        st.cycle = cycle_idx
        snaps.append({"cycle": cycle_idx, "step": step["name"], "t_step": ts, "t_cycle": t_cycle0 + ts,
                      "t_total": state_in.t + ts, "state": st})
    inv0 = model.inventories(y0)
    inv1 = model.inventories(y_end)
    acc_d = {nm: float(acc_end[i]) for i, nm in enumerate(model.acc_names)}
    T_jac = model.T_jacket(t_end) if model.jacket_on else cfg.ambient.T
    sr = StepResult(name=step["name"], type=step["type"], cycle=cycle_idx, t_start_cycle=t_cycle0,
                    duration_planned=dur, duration=t_end, stopped_early=stopped, acc=acc_d,
                    inv_start=inv0, inv_end=inv1, outlet=outlet, snapshots=snaps, solver=solver_stats,
                    wall_time=wall, T_jacket=float(T_jac), jacket_on=model.jacket_on,
                    product_port=step.get("product_port", "V3"), valves=dict(step["valves"]),
                    ads_start={s: inv0[f"ads_{s}"] for s in model.adsorbed},
                    ads_end={s: inv1[f"ads_{s}"] for s in model.adsorbed})
    ba = bed_averages(state_out, grid, cfg)
    log(f"    {step['name']:<14s} {t_end:9.1f} s  steps={solver_stats.get('steps', '?'):>5}  "
        f"wall={wall:6.1f}s  T_mean={ba['T_mean_K']:.2f} K  P_mean={ba['P_mean_Pa'] / 1e3:.2f} kPa  "
        f"q_CO2={ba.get('q_CO2_mean', 0):.4f} mol/kg" + ("  (stopped by criterion)" if stopped else ""))
    if float(state_out.T.max()) > cfg.sorbent.max_temperature:
        log(f"    WARNING: bed temperature {state_out.T.max():.1f} K exceeds sorbent.max_temperature")
    return state_out, sr


def _run_scipy(model, y0, dur, t_out, on_step):
    """Cross-check backend: scipy BDF with the same sparsity pattern.
    Accumulators are integrated by the trapezoidal rule on a dense output grid."""
    from scipy.integrate import solve_ivp
    from numerics.sparsity import pattern_and_colors
    g = model.grid
    P, _ = pattern_and_colors(g.Nr, g.Nth, g.Nz, model.ng, model.na, model.n_w > 0)
    t_dense = np.unique(np.concatenate([np.linspace(0, dur, 2001), t_out]))
    sol = solve_ivp(lambda t, y: model.rhs(t, y, want_acc=False)[0], (0.0, dur), y0, method="BDF",
                    t_eval=t_dense, jac_sparsity=P, rtol=model.cfg.solver.rtol, atol=model.atol,
                    first_step=model.cfg.solver.h_initial)
    if not sol.success:
        raise RuntimeError("scipy BDF failed: " + sol.message)
    accs = np.array([model.rhs(t, sol.y[:, i])[1] for i, t in enumerate(sol.t)])
    acc_int = np.concatenate([np.zeros((1, model.n_acc)),
                              np.cumsum(0.5 * (accs[1:] + accs[:-1]) * np.diff(sol.t)[:, None], axis=0)])
    for i, t in enumerate(sol.t):
        if i % 20 == 0:
            on_step(t, sol.y[:, i], acc_int[i])
    snaps = []
    for ts in t_out:
        i = int(np.argmin(np.abs(sol.t - ts)))
        snaps.append((sol.t[i], sol.y[:, i].copy(), acc_int[i].copy()))
    stats = {"method": "scipy_bdf", "steps": len(sol.t), "nfev": sol.nfev, "njev": sol.njev, "nlu": sol.nlu}
    return sol.y[:, -1], acc_int[-1], sol.t[-1], snaps, stats


def _run_rk4(model, y0, dur, t_out, on_step):
    from core.ode_solver import rk4_step
    from numerics.cfl import explicit_dt
    dt = explicit_dt(model, y0)
    n = int(np.ceil(dur / dt))
    if n > 2_000_000:
        raise RuntimeError(f"RK4 needs {n:.3g} steps (dt={dt:.2e} s) for {dur} s -- use solver.method=trbdf2")
    dt = dur / n
    y = y0.copy()
    acc = np.zeros(model.n_acc)
    snaps = [(0.0, y.copy(), acc.copy())]
    targets = list(t_out[1:])
    t = 0.0
    for i in range(n):
        y, da = rk4_step(model, t, y, dt)
        acc = acc + da
        t += dt
        on_step(t, y, acc)
        while targets and t >= targets[0] - 1e-12:
            snaps.append((t, y.copy(), acc.copy()))
            targets.pop(0)
    return y, acc, t, snaps, {"method": "rk4", "steps": n, "dt": dt, "nfev": model.nfev}


# ---------------------------------------------------------------------------- standalone CLI

def step_cli(step_type: str):
    """`python -m steps.<type> --state state.csv [--config c.yaml] [--set k=v ...] [--out dir]`"""
    import sys
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if root not in sys.path:
        sys.path.insert(0, root)
    from config.parameters import load_config
    from core.grid import Grid
    from core.state import initial_state
    from properties import get_backend
    from simio.checkpoint import load_state_csv, save_state_csv
    from simio.export import write_rows, write_fields_csv
    ap = argparse.ArgumentParser(description=f"run the {step_type} step standalone")
    ap.add_argument("--state", help="input state CSV (default: fresh initial state)")
    ap.add_argument("--config", default=None)
    ap.add_argument("--set", action="append", default=[], help="override key=value")
    ap.add_argument("--name", default=None, help="step name in the cycle list (default: first of this type)")
    ap.add_argument("--out", default=None, help="output directory")
    a = ap.parse_args()
    cfg = load_config(a.config, a.set)
    grid = Grid.from_config(cfg)
    props = get_backend(cfg)
    steps = [s for s in cfg.cycle if s["type"] == step_type and (a.name is None or s["name"] == a.name)]
    if not steps:
        raise SystemExit(f"no step of type {step_type} in the cycle definition")
    step = steps[0]
    species = props.species
    adsorbed = ["CO2", "H2O"] if cfg.components.water else ["CO2"]
    st = load_state_csv(a.state, grid, species, adsorbed) if a.state else initial_state(cfg, grid, props)
    out = a.out or os.path.join(root, cfg.run.results_dir, cfg.run.name, "standalone")
    st2, sr = run_step(st, cfg, grid, props, step, logger=print)
    p = save_state_csv(st2, grid, os.path.join(out, f"state_after_{step['name']}.csv"))
    write_rows(os.path.join(out, f"{step['name']}_outlet.csv"), sr.outlet)
    write_fields_csv(os.path.join(out, f"{step['name']}_fields.csv"), grid, sr.snapshots, species, adsorbed)
    print("wrote", p)
