"""ORCHESTRATOR ONLY -- chains the cycle steps, no physics here.

For every cycle:  for every step in cfg.cycle (default adsorption -> blowdown ->
heating -> cooling):

    state  <- load_state_csv(previous step's CSV)        (run.handoff_via_csv)
    state, result <- steps.<type>.run_<type>(state, ...)
    save_state_csv(state, cycle_NNN/state_after_<step>.csv)

then balances (validation/conservation_checks), KPIs (postprocessing/reporting)
and the CSS test (validation/cyclic_steady_state) decide whether to continue.
"""
from __future__ import annotations

import json
import os
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from core.grid import Grid
from core.state import initial_state
from postprocessing import reporting
from properties import get_backend
from simio.checkpoint import load_state_csv, save_state_csv
from simio.export import write_fields_csv, write_rows
from steps import RUNNERS
from validation.conservation_checks import check_cycle
from validation.cyclic_steady_state import cycle_co2_balance, is_css, state_change


class Logger:
    def __init__(self, path: Optional[str] = None, echo: bool = True):
        self.path = path
        self.echo = echo
        if path:
            os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
            open(path, "w").close()

    def __call__(self, msg: str = ""):
        if self.echo:
            print(msg, flush=True)
        if self.path:
            with open(self.path, "a") as fh:
                fh.write(msg + "\n")


@dataclass
class CycleRecord:
    cycle: int
    kpis: Dict[str, float]
    mb_rows: List[Dict]
    eb_rows: List[Dict]
    breakdown: List[Dict]
    change: Dict[str, float]
    co2_balance: float
    balances_ok: bool
    wall_time: float


@dataclass
class RunResult:
    cfg: object
    grid: Grid
    run_dir: str
    cycles: List[CycleRecord] = field(default_factory=list)
    final_steps: List = field(default_factory=list)
    final_state: object = None
    converged: bool = False


def run_cycle(state, cfg, grid, props, cycle_idx: int, cycle_dir: str, log, T_jac_prev=None):
    steps_out = []
    t_cycle = 0.0
    prev_csv = os.path.join(cycle_dir, "state_cycle_start.csv")
    save_state_csv(state, grid, prev_csv)
    species = list(props.species)
    adsorbed = list(state.adsorbed)
    for step in cfg.cycle:
        if cfg.run.handoff_via_csv:
            state = load_state_csv(prev_csv, grid, species, adsorbed)
        runner = RUNNERS[step["type"]]
        state, sr = runner(state, cfg, grid, props, step=step, cycle_idx=cycle_idx, t_cycle0=t_cycle,
                           T_jacket_prev=T_jac_prev, logger=log)
        prev_csv = os.path.join(cycle_dir, f"state_after_{step['name']}.csv")
        save_state_csv(state, grid, prev_csv)
        t_cycle += sr.duration
        T_jac_prev = sr.T_jacket if sr.jacket_on else T_jac_prev
        steps_out.append(sr)
    return state, steps_out, T_jac_prev


def run_to_css(cfg, run_dir: str, log=None) -> RunResult:
    log = log or Logger(None)
    grid = Grid.from_config(cfg)
    props = get_backend(cfg)
    species = list(props.species)
    state = initial_state(cfg, grid, props)
    res = RunResult(cfg=cfg, grid=grid, run_dir=run_dir)
    log(grid.describe())
    log(f"species {species}; isotherm co-adsorption model: {cfg.isotherm.coadsorption.model if cfg.components.water else 'dry Toth'}")
    log("cycle: " + " -> ".join(f"{s['name']}({s['duration']:.0f}s)" for s in cfg.cycle))
    T_jac_prev = None
    prev_end = state.copy()
    for n in range(1, cfg.run.max_cycles + 1):
        t0 = time.time()
        cdir = os.path.join(run_dir, f"cycle_{n:03d}")
        log(f"\n=== cycle {n} ===")
        state, steps, T_jac_prev = run_cycle(state, cfg, grid, props, n, cdir, log, T_jac_prev)
        mb, eb, ok = check_cycle(steps, species)
        kp = reporting.cycle_kpis(steps, cfg, species)
        bd = [reporting.energy_breakdown_row(sr, cfg) for sr in steps]
        ch = state_change(prev_end, state, grid, cfg)
        co2b = cycle_co2_balance(kp)
        rec = CycleRecord(n, kp, mb, eb, bd, ch, co2b, ok, time.time() - t0)
        res.cycles.append(rec)
        write_rows(os.path.join(cdir, "C_out_vs_t.csv"), [r for sr in steps for r in sr.outlet])
        if cfg.run.write_fields_every_cycle:
            snaps = [s for sr in steps for s in sr.snapshots]
            write_fields_csv(os.path.join(cdir, "fields_q_T_P_r_theta_z_t.csv"), grid, snaps, species,
                             list(state.adsorbed))
        worst_m = max(abs(r["relative_error"]) for r in mb)
        worst_e = max(abs(r["relative_error"]) for r in eb)
        log(f"  recovery {kp['recovery_pct']:.2f} %  purity {kp['purity_dry_pct']:.2f} %  "
            f"WC {kp['working_capacity_mol_kg']:.3f} mol/kg  thermal {kp['specific_thermal_kJ_per_kg']:.0f} kJ/kg  "
            f"| MB err {worst_m:.1e}  EB err {worst_e:.1e} {'OK' if ok else 'CHECK!'}  "
            f"| CSS change {ch['max']:.2e}  CO2 cycle balance {co2b:.2e}  | {rec.wall_time:.0f} s")
        res.final_steps = steps
        res.final_state = state
        if is_css(ch, kp, cfg, n):
            res.converged = True
            log(f"  cyclic steady state reached after {n} cycles")
            break
        prev_end = state.copy()
    if not res.converged:
        log(f"  WARNING: CSS not reached in {cfg.run.max_cycles} cycles (last change "
            f"{res.cycles[-1].change['max']:.2e}); results are from the last cycle")
    return res
