#!/usr/bin/env python3
"""CLI entry point (a.k.a. main.py): run the 3D TVSA DAC cycle to cyclic steady state.

    python run_dac_simulation.py                              # defaults.yaml
    python run_dac_simulation.py --config config/cases/quick_test.yaml
    python run_dac_simulation.py --set grid.Nr=10 --set operation.P_vacuum=15000 --name pvac15

Outputs (results/<run name>/):
    C_out_vs_t.csv                    outlet composition / C_out/C_feed / flows / P / T vs time (final cycle)
    fields_q_T_P_r_theta_z_t.csv      q_CO2, q_H2O, T, P, c, y at every (r, theta, z) cell and snapshot time
    fields_final_cycle.npz            same, binary (used by postprocessing/probe.py)
    kpis.csv / kpis.json              % CO2 recovered, purity, productivity, kJ/kg CO2 (thermal, electric,
                                      total, work-equivalent), working capacity ... per cycle / final
    mass_balance.csv                  per step & species: inventory change vs in - out, residual (MB)
    energy_balance.csv                per step: dE vs enthalpy flows + heat, residual (EB)
    energy_breakdown.csv              where the heating energy went (desorption, sensible, wall, losses)
    css_history.csv                   convergence to cyclic steady state
    cycle_NNN/state_after_<step>.csv  the hand-off state of every step (restartable)
    plots/*.png, summary.txt, vtsa_main_run.log, config_used.yaml
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import yaml  # noqa: E402

from config.parameters import load_config  # noqa: E402
from postprocessing import reporting  # noqa: E402
from simio.export import write_fields_csv, write_npz, write_rows  # noqa: E402
from vtsa_cycle import Logger, run_to_css  # noqa: E402


def write_outputs(res, run_dir, log, make_plots=True):
    cfg, grid = res.cfg, res.grid
    species = list(res.final_state.species)
    adsorbed = list(res.final_state.adsorbed)
    steps = res.final_steps
    last = res.cycles[-1]
    paths = {}
    paths["C_out_vs_t"] = write_rows(os.path.join(run_dir, "C_out_vs_t.csv"), [r for sr in steps for r in sr.outlet])
    snaps = [s for sr in steps for s in sr.snapshots]
    paths["fields"] = write_fields_csv(os.path.join(run_dir, "fields_q_T_P_r_theta_z_t.csv"), grid, snaps,
                                       species, adsorbed)
    write_npz(os.path.join(run_dir, "fields_final_cycle.npz"), grid, snaps)
    rows = [{"cycle": c.cycle, **c.kpis, "css_change": c.change["max"], "co2_cycle_balance": c.co2_balance,
             "balances_ok": c.balances_ok, "wall_time_s": c.wall_time} for c in res.cycles]
    paths["kpis"] = write_rows(os.path.join(run_dir, "kpis.csv"), rows)
    final = dict(last.kpis)
    final.update({"cycles_run": len(res.cycles), "css_reached": res.converged,
                  "max_mass_balance_rel_error": max(abs(r["relative_error"]) for r in last.mb_rows),
                  "max_energy_balance_rel_error": max(abs(r["relative_error"]) for r in last.eb_rows)})
    with open(os.path.join(run_dir, "kpis.json"), "w") as fh:
        json.dump({k: (float(v) if isinstance(v, (int, float)) else v) for k, v in final.items()}, fh, indent=2)
    paths["mass_balance"] = write_rows(os.path.join(run_dir, "mass_balance.csv"),
                                       [r for c in res.cycles for r in c.mb_rows])
    eb_cols = ["cycle", "step", "E_start_J", "E_end_J", "dE_J", "H_in_J", "H_out_J", "Q_jacket_J",
               "Q_internal_hx_J", "Q_ambient_J", "residual_J", "relative_error"]
    paths["energy_balance"] = write_rows(os.path.join(run_dir, "energy_balance.csv"),
                                         [r for c in res.cycles for r in c.eb_rows], eb_cols)
    kg = last.kpis["CO2_product_kg"]
    bd_rows = []
    for r in last.breakdown:
        rr = dict(r)
        for k in list(r.keys()):
            if k.endswith("_J"):
                rr[k.replace("_J", "_kJ_per_kgCO2")] = r[k] / kg / 1e3 if kg > 0 else float("nan")
        bd_rows.append(rr)
    paths["energy_breakdown"] = write_rows(os.path.join(run_dir, "energy_breakdown.csv"), bd_rows)
    paths["css"] = write_rows(os.path.join(run_dir, "css_history.csv"),
                              [{"cycle": c.cycle, **c.change, "co2_cycle_balance": c.co2_balance} for c in res.cycles])
    solver_rows = [{"step": sr.name, "duration_s": sr.duration, "wall_s": sr.wall_time, **sr.solver} for sr in steps]
    write_rows(os.path.join(run_dir, "solver_stats_final_cycle.csv"), solver_rows,
               sorted({k for r in solver_rows for k in r}, key=lambda x: (x != "step", x)))
    if make_plots:
        from postprocessing import plotting as P
        pd = os.path.join(run_dir, "plots")
        P.plot_outlet([r for sr in steps for r in sr.outlet], os.path.join(pd, "outlet_and_bed_history.png"),
                      f"{cfg.run.name}: final cycle")
        for sr in steps:
            sel = sr.snapshots[:: max(1, len(sr.snapshots) // 4)][:4] + [sr.snapshots[-1]]
            for fld in ["q_CO2", "T"] + (["q_H2O"] if "H2O" in adsorbed else []):
                P.plot_zr_slice(grid, sel, fld, os.path.join(pd, f"zr_{sr.name}_{fld}.png"))
            if grid.Nth > 1:
                P.plot_polar(grid, sr.snapshots[-1], "T", os.path.join(pd, f"polar_T_end_{sr.name}.png"))
                P.plot_polar(grid, sr.snapshots[-1], "q_CO2", os.path.join(pd, f"polar_qCO2_end_{sr.name}.png"))
        P.plot_probe_histories(grid, snaps, os.path.join(pd, "T_probes_cycle.png"), "T")
        P.plot_probe_histories(grid, snaps, os.path.join(pd, "qCO2_probes_cycle.png"), "q_CO2")
        P.plot_css(res.cycles, os.path.join(pd, "css_convergence.png"))
        heat = [r for r in last.breakdown if r["heat_input_J"] > 0]
        if heat:
            hb = max(heat, key=lambda r: r["heat_input_J"])
            P.plot_energy_breakdown(hb, os.path.join(pd, "energy_breakdown.png"), kg)
    summary = ["3D TVSA DAC simulation -- " + cfg.run.name,
               grid.describe(),
               f"cycles run: {len(res.cycles)}  (CSS {'reached' if res.converged else 'NOT reached'})",
               reporting.format_kpis(last.kpis),
               f"  mass balance  (max rel. residual, last cycle): {final['max_mass_balance_rel_error']:.2e}",
               f"  energy balance(max rel. residual, last cycle): {final['max_energy_balance_rel_error']:.2e}"]
    with open(os.path.join(run_dir, "summary.txt"), "w") as fh:
        fh.write("\n".join(summary) + "\n")
    log("\n".join(summary))
    return paths


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None, help="YAML file overriding defaults.yaml")
    ap.add_argument("--set", action="append", default=[], help="override, e.g. --set grid.Nr=10")
    ap.add_argument("--name", default=None, help="run name (results/<name>)")
    ap.add_argument("--no-plots", action="store_true")
    a = ap.parse_args(argv)
    sets = list(a.set)
    if a.name:
        sets.append(f"run.name={a.name}")
    cfg = load_config(a.config, sets)
    run_dir = os.path.join(ROOT, cfg.run.results_dir, cfg.run.name)
    os.makedirs(run_dir, exist_ok=True)
    log = Logger(os.path.join(run_dir, "vtsa_main_run.log"))
    with open(os.path.join(run_dir, "config_used.yaml"), "w") as fh:
        raw = cfg.to_dict()
        raw.pop("derived", None)        # recomputed on load; keeps the file re-loadable as --config
        yaml.safe_dump(raw, fh, sort_keys=False)
    log(f"dac_vtsa run '{cfg.run.name}'  ->  {run_dir}")
    t0 = time.time()
    res = run_to_css(cfg, run_dir, log)
    write_outputs(res, run_dir, log, make_plots=cfg.run.make_plots and not a.no_plots)
    log(f"total wall time {time.time() - t0:.0f} s")
    return res


if __name__ == "__main__":
    main()
