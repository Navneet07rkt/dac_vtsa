"""Mass / energy balance closure check for every step.

For each step and species:  (inventory_end - inventory_start) - (in - out) = residual
and for energy:             dE_bed+wall - (H_in - H_out + Q_jacket + Q_hx + Q_ambient) = residual

Because the finite-volume fluxes telescope exactly and the accumulators are
integrated with the solver's own stage formulas, residuals are limited only
by the Newton tolerance (typically 1e-9 .. 1e-6 relative).  Anything larger
indicates a bug or a step that did not converge -- the run log flags it.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

from postprocessing.reporting import energy_balance_row, mass_balance_rows

MASS_TOL = 1e-4
ENERGY_TOL = 1e-4
FLOOR = 1e-3          # see check_cycle


def check_step(sr, species, floors=None, e_floor: float = 0.0) -> Tuple[List[Dict], Dict, bool]:
    mb = mass_balance_rows(sr, species, floors)
    eb = energy_balance_row(sr, e_floor)
    ok = all(abs(r["relative_error"]) < MASS_TOL for r in mb) and abs(eb["relative_error"]) < ENERGY_TOL
    return mb, eb, ok


def check_cycle(step_results, species):
    mb_rows, eb_rows, all_ok = [], [], True
    # Relative errors of a step are taken w.r.t. max(the step's own flows, FLOOR x the cycle
    # throughput of that species / of energy).  Without the floor a step that moves almost nothing
    # (N2 in heating, energy in a 60 s blowdown) reports a large *relative* error for a residual
    # that is ~1e-8 of what the cycle moves.  The absolute residuals are always in the CSVs.
    floors = {}
    for s in species:
        thr = sum(sum(sr.acc.get(f"n_in_{v}_{s}", 0.0) + sr.acc.get(f"n_out_{v}_{s}", 0.0)
                      for v in ("V1", "V2", "V3")) for sr in step_results)
        floors[s] = FLOOR * max(thr, 1e-30)
    e_thr = 0.0
    for sr in step_results:
        a = sr.acc
        e_thr += (sum(abs(a.get(f"H_in_{v}", 0.0)) + abs(a.get(f"H_out_{v}", 0.0)) for v in ("V1", "V2", "V3"))
                  + abs(a["Q_jac_pos"]) + abs(a["Q_jac_neg"]) + abs(a["Q_hx_pos"]) + abs(a["Q_hx_neg"])
                  + abs(a["Q_amb"]))
    e_floor = FLOOR * e_thr
    for sr in step_results:
        mb, eb, ok = check_step(sr, species, floors, e_floor)
        mb_rows += mb
        eb_rows.append(eb)
        all_ok &= ok
    # whole-cycle rows
    for s in species:
        a = step_results[0].inv_start[f"tot_{s}"]
        b = step_results[-1].inv_end[f"tot_{s}"]
        fin = sum(r["in_mol"] for r in mb_rows if r["species"] == s)
        fout = sum(r["out_mol"] for r in mb_rows if r["species"] == s)
        res = (b - a) - (fin - fout)
        mb_rows.append({"cycle": step_results[0].cycle, "step": "CYCLE", "species": s,
                        "inventory_start_mol": a, "inventory_end_mol": b, "in_mol": fin, "out_mol": fout,
                        "residual_mol": res, "relative_error": res / max(abs(fin), abs(fout), abs(a), 1e-30)})
    E0 = step_results[0].inv_start["E_total"]
    E1 = step_results[-1].inv_end["E_total"]
    keys = ("H_in_J", "H_out_J", "Q_jacket_J", "Q_internal_hx_J", "Q_ambient_J")
    tot = {k: sum(r[k] for r in eb_rows) for k in keys}
    res = (E1 - E0) - (tot["H_in_J"] - tot["H_out_J"] + tot["Q_jacket_J"] + tot["Q_internal_hx_J"]
                       + tot["Q_ambient_J"])
    scale = max(abs(E1 - E0), *[abs(v) for v in tot.values()], 1e-30)
    eb_rows.append({"cycle": step_results[0].cycle, "step": "CYCLE", "E_start_J": E0, "E_end_J": E1,
                    "dE_J": E1 - E0, **tot, "residual_J": res, "relative_error": res / scale})
    return mb_rows, eb_rows, all_ok
