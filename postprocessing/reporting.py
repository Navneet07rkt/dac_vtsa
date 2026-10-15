"""Cycle performance indicators, mass/energy balance tables and log printing.

All quantities are built from the integrated boundary/heat-flow accumulators
of each step (exact time integrals with the solver's own quadrature) and from
the bed inventories at step start/end.

Definitions (per cycle; at CSS these are the process KPIs)
---------------------------------------------------------
fed_i        moles of i entering the bed through any valve (feed + vent back-flow)
product_i    moles of i leaving through the step's product port (V3 by default)
vented_i     all other outflows (adsorption breakthrough, blowdown)
recovery     100 * product_CO2 / fed_CO2                     ("% CO2 recovered")
purity_dry   100 * product_CO2 / (product_CO2 + product_N2)  (water condensed, Young 2021 / Zhang 2025 eq. 9)
productivity kg CO2 / (m3 bed * day) and mol / (kg sorbent * day)
working cap. product_CO2 / sorbent mass                      [mol/kg]
Q_th         heat delivered by the jacket (+ internal exchanger) in steps whose
             jacket is hotter than ambient                   [kJ/kg CO2]
W_el         vacuum pump + feed blower isentropic work / efficiency  [kJ/kg CO2]
W_eq         W_el + eta_turb (1 - T_L/T_H) Q_th  (work-equivalent; Young 2021, Zhang 2025 eq. 14)
"""
from __future__ import annotations

from typing import Dict, List

import numpy as np

M_CO2 = 0.04401
M_H2O = 0.018015
VALVES = ("V1", "V2", "V3")


def _sum(acc, prefix, species):
    return sum(acc.get(f"{prefix}_{v}_{species}", 0.0) for v in VALVES)


def step_flows(sr, species):
    """fed / product / vented moles per species for one step."""
    pp = sr.product_port
    out = {}
    for s in species:
        fed = _sum(sr.acc, "n_in", s)
        allout = _sum(sr.acc, "n_out", s)
        prod = sr.acc.get(f"n_out_{pp}_{s}", 0.0) if (pp in VALVES and sr.valves.get(pp, "CLOSED") != "CLOSED") else 0.0
        out[s] = {"fed": fed, "product": prod, "vented": allout - prod, "out": allout}
    return out


def is_heating_step(sr, cfg) -> bool:
    return sr.jacket_on and sr.T_jacket > cfg.ambient.T + 1.0


def cycle_kpis(steps: List, cfg, species) -> Dict[str, float]:
    tot = {s: {"fed": 0.0, "product": 0.0, "vented": 0.0} for s in species}
    for sr in steps:
        fl = step_flows(sr, species)
        for s in species:
            for k in ("fed", "product", "vented"):
                tot[s][k] += fl[s][k]
    t_cycle = sum(sr.duration for sr in steps)
    m_sorb = cfg.derived.sorbent_mass
    V_bed = cfg.derived.bed_volume
    n_co2 = tot["CO2"]["product"]
    kg = n_co2 * M_CO2
    Q_th = sum(sr.acc["Q_jac_pos"] + sr.acc["Q_hx_pos"] for sr in steps if is_heating_step(sr, cfg))
    Q_cool = -sum(sr.acc["Q_jac_neg"] + sr.acc["Q_hx_neg"] for sr in steps)
    W_vac = sum(sr.acc.get(f"W_pump_{v}", 0.0) for sr in steps for v in VALVES)
    W_fan = sum(sr.acc.get(f"W_fan_{v}", 0.0) for sr in steps for v in VALVES)
    W_el = W_vac + W_fan
    T_H = max([sr.T_jacket for sr in steps if is_heating_step(sr, cfg)] or [cfg.operation.T_heat])
    T_L = cfg.energy.T_low_equiv
    eta_t = cfg.energy.eta_turbine
    W_eq = W_el + eta_t * (1.0 - T_L / T_H) * Q_th
    prod_all = sum(tot[s]["product"] for s in species)
    nn = tot["N2"]["product"]
    # specific energies are undefined when (almost) nothing was produced, e.g. the first cycle
    # from a clean bed: report NaN instead of a meaningless 1e10 kJ/kg
    safe = kg if n_co2 > 1e-6 * max(tot["CO2"]["fed"], 1e-30) else np.nan
    k = {
        "cycle_time_s": t_cycle,
        "CO2_fed_mol": tot["CO2"]["fed"],
        "CO2_product_mol": n_co2,
        "CO2_vented_mol": tot["CO2"]["vented"],
        "CO2_product_kg": kg,
        "recovery_pct": 100.0 * n_co2 / tot["CO2"]["fed"] if tot["CO2"]["fed"] > 0 else np.nan,
        "purity_dry_pct": 100.0 * n_co2 / (n_co2 + nn) if (n_co2 + nn) > 0 else np.nan,
        "purity_wet_pct": 100.0 * n_co2 / prod_all if prod_all > 0 else np.nan,
        "N2_product_mol": nn,
        "productivity_kg_m3_day": kg / V_bed / (t_cycle / 86400.0),
        "productivity_mol_kg_day": n_co2 / m_sorb / (t_cycle / 86400.0),
        "capture_rate_kg_day": kg / (t_cycle / 86400.0),
        "working_capacity_mol_kg": n_co2 / m_sorb,
        "Q_thermal_J": Q_th,
        "Q_cooling_removed_J": Q_cool,
        "W_vacuum_J": W_vac,
        "W_blower_J": W_fan,
        "specific_thermal_kJ_per_kg": Q_th / safe / 1e3,
        "specific_electric_kJ_per_kg": W_el / safe / 1e3,
        "specific_vacuum_kJ_per_kg": W_vac / safe / 1e3,
        "specific_blower_kJ_per_kg": W_fan / safe / 1e3,
        "specific_total_kJ_per_kg": (Q_th + W_el) / safe / 1e3,
        "work_equivalent_kJ_per_kg": W_eq / safe / 1e3,
        "T_heat_K": T_H,
        "sorbent_mass_kg": m_sorb,
        "bed_volume_m3": V_bed,
    }
    if "H2O" in species:
        k["H2O_product_mol"] = tot["H2O"]["product"]
        k["H2O_per_CO2_kg_kg"] = tot["H2O"]["product"] * M_H2O / safe
        k["H2O_fed_mol"] = tot["H2O"]["fed"]
    return k


def mass_balance_rows(sr, species, floors=None) -> List[Dict]:
    rows = []
    floors = floors or {}
    for s in species:
        a, b = sr.inv_start[f"tot_{s}"], sr.inv_end[f"tot_{s}"]
        fin = _sum(sr.acc, "n_in", s)
        fout = _sum(sr.acc, "n_out", s)
        res = (b - a) - (fin - fout)
        scale = max(abs(fin), abs(fout), abs(a), abs(b), floors.get(s, 1e-30))
        rows.append({"cycle": sr.cycle, "step": sr.name, "species": s, "inventory_start_mol": a,
                     "inventory_end_mol": b, "in_mol": fin, "out_mol": fout, "residual_mol": res,
                     "relative_error": res / scale})
    return rows


def energy_balance_row(sr, floor: float = 0.0) -> Dict:
    acc = sr.acc
    E0, E1 = sr.inv_start["E_total"], sr.inv_end["E_total"]
    H_in = sum(acc.get(f"H_in_{v}", 0.0) for v in VALVES)
    H_out = sum(acc.get(f"H_out_{v}", 0.0) for v in VALVES)
    Q_j = acc["Q_jac_pos"] + acc["Q_jac_neg"]
    Q_hx = acc["Q_hx_pos"] + acc["Q_hx_neg"]
    Q_amb = acc["Q_amb"]
    res = (E1 - E0) - (H_in - H_out + Q_j + Q_hx + Q_amb)
    scale = max(abs(E1 - E0), abs(H_in), abs(H_out), abs(Q_j), abs(Q_hx), abs(Q_amb), floor, 1e-30)
    return {"cycle": sr.cycle, "step": sr.name, "E_start_J": E0, "E_end_J": E1, "dE_J": E1 - E0,
            "H_in_J": H_in, "H_out_J": H_out, "Q_jacket_J": Q_j, "Q_internal_hx_J": Q_hx,
            "Q_ambient_J": Q_amb, "residual_J": res, "relative_error": res / scale}


def energy_breakdown_row(sr, cfg) -> Dict:
    """Where the heat of one step went (sums exactly to the jacket + exchanger input)."""
    acc = sr.acc
    dH = {"CO2": cfg.thermo.dH_ads_co2, "H2O": cfg.thermo.dH_ads_h2o}
    q_des = {s: dH[s] * (sr.ads_end[s] - sr.ads_start[s]) for s in sr.ads_start}   # >0 when desorbing
    dE_bed = sr.inv_end["E_bed"] - sr.inv_start["E_bed"]
    dE_wall = sr.inv_end["E_wall"] - sr.inv_start["E_wall"]
    H_net_out = (sum(acc.get(f"H_out_{v}", 0.0) for v in VALVES)
                 - sum(acc.get(f"H_in_{v}", 0.0) for v in VALVES))
    Q_in = acc["Q_jac_pos"] + acc["Q_jac_neg"] + acc["Q_hx_pos"] + acc["Q_hx_neg"]
    row = {"cycle": sr.cycle, "step": sr.name, "heat_input_J": Q_in,
           "desorption_CO2_J": q_des.get("CO2", 0.0), "desorption_H2O_J": q_des.get("H2O", 0.0),
           "bed_sensible_J": dE_bed - sum(q_des.values()), "wall_sensible_J": dE_wall,
           "gas_enthalpy_out_J": H_net_out, "ambient_loss_J": -acc["Q_amb"]}
    row["closure_J"] = Q_in - (row["desorption_CO2_J"] + row["desorption_H2O_J"] + row["bed_sensible_J"]
                               + row["wall_sensible_J"] + row["gas_enthalpy_out_J"] + row["ambient_loss_J"])
    return row


def format_kpis(k: Dict[str, float]) -> str:
    lines = [
        f"  CO2 recovered            : {k['recovery_pct']:8.2f} %",
        f"  CO2 purity (dry / wet)   : {k['purity_dry_pct']:8.2f} % / {k['purity_wet_pct']:.2f} %",
        f"  working capacity         : {k['working_capacity_mol_kg']:8.4f} mol/kg",
        f"  productivity             : {k['productivity_kg_m3_day']:8.2f} kg/m3/day   "
        f"({k['productivity_mol_kg_day']:.3f} mol/kg/day)",
        f"  specific thermal energy  : {k['specific_thermal_kJ_per_kg']:8.0f} kJ/kg CO2",
        f"  specific electric energy : {k['specific_electric_kJ_per_kg']:8.0f} kJ/kg CO2 "
        f"(vacuum {k['specific_vacuum_kJ_per_kg']:.0f}, blower {k['specific_blower_kJ_per_kg']:.0f})",
        f"  total (thermal+electric) : {k['specific_total_kJ_per_kg']:8.0f} kJ/kg CO2",
        f"  work equivalent          : {k['work_equivalent_kJ_per_kg']:8.0f} kJ/kg CO2",
        f"  cycle time               : {k['cycle_time_s']:8.0f} s",
    ]
    if "H2O_per_CO2_kg_kg" in k:
        lines.append(f"  water co-desorbed        : {k['H2O_per_CO2_kg_kg']:8.2f} kg H2O / kg CO2")
    return "\n".join(lines)
