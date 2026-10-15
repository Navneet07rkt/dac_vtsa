#!/usr/bin/env python3
"""Regenerate PARAMETERS.md: the hand-written guide below + the full table from config/schema.py.

    python tools/make_parameters_md.py
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from config import schema  # noqa: E402

HEADER = r"""# PARAMETERS.md — every parameter you can change, its default, and where it comes from

All parameters live in `config/defaults.yaml` (SI units: m, s, K, Pa, mol, kg, J, W).
Every key is declared in `config/schema.py` with its type, unit, allowed range/choices,
meaning and **source**; the table at the end of this file is generated from it
(`python tools/make_parameters_md.py`). Unknown keys and out-of-range values are
rejected when a run starts, so a typo fails loudly instead of silently using a default.

## How to change a parameter

| way | example | use it for |
|---|---|---|
| command line | `python run_dac_simulation.py --set operation.T_heat=383.15 --set grid.Nr=10 --name hot` | one-off experiments |
| case file | `python run_dac_simulation.py --config config/cases/asymmetric_jacket.yaml` (the file lists only the keys you change) | scenarios you want to keep |
| `config/defaults.yaml` | edit the value | changing the baseline for every run |
| sensitivity | list the key under `sensitivity.parameters` and run `python run_sensitivity.py` | how much a KPI depends on it |

`--set` and case files can be combined; `--set` wins. The `config_used.yaml` written
into every results folder is a complete, re-loadable record of what was run
(`--config results/<run>/config_used.yaml` reproduces it).

## The cycle

`cycle:` is an ordered list of steps. A case file that contains `cycle:` replaces the
whole default list. Each step:

```yaml
- name: heating                 # unique; used in file names (state_after_heating.csv)
  type: heating                 # adsorption | blowdown | heating | cooling | pressurization
  duration: 22400.0             # s (maximum duration if stop_on is given)
  valves: {V1: CLOSED, V2: CLOSED, V3: OPEN_VACUUM}   # CLOSED | OPEN (feed, V1) | OPEN_ATMOSPHERIC | OPEN_VACUUM
  jacket: {mode: temperature, T: operation.T_heat}   # number in K, or a reference to another key; mode: adiabatic
  stop_on: {criterion: end_saturation, value: 0.9895}  # optional: outlet_ratio | end_saturation | min_temperature | max_temperature
  product_port: V3              # optional (default V3): what leaves here counts as product
```

`ports:` says where each valve sits (V1 top, V2 and V3 bottom by default).

## The parameters that matter most — start here

**Process design** (what an engineer actually chooses)

| key | default | why it matters |
|---|---|---|
| `cycle[i].duration` | 8000 / 60 / 22400 / 400 s | adsorption time sets how much CO2 slips (recovery); heating time sets how much is desorbed (working capacity) and the heat duty |
| `operation.T_heat` | 373.15 K | regeneration temperature; Lewatit oxidises above ~100 C (Young 2021) |
| `operation.P_vacuum` | 20 kPa | lower P = more desorption and purity, more pump work |
| `operation.pump_speed_per_area` | 5e-3 m/s | pump size per bed cross-section; sets blowdown speed and the pressure during desorption |
| `feed.velocity` | 0.0706 m/s | superficial feed velocity (Young: 90 % of minimum fluidisation); productivity vs. pressure drop |
| `jacket.h_jacket` | 26 W/m2/K | jacket-to-wall coefficient (Galanti 2025); with a rim jacket it is often THE bottleneck of heating |
| `geometry.bed_radius`, `geometry.bed_length` | 0.025 m, 0.01 m | conduction length for rim heating; bed inventory |
| `internal_hx.*` | off | embedded exchanger — the usual engineering fix for conduction-limited beds |
| `jacket.theta_coverage_deg`, `jacket.z_coverage`, `feed.profile*` | full, uniform | make the problem genuinely 3D (partial jackets, maldistributed feed) |
| `feed.T`, `feed.RH`, `feed.y_co2_dry` | 288 K, 0.55, 400 ppm | climate of the site; humidity changes both capacity and water co-desorption |

**Sorbent (Lewatit VP OC 1065)** — fitted by Young et al. 2021 to their own measurements;
change these only together with a new data set: `isotherm.co2.*` (Toth), `isotherm.h2o.*`
(GAB), `isotherm.coadsorption.*`, `kinetics.k_ldf_*` (the LDF constants are lumped fits —
the most uncertain numbers in the model, so run them through the sensitivity analysis),
`thermo.dH_ads_*`, `sorbent.*`.

**Transport and heat-transfer correlations** — literature correlations with tunable
constants: `transport.*` (dispersion, effective conductivity, wall coefficient, Ergun).

**Energy accounting** — `energy.*` (pump, blower and heat-to-work efficiencies; they
change the kJ/kg numbers, not the process).

**Numerics** — `grid.*` (resolution), `solver.*` (tolerances, Newton/Krylov settings).
Defaults are tested; refine the grid to check that a result is grid-independent before
you trust it (halve the cells in each direction and compare the KPIs).

## Sources

* **[Y21]** Young, García-Díez, Garcia, van der Spek (2021) *Energy Environ. Sci.* 14, 5377–5394,
  doi:10.1039/D1EE01272J, and its ESI (Tables S5, S9, S10) — isotherms, kinetics, heats,
  feed, cycle, efficiencies; validation case.
* **[G25]** Galanti et al. (2025) *RSC Energy Advances*, doi:10.1039/D5YA00336A, SI Tables S7–S8 —
  wall and jacket properties, heat-transfer correlations.
* **[Z25]** Zhang et al. (2025) *Energy* 334, 137733 — sorbent thermal conductivity.
  (Its Table A.1 contains two typos for the Lewatit isotherm — see MODEL.md §13.)
* **[ER68]** Edwards & Richardson (1968) *Chem. Eng. Sci.* 23, 109 — axial dispersion.
* **[ZS70]** Zehner & Schlünder (1970) *Chem.-Ing.-Tech.* 42, 933 — stagnant bed conductivity.
* **[MN93]** Martin & Nilles (1993) *Chem.-Ing.-Tech.* 65, 1468 — bed-to-wall heat transfer.
* "modelling assumption / user choice" — no literature value exists or it is a design
  decision (grid, wall thickness, solver settings, partial jackets ...).

## Full list (generated from `config/schema.py`)

"""

if __name__ == "__main__":
    out = os.path.join(ROOT, "PARAMETERS.md")
    with open(out, "w") as fh:
        fh.write(HEADER)
        fh.write(schema.markdown_table())
        fh.write("\n")
    print("wrote", out)
