# dac_vtsa — 3D temperature–vacuum swing adsorption (TVSA) simulator for direct air capture

A finite-volume simulator of a jacketed packed bed of amine sorbent (Lewatit VP OC 1065)
capturing CO2 from humid air, resolved in **(r, θ, z)**. It runs the cycle
**adsorption → blowdown → heating → cooling** to cyclic steady state, each step handing its
full state to the next through a CSV file, and reports CO2 recovery, purity, productivity,
specific energy (kJ/kg CO2), and **mass and energy balances** for every step.

Every parameter is sourced (mostly Young et al. 2021, *Energy Environ. Sci.* 14, 5377) and
documented in [PARAMETERS.md](PARAMETERS.md); every equation, boundary condition and
numerical choice is in [MODEL.md](MODEL.md).

---------------------------------------------------------------------------------------

## Quickstart

```bash
pip install -r requirements.txt            # numpy, scipy, pyyaml, matplotlib (+ pytest)

python -m pytest -q tests/                  # 19 tests, ~35 s (one integrates 2 ms with explicit RK4)

python run_dac_simulation.py --config config/cases/quick_test.yaml     # ~1 min smoke test
python run_dac_simulation.py                                           # default 3D case
python run_dac_simulation.py --config config/cases/young2021_1d_equivalent.yaml
python -m validation.young2021_compare                                 # vs Young 2021 Table S10
```

Change anything from the command line or with a case file (see PARAMETERS.md):

```bash
python run_dac_simulation.py --set operation.T_heat=383.15 --set operation.P_vacuum=15000 --name hot_deep_vacuum
python run_dac_simulation.py --config config/cases/asymmetric_jacket.yaml
```

## What a run produces — `results/<run name>/`

| file | content |
|---|---|
| `C_out_vs_t.csv` | final cycle, every time step: valve state, outlet flow, **C_CO2,out / C_CO2,feed**, outlet mole fractions, P and T at both ends, bed-mean T, P, loadings, cumulative CO2 out per valve, jacket heat |
| `fields_q_T_P_r_theta_z_t.csv` | **q_CO2, q_H2O, T, P, c_k, y_k at every (r, θ, z) cell** and snapshot time of the final cycle (long format, one row per cell and time) |
| `fields_final_cycle.npz` | the same fields in binary (used by the probe, below) |
| `kpis.csv`, `kpis.json` | per cycle / final: **% CO2 recovered**, purity (dry/wet), working capacity, productivity, **kJ/kg CO2** (thermal, electric = vacuum + blower, total, work-equivalent), water co-desorbed, balance errors |
| `mass_balance.csv` | **MB** per step, per species and per cycle: inventory start/end, in, out, residual, relative error |
| `energy_balance.csv` | **EB** per step and per cycle: ΔE (bed + wall), enthalpy in/out, jacket/exchanger/ambient heat, residual |
| `energy_breakdown.csv` | where the heat went: CO2 and H2O desorption, bed and wall sensible heat, gas enthalpy out, losses (J and kJ/kg CO2) |
| `css_history.csv` | approach to cyclic steady state |
| `cycle_NNN/state_after_<step>.csv` | the hand-off state of every step (restartable, exact 17-digit round trip) |
| `plots/*.png` | outlet history, (z, r) maps of q and T, polar maps (3D), probe histories, CSS convergence, energy breakdown |
| `summary.txt`, `vtsa_main_run.log`, `config_used.yaml`, `solver_stats_final_cycle.csv` | run record (`config_used.yaml` reloads with `--config`) |

### Any quantity at any (r, θ, z, t)

```bash
python -m postprocessing.probe results/default_3d --field q_CO2 --r 0.012 --theta 90 --z 0.005 --t-cycle 3600
python -m postprocessing.probe results/default_3d --field T --r 0 --z 0.008 --step heating --t-step 5000
python -m postprocessing.probe results/default_3d --field RH --r 0.02 --z 0.002 --series      # whole cycle
```

Fields: `T`, `P`, `Tw` (wall), `c_<gas>`, `y_<gas>`, `p_<gas>` (partial pressure), `q_CO2`,
`q_H2O`, `q_star_CO2`, `q_star_H2O` (local equilibrium loading), `RH`. Values are
interpolated linearly in r, z, periodically in θ, and in time between snapshots
(`run.snapshots_per_step` sets their density). From Python:

```python
from postprocessing.probe import FieldProbe
pr = FieldProbe("results/default_3d")
pr.value("q_CO2", r=0.012, theta_deg=90, z=0.005, t_cycle=3600)
t, values, steps = pr.series("T", r=0.0, z=0.005)
```

## The cycle and the CSV hand-off

`vtsa_cycle.py` runs the steps listed under `cycle:` in `config/defaults.yaml`
(default, Young 2021 timing):

| step | V1 (top) | V2 (bottom) | V3 (bottom, product) | jacket | duration |
|---|---|---|---|---|---|
| adsorption | OPEN (feed blower) | OPEN_ATMOSPHERIC | CLOSED | 288 K | 8000 s |
| blowdown | CLOSED | OPEN_VACUUM | CLOSED | 288 K | 60 s |
| heating | CLOSED | CLOSED | OPEN_VACUUM | 373 K | 22 400 s |
| cooling | CLOSED | CLOSED | CLOSED | 288 K | 400 s |

After each step the complete state (c_k, q_j, T in every cell, wall T) is written to
`cycle_NNN/state_after_<step>.csv`; with `run.handoff_via_csv: true` (default) the next
step **reads that file** — so any step can be rerun or inspected on its own:

```bash
python -m steps.heating --state results/default_3d/cycle_003/state_after_blowdown.csv   # -> results/default_3d/standalone/
```

Steps can be added, removed, re-timed, or ended by a criterion (`stop_on`: outlet
CO2 ratio, end-of-column saturation, temperature); see PARAMETERS.md → "The cycle".

## Folder map

```
run_dac_simulation.py        main entry point (runs to CSS, writes all outputs)
run_sensitivity.py           sensitivity analysis entry point
vtsa_cycle.py                orchestrator: steps, CSV hand-off, balances, KPIs, CSS test
config/                      defaults.yaml (all parameters), schema.py (units, ranges, sources), cases/
core/                        grid, state, residual (the semi-discrete model), Jacobian, implicit solver
physics/                     isotherms, kinetics, dispersion, heat transfer, Ergun, mass and energy fluxes
boundaries/, valves/         valve states -> boundary fluxes (feed, vents, vacuum pump, closed ends)
properties/                  gas properties (ideal-gas correlations; optional CoolProp)
numerics/                    stencils/limiters, sparsity + colouring, explicit-limit (CFL) analysis
steps/                       one runner per step type (adsorption, blowdown, heating, cooling, pressurization)
postprocessing/              KPIs and balance tables, plots, the (r, theta, z, t) probe
validation/                  conservation checks, CSS test, Young-2021 comparison, regression baseline
analysis/                    one-at-a-time sensitivity analysis
simio/                       state CSV checkpoints, field exports
tests/                       unit tests
MODEL.md / PARAMETERS.md     equations and numerics / every parameter with its source
```

## Sensitivity analysis

```bash
python run_sensitivity.py                       # parameters and KPIs listed under `sensitivity:` in defaults.yaml
python run_sensitivity.py --set sensitivity.processes=4
```

Each parameter is perturbed ±10 % (or by an absolute step for temperatures and humidity)
around the base case on a reduced grid, every perturbed run warm-started from the base
case's cyclic steady state. Output: `results/<run name>_sensitivity/sensitivity.csv` with normalised
sensitivities S = ∂ln KPI/∂ln p (|S| = 1: a 10 % change moves the KPI by 10 %), and tornado
plots. These are *local* derivatives around one operating point.

## Results

### Validation against Young et al. 2021 (1D-equivalent case)

`config/cases/young2021_1d_equivalent.yaml` reproduces the assumptions of Young's 1D model
(one radial cell, lumped jacket U = 14 W/m²/K, no wall heat capacity, their dispersion
correlation, the six-step cycle of their minimum-work optimum). Cyclic steady state after
3 cycles, ~2 min; balances close to 4·10⁻⁷ (mass) and 4·10⁻⁶ (energy).
`python -m validation.young2021_compare` prints:

| quantity | this model | Young 2021 (ESI Table S10, WADST) | difference |
|---|---|---|---|
| CO2 purity (dry) | 99.49 % | 97.45 % | +2.0 points |
| working capacity | 0.984 mol/kg | 1.00 mol/kg | −1.6 % |
| specific heat duty | 8.92 MJ/kg | 9.19 MJ/kg | −2.9 % |
| productivity | 64.2 kg/m³/day | 64.4 kg/m³/day | −0.2 % |
| mechanical work (vacuum + blower) | 0.42 MJ/kg | 0.80 MJ/kg | see below |
| work equivalent | 1.94 MJ/kg | 2.37 MJ/kg | follows the line above |
| CO2 recovery | 69.1 % | not tabulated | |

Capacity, heat duty and productivity agree within 3 %. Two honest caveats:

* **Mechanical work** depends on a convention the paper's text does not settle (its
  equations are images): whether the vacuum pump compresses the water vapour. With
  `energy.pump_handles_water: false` (default: water condensed before compression) the
  case gives 0.42 MJ/kg; with `true` it gives 1.13 MJ/kg; Young's 0.80 lies in between.
  The thermal side is unaffected.
* **Purity** is 2 points higher, most likely because the blowdown here is pressure-driven
  through a finite-speed pump and removes more N2 before the product step than Young's
  prescribed pressure history. Treat purity differences of this size as model-form
  differences, not as validation failures.

### Default 3D case

`python run_dac_simulation.py` — 8 × 8 × 12 cells in (r, θ, z), d = 5 cm, L = 1 cm, full
rim jacket, Young's sorbent, feed and cycle timing. Cyclic steady state after 3 cycles,
~3.4 min per cycle on one CPU core; balances close to 7·10⁻⁶ (mass) and 4·10⁻⁶ (energy).

| KPI | value |
|---|---|
| CO2 recovered | 54.2 % |
| purity (dry / wet) | 99.2 % / 20.8 % |
| working capacity | 0.97 mol/kg |
| productivity | 63.3 kg/m³/day |
| specific thermal energy | 11.5 MJ/kg CO2 |
| specific electric energy | 0.47 MJ/kg CO2 (vacuum 0.20, blower 0.27) |
| work equivalent | 2.43 MJ/kg CO2 |
| water co-desorbed | 1.56 kg H2O/kg CO2 |

Where the heat goes (`energy_breakdown.csv`, heating step, per kg CO2): water desorption
4.0 MJ, bed sensible heat 3.2 MJ, **steel wall 2.5 MJ**, CO2 desorption 1.6 MJ, gas enthalpy
0.2 MJ. Two things to think about: water costs 2.5× more heat than CO2 itself, and for a
lab-scale bed the 2 mm steel wall has *more* heat capacity than the sorbent it holds — at
industrial scale that term disappears. Recovery is 54 % rather than Young's 69 % because
this cycle keeps Young's fixed 8000 s adsorption time, which is longer than this bed needs
(CO2 breaks through after ~1.2 h; see `plots/outlet_and_bed_history.png`) — shortening it,
or ending it with `stop_on`, is the first optimisation to try.

**A genuinely 3D case** — `config/cases/asymmetric_jacket.yaml` (jacket on half the
circumference, ±30 % feed maldistribution): recovery falls to 17 % and the heat duty rises
to 59 MJ/kg, of which 37 MJ/kg is lost to ambient — the steel wall conducts the jacket heat
around to the uncovered side, which loses it. At the end of heating, at r = 1.25 cm,
z = 5 mm: q_CO2 = 1.07 mol/kg and T = 354 K under the jacket (θ = 90°) but 1.45 mol/kg and
345 K on the other side (θ = 270°); the fully jacketed default reaches 0.57 mol/kg there.

## Numerics in one paragraph

Conservative finite volumes on an (r, θ, z) grid (the axis needs no special treatment:
the inner wedge faces have zero area). State per cell: gas concentrations, loadings and
volumetric internal energy; T and P are algebraic. Pressure diffuses on a ~10⁻⁷ s time
scale while CO2 fronts move over hours, so explicit methods would need ~10¹¹ steps per
adsorption step (`numerics/cfl.py`); the production solver is an L-stable SDIRK2 with
Newton–Krylov (GMRES on finite-difference Jacobian-vector products, preconditioned by a
sparse LU of a coloured finite-difference Jacobian), Newton-aware step control, and a
"conservation polish" that makes every stage close the global species and energy balances.
Details and the reasons behind each choice: MODEL.md §9.

## Honest notes — read before using a number

* **Deviation from the roadmap:** OPEN_VACUUM is implemented as a vacuum *pump* (it fixes
  the pressure-driven outflow), not as "C = 0 at the outlet". Imposing C = 0 creates an
  artificial dispersive sink and over-predicts desorption (MODEL.md §8).
* **Default geometry:** L = 1 cm as Young's flat bed but d = 5 cm instead of 10 cm. Heated
  only through the rim, a 10 cm disc is conduction-limited (centre ≈ 30 K below the wall
  during desorption); that case is kept as `config/cases/young2021_flat_bed_rim_jacket.yaml`
  and is worth running once to see why real contactors use internal heat exchangers
  (`config/cases/internal_hx.yaml`).
* **No condensation model:** relative humidity fed to the water isotherm is capped at 0.97.
* **The LDF constants are lumped fits** to adsorption breakthroughs; desorption kinetics at
  100 °C are extrapolated. They are the most uncertain inputs — run the sensitivity analysis.
* **Literature errata** found while sourcing parameters (Zhang et al. 2025, Table A.1 and
  eq. A.1; the project's equation notes) are listed in MODEL.md §13; the code uses the
  corrected values.
* **Numerical lessons** that cost real time and are documented in MODEL.md: a
  velocity-dependent dispersion coefficient evaluated with the pressure-driven face velocity
  makes the discrete system unstable; a blended valve stream is anti-dissipative near zero
  flow; step-size control must respect Newton, not only the error estimate.
