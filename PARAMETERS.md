# PARAMETERS.md — every parameter you can change, its default, and where it comes from

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

| key | default | unit | allowed | meaning | source |
|---|---|---|---|---|---|
| `run.name` | `default_3d` | - | str | run name; outputs go to <results_dir>/<name>/ | modelling assumption / user choice |
| `run.results_dir` | `results` | - | str | output root directory (relative to project root) | modelling assumption / user choice |
| `run.max_cycles` | `12` | - | [1, 500] | maximum number of TVSA cycles | modelling assumption / user choice |
| `run.min_cycles` | `2` | - | [1, 500] | minimum number of cycles before CSS may be declared | modelling assumption / user choice |
| `run.css_tolerance` | `0.001` | - | [1e-08, 1.0] | cyclic-steady-state tolerance on normalised end-of-cycle state change | modelling assumption / user choice |
| `run.css_mb_tolerance` | `0.005` | - | [1e-08, 1.0] | CSS tolerance on the cycle CO2 balance \|in-out\|/in | Young 2021 used 0.5 % |
| `run.handoff_via_csv` | `True` | - | bool | each step reloads the state CSV written by the previous step | modelling assumption / user choice |
| `run.snapshots_per_step` | `12` | - | [2, 10000] | number of field snapshots stored per step | modelling assumption / user choice |
| `run.write_fields_every_cycle` | `False` | - | bool | write field CSVs for every cycle (else only the last) | modelling assumption / user choice |
| `run.initial_state` | `clean` | - | clean / feed_equilibrium / csv | initial bed: clean sorbent in feed air, saturated with feed, or loaded from CSV | modelling assumption / user choice |
| `run.initial_state_csv` | `None` | - | str | state CSV used when initial_state == csv | modelling assumption / user choice |
| `run.make_plots` | `True` | - | bool | write PNG plots at the end of the run | modelling assumption / user choice |
| `components.water` | `True` | - | bool | include H2O (humid air, co-adsorption). false = dry CO2/N2 model | modelling assumption / user choice |
| `geometry.bed_length` | `0.01` | m | [0.0001, 50.0] | packed length L (flow direction z) | Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9 |
| `geometry.bed_radius` | `0.025` | m | [0.0001, 10.0] | inner radius R of the bed | design choice: d = 5 cm; Young 2021 ESI S9 flat bed has d = 10 cm, which is conduction-limited with a rim jacket (cases/young2021_flat_bed_rim_jacket.yaml) |
| `geometry.wall_thickness` | `0.002` | m | [1e-05, 0.5] | metal wall thickness | modelling assumption / user choice |
| `grid.Nr` | `8` | - | [1, 200] | radial finite volumes | modelling assumption / user choice |
| `grid.Ntheta` | `8` | - | [1, 360] | angular finite volumes (1 = axisymmetric 2D) | modelling assumption / user choice |
| `grid.Nz` | `12` | - | [2, 2000] | axial finite volumes | modelling assumption / user choice |
| `grid.radial_spacing` | `uniform` | - | uniform / equal_area | radial face placement | modelling assumption / user choice |
| `sorbent.name` | `Lewatit VP OC 1065` | - | str | label only | modelling assumption / user choice |
| `sorbent.particle_diameter` | `0.00052` | m | [1e-06, 0.05] | bead diameter d_p | Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9 (R_p = 0.26 mm) |
| `sorbent.particle_density` | `880.0` | kg/m3 | [50, 5000] | pellet density rho_p | Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9 |
| `sorbent.bed_voidage` | `0.4` | - | [0.2, 0.95] | inter-particle voidage eps_b | Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9 |
| `sorbent.total_voidage` | `0.54` | - | [0.2, 0.99] | total voidage eps_t = eps_b + eps_p(1-eps_b) | Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9 |
| `sorbent.cp_solid` | `1580.0` | J/kg/K | [100, 5000] | sorbent heat capacity | Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9 |
| `sorbent.particle_conductivity` | `0.16` | W/m/K | [0.01, 50] | solid (bead) thermal conductivity | Zhang et al. 2025, Energy 334:137733, Table 1 (Galanti 2025 SI uses 0.43) |
| `sorbent.max_temperature` | `393.15` | K | [300, 700] | temperature above which a warning is logged | Young 2021 (max 100 C regeneration for oxidative stability) |
| `isotherm.co2.q_inf0` | `4.86` | mol/kg | [0.01, 50] | Toth saturation capacity at T0 | Young et al. 2021 EES, Table 2 (Toth fit, Lewatit VP OC 1065) |
| `isotherm.co2.chi` | `0.0` | - | [-50, 50] | Toth capacity T-dependence | Young et al. 2021 EES, Table 2 (Toth fit, Lewatit VP OC 1065) |
| `isotherm.co2.b0` | `2.85e-21` | 1/Pa | [0, 10000000000.0] | Toth affinity pre-factor (b = b0 exp(-dH0/RT)) | Young et al. 2021 EES, Table 2 (Toth fit, Lewatit VP OC 1065) |
| `isotherm.co2.dH0` | `-117798.0` | J/mol | [-500000.0, 0] | Toth affinity energy (negative) | Young et al. 2021 EES, Table 2 (Toth fit, Lewatit VP OC 1065) |
| `isotherm.co2.tau0` | `0.209` | - | [0.01, 5] | Toth heterogeneity at T0 | Young et al. 2021 EES, Table 2 (Toth fit, Lewatit VP OC 1065) |
| `isotherm.co2.alpha` | `0.523` | - | [-10, 10] | Toth heterogeneity T-dependence | Young et al. 2021 EES, Table 2 (Toth fit, Lewatit VP OC 1065) |
| `isotherm.co2.T0` | `298.15` | K | [200, 500] | Toth reference temperature | Young et al. 2021 EES, Table 2 (Toth fit, Lewatit VP OC 1065) |
| `isotherm.co2.p_regularization` | `0.01` | Pa | [0.0, 10.0] | trace-level regularisation q* -> q* p/(p+p_reg); 0 disables (numerically harsh with WADST) | modelling assumption / user choice |
| `isotherm.h2o.q_m` | `3.63` | mol/kg | [0.01, 100] | GAB monolayer capacity | Young et al. 2021 EES, Table 3 (GAB fit) |
| `isotherm.h2o.C` | `47110.0` | J/mol | [0, 200000.0] | GAB: E1 = C - exp(D T) | Young et al. 2021 EES, Table 3 (GAB fit) |
| `isotherm.h2o.D` | `0.023744` | 1/K | [0, 1] | GAB: E1 = C - exp(D T) | Young et al. 2021 EES, Table 3 (GAB fit) |
| `isotherm.h2o.F` | `57706.0` | J/mol | [0, 200000.0] | GAB: E2-9 = F + G T | Young et al. 2021 EES, Table 3 (GAB fit) |
| `isotherm.h2o.G` | `-47.814` | J/mol/K | [-1000, 1000] | GAB: E2-9 = F + G T | Young et al. 2021 EES, Table 3 (GAB fit) |
| `isotherm.h2o.E10_slope` | `-44.38` | J/mol/K | [-1000, 1000] | E10+ = slope*T + intercept | Young 2021 eq. 8 |
| `isotherm.h2o.E10_intercept` | `57220.0` | J/mol | [0, 200000.0] | E10+ = slope*T + intercept | Young 2021 eq. 8 |
| `isotherm.h2o.rh_cap` | `0.97` | - | [0.5, 0.999] | numerical cap on relative humidity passed to GAB (no condensation model) | modelling assumption / user choice |
| `isotherm.coadsorption.model` | `wadst` | - | none / wadst / mechanistic / stampi_bombelli | effect of adsorbed water on CO2 isotherm | Young 2021 EES |
| `isotherm.coadsorption.water_loading_basis` | `state` | - | state / equilibrium | use local adsorbed water q_H2O (state) or q*_H2O in the CO2 isotherm | modelling assumption / user choice |
| `isotherm.coadsorption.wadst.b0_wet` | `1.23e-18` | 1/Pa | [0, 10000000000.0] | WADST wet-site affinity pre-factor | Young et al. 2021 EES, Table 4 (co-adsorption fits) |
| `isotherm.coadsorption.wadst.q_inf0_wet` | `9.035` | mol/kg | [0, 100] | WADST wet-site capacity | Young et al. 2021 EES, Table 4 (co-adsorption fits) |
| `isotherm.coadsorption.wadst.tau0_wet` | `0.053` | - | [0.001, 5] | WADST wet-site Toth tau0 | Young et al. 2021 EES, Table 4 (co-adsorption fits) |
| `isotherm.coadsorption.wadst.chi_wet` | `0.0` | - | [-50, 50] | WADST wet-site chi | Young et al. 2021 EES, Table 4 (co-adsorption fits) |
| `isotherm.coadsorption.wadst.alpha_wet` | `0.053` | - | [-10, 10] | WADST wet-site alpha | Young et al. 2021 EES, Table 4 (co-adsorption fits) |
| `isotherm.coadsorption.wadst.dH_wet` | `-203687.0` | J/mol | [-500000.0, 0] | WADST wet-site affinity energy | Young et al. 2021 EES, Table 4 (co-adsorption fits) |
| `isotherm.coadsorption.wadst.A` | `1.532` | mol/kg | [1e-06, 100] | WADST critical water loading | Young et al. 2021 EES, Table 4 (co-adsorption fits) |
| `isotherm.coadsorption.mechanistic.f_blocked_max` | `0.433` | - | [0, 1] | max blocked amine fraction | Young et al. 2021 EES, Table 4 (co-adsorption fits) |
| `isotherm.coadsorption.mechanistic.k` | `0.795` | kg/mol | [0, 100] | Avrami blocking constant | Young et al. 2021 EES, Table 4 (co-adsorption fits) |
| `isotherm.coadsorption.mechanistic.phi_dry` | `1.0` | - | [0.01, 2] | dry amine efficiency | Young et al. 2021 EES, Table 4 (co-adsorption fits) |
| `isotherm.coadsorption.mechanistic.A` | `1.535` | mol/kg | [1e-06, 100] | critical water loading | Young et al. 2021 EES, Table 4 (co-adsorption fits) |
| `isotherm.coadsorption.mechanistic.dH_wet` | `-130155.0` | J/mol | [-500000.0, 0] | wet affinity energy | Young et al. 2021 EES, Table 4 (co-adsorption fits) |
| `isotherm.coadsorption.mechanistic.n` | `1.425` | - | [0.1, 10] | Avrami exponent | Young et al. 2021 EES, Table 4 (co-adsorption fits) |
| `isotherm.coadsorption.stampi_bombelli.gamma` | `-0.137` | kg/mol | [-10, 10] | capacity factor | Young et al. 2021 EES, Table 4 (co-adsorption fits) |
| `isotherm.coadsorption.stampi_bombelli.beta` | `5.612` | kg/mol | [-10, 100] | affinity factor | Young et al. 2021 EES, Table 4 (co-adsorption fits) |
| `kinetics.k_ldf_co2` | `0.003` | 1/s | [1e-07, 100] | LDF coefficient CO2 at T_ref | Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9 (breakthrough fit) |
| `kinetics.k_ldf_h2o` | `0.0086` | 1/s | [1e-07, 100] | LDF coefficient H2O at T_ref | Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9 (breakthrough fit) |
| `kinetics.Ea_co2` | `0.0` | J/mol | [0, 200000.0] | Arrhenius activation energy of k_LDF,CO2 (0 = constant) | Young 2021 / Galanti 2025 use constant k |
| `kinetics.Ea_h2o` | `0.0` | J/mol | [0, 200000.0] | Arrhenius activation energy of k_LDF,H2O | modelling assumption / user choice |
| `kinetics.T_ref` | `298.15` | K | [200, 500] | reference temperature of k_LDF | modelling assumption / user choice |
| `thermo.dH_ads_co2` | `-70000.0` | J/mol | [-300000.0, 0] | heat of adsorption CO2 (energy balance) | Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9 |
| `thermo.dH_ads_h2o` | `-46000.0` | J/mol | [-300000.0, 0] | heat of adsorption H2O (energy balance) | Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9 |
| `thermo.cp_ads_co2` | `37.1` | J/mol/K | [0, 500] | adsorbed-phase heat capacity CO2 | modelling assumption / user choice (gas-like) |
| `thermo.cp_ads_h2o` | `75.3` | J/mol/K | [0, 500] | adsorbed-phase heat capacity H2O | modelling assumption / user choice (liquid water 75.3; gives Kirchhoff T-dependence of dH) |
| `thermo.T_ref` | `298.15` | K | [200, 500] | reference temperature of enthalpies / dH_ads | modelling assumption / user choice |
| `properties.backend` | `ideal_gas` | - | ideal_gas / coolprop | gas property backend (coolprop needs `pip install CoolProp`) | modelling assumption / user choice |
| `properties.cp_reference_T` | `330.0` | K | [200, 600] | temperature at which constant gas cp are evaluated for the energy balance | modelling assumption / user choice |
| `transport.dispersion_model` | `edwards_richardson` | - | edwards_richardson / young2021 / wakao_kagei | axial dispersion correlation | Edwards & Richardson 1968, Chem. Eng. Sci. 23:109 |
| `transport.er_gamma1` | `0.73` | - | [0, 5] | E-R tortuosity factor gamma1 | Edwards & Richardson 1968, Chem. Eng. Sci. 23:109 |
| `transport.er_gamma2` | `0.5` | - | [0, 5] | E-R gamma2 | Edwards & Richardson 1968, Chem. Eng. Sci. 23:109 |
| `transport.er_beta` | `13.0` | - | [0, 100] | E-R beta (beta*gamma1 = 9.49) | Edwards & Richardson 1968, Chem. Eng. Sci. 23:109 |
| `transport.Pe_mass_radial_inf` | `10.0` | - | [1, 100] | asymptotic transverse mass Peclet number | Delgado 2006 review (8-12) |
| `transport.Pe_heat_radial_inf` | `10.0` | - | [1, 100] | asymptotic radial thermal Peclet number | Galanti et al. 2025, RSC Energy Adv. (doi:10.1039/D5YA00336A), SI Tables S7-S8 |
| `transport.Pe_heat_axial_inf` | `2.0` | - | [0.1, 100] | asymptotic axial thermal Peclet number | Yagi & Kunii (0.5 Re Pr) |
| `transport.bed_conductivity_model` | `zehner_schlunder` | - | zehner_schlunder / parallel / constant | stagnant effective bed conductivity | Zehner & Schlünder 1970, Chem.-Ing.-Tech. 42:933 (VDI Heat Atlas D6.3) |
| `transport.bed_conductivity_constant` | `0.12` | W/m/K | [0.001, 50] | used if model == constant | modelling assumption / user choice |
| `transport.zs_shape_factor` | `1.25` | - | [0.5, 5] | Zehner-Schluender shape factor (1.25 spheres) | Zehner & Schlünder 1970, Chem.-Ing.-Tech. 42:933 (VDI Heat Atlas D6.3) |
| `transport.wall_htc_model` | `martin_nilles` | - | martin_nilles / constant | bed-to-wall heat transfer coefficient | Martin & Nilles 1993, Chem.-Ing.-Tech. 65:1468 (as used in G25 SI) |
| `transport.wall_htc_constant` | `300.0` | W/m2/K | [0.1, 100000.0] | used if wall_htc_model == constant | modelling assumption / user choice |
| `transport.ergun_A` | `150.0` | - | [0, 10000.0] | Ergun viscous constant | modelling assumption / user choice |
| `transport.ergun_B` | `1.75` | - | [0, 100] | Ergun inertial constant | modelling assumption / user choice |
| `transport.heat_transfer_mode` | `resolved` | - | resolved / lumped | resolved: conduction resistances in the bed; lumped: 1D-style overall U at the wall | modelling assumption / user choice |
| `transport.lumped_U` | `14.0` | W/m2/K | [0.01, 100000.0] | overall coefficient for lumped mode | Young 2021 ESI Table S9 (h_j = 14) |
| `wall.enabled` | `True` | - | bool | resolve wall heat capacity and conduction | modelling assumption / user choice |
| `wall.density` | `7800.0` | kg/m3 | [100, 30000] | wall density | Galanti et al. 2025, RSC Energy Adv. (doi:10.1039/D5YA00336A), SI Tables S7-S8 |
| `wall.cp` | `513.0` | J/kg/K | [100, 5000] | wall heat capacity | Galanti et al. 2025, RSC Energy Adv. (doi:10.1039/D5YA00336A), SI Tables S7-S8 |
| `wall.conductivity` | `16.0` | W/m/K | [0.01, 500] | wall conductivity | Galanti et al. 2025, RSC Energy Adv. (doi:10.1039/D5YA00336A), SI Tables S7-S8 |
| `jacket.h_jacket` | `26.0` | W/m2/K | [0, 100000.0] | jacket fluid to wall heat transfer coefficient | Galanti et al. 2025, RSC Energy Adv. (doi:10.1039/D5YA00336A), SI Tables S7-S8 (h_ext = 26) |
| `jacket.theta_coverage_deg` | `[0.0, 360.0]` | deg | list | [start, end] angular extent of the jacket (0-360) | modelling assumption / user choice |
| `jacket.z_coverage` | `[0.0, 1.0]` | - | list | [start, end] axial extent of the jacket as fractions of L | modelling assumption / user choice |
| `jacket.h_ambient` | `5.0` | W/m2/K | [0, 10000.0] | uncovered wall to ambient coefficient | modelling assumption / user choice |
| `jacket.ramp_time` | `0.0` | s | [0, 100000.0] | first-order lag of jacket fluid temperature (0 = step) | modelling assumption / user choice |
| `internal_hx.enabled` | `False` | - | bool | embedded volumetric heat exchanger | modelling assumption / user choice |
| `internal_hx.hA_per_volume` | `2000.0` | W/m3/K | [0, 10000000.0] | volumetric hA of embedded exchanger | modelling assumption / user choice |
| `internal_hx.follows_jacket` | `True` | - | bool | exchanger fluid temperature = jacket temperature | modelling assumption / user choice |
| `feed.T` | `288.15` | K | [200, 400] | feed temperature | Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9 |
| `feed.P` | `101325.0` | Pa | [10000.0, 1000000.0] | feed pressure (also outlet pressure in adsorption) | Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9 |
| `feed.RH` | `0.55` | - | [0, 1] | feed relative humidity | Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9 |
| `feed.y_co2_dry` | `0.0004` | - | [0, 1] | CO2 mole fraction in dry air | Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9 |
| `feed.velocity` | `0.0706` | m/s | [1e-06, 10] | superficial feed velocity at feed T,P | Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9 |
| `feed.profile` | `uniform` | - | uniform / cosine_theta / linear_r / central_disc / sector | inlet flow distribution (maldistribution / partial-inlet studies) | modelling assumption / user choice |
| `feed.profile_amplitude` | `0.0` | - | [0, 0.99] | amplitude of the inlet maldistribution | modelling assumption / user choice |
| `feed.profile_theta0_deg` | `0.0` | deg | [0, 360] | angle of maximum flow (cosine_theta, sector) | modelling assumption / user choice |
| `feed.profile_radius_fraction` | `0.5` | - | [0.05, 1.0] | inlet radius / R for central_disc | modelling assumption / user choice |
| `feed.profile_sector_deg` | `90.0` | deg | [1.0, 360.0] | angular width of the inlet for sector | modelling assumption / user choice |
| `ambient.T` | `288.15` | K | [200, 400] | ambient temperature (vent back-flow, heat loss) | modelling assumption / user choice |
| `ambient.P` | `101325.0` | Pa | [10000.0, 1000000.0] | ambient pressure (atmospheric vent, pump discharge) | modelling assumption / user choice |
| `ambient.RH` | `0.55` | - | [0, 1] | ambient relative humidity (vent back-flow) | modelling assumption / user choice |
| `ambient.y_co2_dry` | `0.0004` | - | [0, 1] | ambient CO2 dry mole fraction | modelling assumption / user choice |
| `operation.P_vacuum` | `2.0e4` | Pa | [100, 100000.0] | vacuum level (blowdown / desorption) | Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9 |
| `operation.T_heat` | `373.15` | K | [273, 600] | jacket temperature in heating | Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9 |
| `operation.T_cool` | `288.15` | K | [200, 600] | jacket temperature in cooling | Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9 |
| `operation.T_adsorption_jacket` | `288.15` | K | [200, 600] | jacket temperature in adsorption | Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9 |
| `operation.valve_time_constant` | `5.0` | s | [0.001, 10000.0] | time constant of the exponential pressure ramp at an opened valve | modelling assumption / user choice |
| `operation.vent_model` | `valve` | - | valve / ramp | atmospheric vent / feed header: finite-conductance valve or imposed exponential P ramp | modelling assumption / user choice |
| `operation.vent_conductance` | `500.0` | m/s | [1e-06, 10000.0] | vent valve conductance per unit cross-section, N = k_v (P_res - P)/(R T) | modelling assumption / user choice |
| `operation.vacuum_model` | `pump` | - | pump / ramp | vacuum port: finite-speed pump (P relaxes to P_vacuum) or imposed exponential P ramp | modelling assumption / user choice |
| `operation.pump_speed_per_area` | `0.005` | m/s | [1e-07, 100.0] | vacuum pump volumetric speed per unit bed cross-section (S/A); sets the blowdown time constant eps_t L/(S/A) | modelling assumption / user choice |
| `operation.vacuum_check_valve` | `True` | - | bool | prevent back-flow from vacuum lines | modelling assumption / user choice |
| `ports.V1` | `top` | - | top / bottom | end of the column where V1 (inlet) sits | modelling assumption / user choice |
| `ports.V2` | `bottom` | - | top / bottom | end of V2 (vent / blowdown) | modelling assumption / user choice |
| `ports.V3` | `bottom` | - | top / bottom | end of V3 (product) | Young 2021 (z = L) |
| `energy.eta_vacuum` | `0.7` | - | [0.01, 1] | vacuum pump isentropic efficiency | Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9 |
| `energy.vacuum_efficiency_model` | `constant` | - | constant / krishnamurthy | constant, or 0.8*19.55p/(1+19.55p) (p in bar) as in G25 / Krishnamurthy et al. 2014 | modelling assumption / user choice |
| `energy.eta_blower` | `0.5` | - | [0.01, 1] | feed blower efficiency | Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9 |
| `energy.gamma_gas` | `1.4` | - | [1.0, 1.7] | cp/cv for compression work | Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9 |
| `energy.pump_handles_water` | `False` | - | bool | include H2O in the vacuum-pump molar flow (false: water condensed before compression) | assumption; Young 2021 condense water 'in the compression process' without giving the work convention |
| `energy.eta_turbine` | `0.75` | - | [0, 1] | efficiency for heat -> work equivalent | Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9 |
| `energy.T_low_equiv` | `288.15` | K | [200, 400] | T_L in work equivalent | Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9 |
| `solver.method` | `trbdf2` | - | trbdf2 / scipy_bdf / rk4 | time integrator: implicit Newton TR-BDF2 (default), scipy BDF cross-check, explicit RK4 | modelling assumption / user choice |
| `solver.scheme` | `sdirk2` | - | sdirk2 / trbdf2 | 2-stage implicit RK used by the Newton integrator: L-stable SDIRK2 or TR-BDF2 | modelling assumption / user choice |
| `solver.rtol` | `0.0001` | - | [1e-10, 0.1] | relative tolerance | modelling assumption / user choice |
| `solver.atol_c_co2` | `1e-06` | mol/m3 | [1e-14, 1] | absolute tolerance CO2 concentration | modelling assumption / user choice |
| `solver.atol_c_h2o` | `1e-05` | mol/m3 | [1e-14, 1] | absolute tolerance H2O concentration | modelling assumption / user choice |
| `solver.atol_c_inert` | `0.0001` | mol/m3 | [1e-14, 10] | absolute tolerance inert concentration | modelling assumption / user choice |
| `solver.atol_q` | `1e-06` | mol/kg | [1e-14, 1] | absolute tolerance loadings | modelling assumption / user choice |
| `solver.atol_T` | `0.001` | K | [1e-08, 10] | absolute tolerance temperature | modelling assumption / user choice |
| `solver.h_initial` | `0.0001` | s | [1e-10, 100] | initial time step | modelling assumption / user choice |
| `solver.h_max_fraction` | `0.02` | - | [1e-05, 1] | max time step as a fraction of step duration | modelling assumption / user choice |
| `solver.max_newton_iter` | `7` | - | [2, 50] | Newton iterations per stage | modelling assumption / user choice |
| `solver.preconditioner` | `theta_block` | - | theta_block / full | LU used by Newton: full Jacobian, or theta-sector-decoupled (exact for axisymmetric states) | modelling assumption / user choice |
| `solver.precond_theta_rings` | `3` | - | [0, 1000] | number of innermost rings whose theta couplings are kept in the theta_block preconditioner | modelling assumption / user choice |
| `solver.refactor_ratio` | `3.0` | - | [1.05, 100.0] | h change factor that triggers re-factorisation | modelling assumption / user choice |
| `solver.gmres_refresh` | `20` | - | [1, 100000] | GMRES iterations per step that trigger a new Jacobian | modelling assumption / user choice |
| `solver.newton_mode` | `krylov` | - | hybrid / krylov / modified | stage solver: modified Newton with stored LU, Newton-Krylov (GMRES, LU-preconditioned), or both | modelling assumption / user choice |
| `solver.newton_kappa` | `0.01` | - | [0.0001, 1] | Newton convergence safety factor | modelling assumption / user choice |
| `solver.conservation_tol` | `0.001` | - | [0.0, 1.0] | stage balance tolerance: \|W R\| <= tol * W (atol + rtol\|y\|) per species and energy (0 = off) | modelling assumption / user choice |
| `solver.jacobian_max_age` | `40` | - | [1, 10000] | accepted steps before forced Jacobian refresh | modelling assumption / user choice |
| `solver.limiter` | `van_albada` | - | van_albada / van_albada_tvd / van_leer / minmod / upwind | axial TVD flux limiter | modelling assumption / user choice |
| `solver.max_steps_per_step` | `200000` | - | [10, 100000000] | safety cap on time steps per cycle step | modelling assumption / user choice |
| `solver.outlet_record_every` | `1` | - | [1, 1000000] | record outlet time series every n steps | modelling assumption / user choice |
| `sensitivity.relative_step` | `0.1` | - | [0.0001, 0.9] | relative perturbation for sensitivities | modelling assumption / user choice |
| `sensitivity.cycles` | `6` | - | [1, 200] | cycles per perturbed run | modelling assumption / user choice |
| `sensitivity.grid` | `{'Nr': 6, 'Ntheta': 1, 'Nz': 10}` | - | any | grid override for sensitivity runs {Nr, Ntheta, Nz} | modelling assumption / user choice |
| `sensitivity.processes` | `2` | - | [1, 256] | parallel worker processes | modelling assumption / user choice |
| `sensitivity.parameters` | `['kinetics.k_ldf_co2', 'kinetics.k_ldf_h2o', 'thermo.dH_ads_co2', 'thermo.dH_ads_h2o', 'sorbent.cp_solid', 'sorbent.particle_conductivity', 'jacket.h_jacket', 'operation.T_heat', 'operation.P_vacuum', 'feed.velocity', 'feed.RH']` | - | list | dotted parameter keys to perturb | modelling assumption / user choice |
| `sensitivity.kpis` | `['recovery_pct', 'purity_dry_pct', 'productivity_kg_m3_day', 'specific_thermal_kJ_per_kg', 'specific_electric_kJ_per_kg', 'specific_total_kJ_per_kg']` | - | list | KPI names to report | modelling assumption / user choice |
