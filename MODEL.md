# MODEL.md — equations, boundary conditions and numerics

This is the technical reference of `dac_vtsa`. Every equation below is implemented
in the module named next to it; the parameter names in `monospace` are the keys
of `config/defaults.yaml` (all documented with unit, range and source in
`PARAMETERS.md`).

---------------------------------------------------------------------------------

## 1. What is being modelled

A packed bed of amine sorbent (Lewatit VP OC 1065) in a cylindrical, jacketed
vessel, operated as a temperature–vacuum swing adsorption (TVSA) cycle on humid
air. Default cycle (Young et al. 2021, *Energy Environ. Sci.* 14, 5377):

| step | V1 (top, feed) | V2 (bottom) | V3 (product, bottom) | jacket | what happens |
|---|---|---|---|---|---|
| adsorption (8000 s) | OPEN (blower) | OPEN_ATMOSPHERIC | CLOSED | 288 K | air in, CO2 + H2O adsorb, N2 vented |
| blowdown (60 s) | CLOSED | OPEN_VACUUM | CLOSED | 288 K | gas phase pumped out (mostly N2) |
| heating (22 400 s) | CLOSED | CLOSED | OPEN_VACUUM | 373 K | CO2 + H2O desorb into the vacuum line = product |
| cooling (400 s) | CLOSED | CLOSED | CLOSED | 288 K | bed cools before the next adsorption |

Assumptions (standard for this class of model; Young 2021, Galanti 2025, Zhang 2025):

* ideal gas; three gas species CO2, H2O, N2 (N2 represents all inert air: N2 + O2 + Ar);
  `components.water: false` gives a dry CO2/N2 model;
* pseudo-homogeneous energy balance (gas and solid at the same local T);
* linear-driving-force (LDF) uptake, CO2 and H2O adsorb, N2 is inert;
* Ergun momentum balance (no inertia term in time — pressure waves are quasi-steady);
* no condensation of water in the bed (RH fed to GAB is capped at `isotherm.h2o.rh_cap`);
* thin metal wall resolved in (θ, z), lumped in r; jacket = fluid at a set temperature.

## 2. Geometry and grid — `core/grid.py`

Cylinder of radius R = `geometry.bed_radius`, length L = `geometry.bed_length`.
**z = 0 is the top (feed end, V1), z = L the bottom (V2/V3).** Conservative finite
volumes on an (r, θ, z) tensor grid of Nr × Nθ × Nz cells
(`grid.Nr`, `grid.Ntheta`, `grid.Nz`), uniform or equal-area rings
(`grid.radial_spacing`).

* Cell (i, j, k): r ∈ [r_i, r_{i+1}], θ ∈ [θ_j, θ_{j+1}], z ∈ [z_k, z_{k+1}].
  Volume V_i = ½(r²_{i+1} − r²_i) Δθ Δz.
* The innermost ring is a set of wedges whose inner face has **zero area**: no flux
  crosses r = 0 and no 1/r term ever appears (the axis singularity of
  finite-difference cylindrical schemes does not exist here).
* θ is periodic. Nθ = 1 gives the axisymmetric 2D (r, z) model; Nr = Nθ = 1 the 1D model.
* Wall ring (thickness `geometry.wall_thickness`) discretised in (θ, z).

## 3. State vector

Per cell: gas concentrations c_k [mol/m³] (k = CO2, H2O, N2), loadings q_j [mol/kg]
(j = CO2, H2O), and the volumetric internal energy e [J/m³]; per wall cell: T_w.
Temperature and pressure are **algebraic** functions of the state:

    P = (Σ_k c_k) R T

and T from e (section 6). Boundary flows and heat duties are integrated alongside
the state as *accumulators* (they are outputs, not unknowns of the Newton system).

## 4. Gas-phase species balances — `physics/mass.py`, `core/residual.py`

For every gas species k and every cell (exact flux form, so the sum over cells
telescopes and species inventories change only through the boundaries):

    ε_t V dc_k/dt = − Σ_faces (N_k · n) A  −  ρ_b V dq_k/dt
    N_k = u c_k − ε_b D C ∇y_k                      (convection + dispersion)

* ε_t = `sorbent.total_voidage`, ε_b = `sorbent.bed_voidage`,
  ρ_b = (1 − ε_b) ρ_p (`sorbent.particle_density`).
* Dispersion acts on **mole-fraction** gradients with one coefficient D for all
  species, so Σ_k J_k = 0 exactly: dispersion never creates net molar flow.
* Face velocities u in z, r and θ come from Ergun (section 5). Convection is
  MUSCL with a limiter (`solver.limiter`, default unclipped van Albada) along z,
  first-order upwind in r and θ (transverse velocities are small).
* Upwinding uses a smooth switch w(u) = ½(1 + tanh(u/10⁻⁴ m/s)) instead of a hard
  `if u > 0`: continuous derivatives, which Newton needs wherever a velocity changes
  sign. Below 10⁻⁴ m/s the face value blends to the central average, which is accurate
  and stable there (cell Péclet number < 0.1: dispersion dominates).
* **Dispersion coefficients use the nominal (feed) molar flux, not the instantaneous face
  velocity.** D_L = γ₁D_m + γ₂d_p v with v = N_nom/(C ε_b) at local conditions (N_nom =
  feed molar flux of the column in feed steps, 0 otherwise); the same for the axial thermal
  dispersion. Evaluating D with the pressure-driven face velocity couples D to ∇P, whose
  response is extremely stiff (∂u/∂P ~ 1/(a Δz)). That coupling has no physical meaning — the
  correlation describes mixing in the mean flow, not in microsecond pressure transients — but
  it gives the semi-discrete system spurious **unstable** eigenvalues (we found 29 modes with
  Re λ ≈ +1.2·10³ s⁻¹ at a near-stagnant open end), after which no time integrator behaves.
  Worth remembering as a dynamical-systems lesson: a constitutive law that is harmless in the
  continuum can destabilise the discretised system through a stiff fast variable.

## 5. Momentum — Ergun in every direction — `physics/pressure.py`

    −∂P/∂n = a u + b |u| u,   a = A μ (1−ε_b)²/(ε_b³ d_p²),   b = B ρ (1−ε_b)/(ε_b³ d_p)

(A = `transport.ergun_A` = 150, B = `transport.ergun_B` = 1.75), inverted per face as
u = −2g/(a + √(a² + 4b|g|)) with g = ∂P/∂n. This makes pressure a *diffusing*
field with diffusivity D_P ≈ K P/(ε_t μ) ~ 1–4 m²/s — the stiffest process in the
model (section 9).

## 6. Energy — conservative internal-energy form — `physics/energy.py`

State variable e = internal energy per unit bed volume:

    e = ε_t Σ_k c_k (cp_k (T − T_ref) − R T)
      + ρ_b [ cp_s (T − T_ref) + Σ_j q_j (cp_a,j (T − T_ref) + ΔH_j) ]

    V de/dt = − Σ_faces [ Σ_k h_k(T_f) N_k − λ ∇T ] · n A + Q_wall + Q_hx,   h_k = cp_k (T − T_ref)

* The heat of adsorption is **not** a source term: an adsorbed mole carries
  enthalpy h_k + ΔH_j (ΔH_j < 0), so adsorption releases heat automatically;
  pressure work and the (−ΔH − RT) constant-volume correction are contained in e.
* e is linear in T at fixed (c, q) ⇒ T = (e − ρ_b Σ q_j ΔH_j + A T_ref)/(A − ε_t R Σc_k),
  A = ε_t Σ c_k cp_k + ρ_b (cp_s + Σ q_j cp_a,j).
* Because e is conserved by construction, the energy balance reported in
  `energy_balance.csv` is a real check of the solver (it closes to ~1e-5–1e-7 relative).
* `thermo.cp_ads_h2o` = 75.3 J/mol/K (liquid-like) makes ΔH_H2O(T) follow Kirchhoff's law.

Conductivity: stagnant λ₀ from Zehner–Schlünder (`transport.bed_conductivity_model`,
particle conductivity `sorbent.particle_conductivity`), plus dispersive terms
C cp |u| d_p / Pe_h (`transport.Pe_heat_axial_inf`, `transport.Pe_heat_radial_inf`), with |u| from
the nominal molar flux (§4).

### Wall, jacket, ambient, internal exchanger

    C_w dT_w/dt = Q_jacket + Q_ambient − Q_wall→bed + (axial/θ conduction in the wall)
    Q_jacket   = cov(θ,z) h_jacket A_out (T_jacket − T_w)
    Q_ambient  = (1 − cov) h_ambient A_out (T_amb − T_w)
    Q_wall→bed = A_in (T_w − T_N) / (1/α_w + (R − r_N)/λ_r)

α_w from Martin & Nilles (`transport.wall_htc_model`), h_jacket = `jacket.h_jacket`.
The jacket can cover only part of the wall (`jacket.theta_coverage_deg`,
`jacket.z_coverage`) — this is what makes a run genuinely 3D. Optional volumetric
internal exchanger (`internal_hx.*`). `transport.heat_transfer_mode: lumped`
replaces bed-side resistance by a single U (`transport.lumped_U`) for 1D-style
comparisons with Young 2021.

## 7. Adsorption equilibrium and kinetics — `physics/isotherms.py`, `physics/kinetics.py`

LDF: dq_j/dt = k_j(T) (q*_j − q_j), k_j(T) = k_ref exp(−Ea/R (1/T − 1/T_ref))
(`kinetics.*`; Young 2021 fitted Ea = 0).

CO2 (dry), temperature-dependent Toth (Young 2021 eqs 1–4, `isotherm.co2.*`):

    q* = q_∞(T) b(T) p / (1 + (b p)^τ)^(1/τ),   b = b₀ exp(−ΔH₀/RT),
    q_∞ = q_∞0 exp(χ(1 − T/T₀)),  τ = τ₀ + α(1 − T₀/T)

H2O: GAB with Anderson temperature dependence (Young eqs 5–8, `isotherm.h2o.*`).

Co-adsorption (`isotherm.coadsorption.model`): `wadst` (default, weighted dual-site
Toth, Young eq. 14), `mechanistic` (Young eqs 9–13), `stampi_bombelli`
(Young eqs 15–16), or `none`.

**Trace regularisation** — q*_CO2 is multiplied by p/(p + p_reg), p_reg =
`isotherm.co2.p_regularization` = 0.01 Pa. The WADST wet site loads ~0.6 mol/kg
already at 1e-3 Pa (dq*/dp ~ 1e20 at p → 0): an artefact of the fit far outside its
data range that makes Newton ill-posed at the leading edge of the CO2 front. The
factor changes q* by < 0.03 % at the feed partial pressure (40 Pa).

## 8. Boundary conditions — `valves/valve_logic.py`, `boundaries/conditions.py`

Each valve has one of four states; the state of the open valve at each end
determines the end-face fluxes:

| state | physical meaning | face fluxes |
|---|---|---|
| CLOSED | sealed | N_k = 0, adiabatic |
| OPEN (V1) | feed blower, fixed molar flow | N_k = F y_feed,k (Danckwerts), F = u_feed P_feed/(R T_feed) × profile(r, θ) |
| OPEN_ATMOSPHERIC | vent / feed header | linear valve N = k_v (P_res − P)/(R T_up); out-flow carries bed gas, back-flow ambient (or feed) gas |
| OPEN_VACUUM | vacuum pump | N_out = k_p (P − P_vac)/(R T), check valve; out-flow carries bed gas |

* k_v = `operation.vent_conductance`, k_p = `operation.pump_speed_per_area`
  (= pump speed S / cross-section). The pump is what produces the blowdown pressure
  history (τ ≈ ε_t L / k_p) — it is derived, not imposed.
* Valves use **flux splitting**: N = N_in − N_out with N_in = k_v softplus(Δp)/(R T_ext)
  carrying reservoir gas and N_out = k_v softplus(−Δp)/(R T) carrying bed gas at bed
  temperature (smoothing 0.2 Pa). Each stream is upwinded exactly. A single blended stream
  (composition and temperature switched by a tanh of Δp) lets out-flowing gas carry less
  enthalpy than h(T_bed) near Δp = 0; a warmer cell then pushes out gas and heats up further
  — a spurious instability (Re λ > 0) at the end of pressurisation. The price of splitting is a
  small two-way exchange through an open valve whose pressure difference is below ~1 Pa
  (≈ 1 % of the feed flow; it only exists while the valve is idle).
* The vacuum pump with check valve uses N_out = k_p softplus(P − P_vac)/(R T): always ≥ 0,
  always bed composition (a velocity-based upwind switch would blend in the vacuum-line gas,
  because pump velocities k_p ΔP/P are only 10⁻⁵–10⁻³ m/s).
* Vents open over `operation.valve_time_constant`: the reservoir pressure seen by the
  bed ramps from the bed pressure to P_atm (a valve-opening time; it also prevents a
  pressure shock when a vent opens on a bed under vacuum).
* `operation.vent_model: ramp` / `operation.vacuum_model: ramp` give the classical
  Dirichlet-pressure ports instead (stiffer; kept for comparison).
* Feed distribution `feed.profile`: uniform, cosine_theta, linear_r, central_disc, sector
  (area-weighted mean 1, so the total feed is unchanged).
* Outer wall r = R: Robin (section 6). Axis: zero face area. θ: periodic.

**Deviation from the roadmap.** The roadmap maps OPEN_VACUUM to "C = 0 at the
outlet". A vacuum pump fixes the *total pressure*, not the composition: during
heating the gas leaving V3 is nearly pure CO2 + H2O, and imposing C_CO2 = 0 on the
face would create an artificial dispersive sink, break the CO2 balance and
over-predict desorption. Here outlets convect the bed composition
(zero-gradient/Danckwerts outlet) and only the pressure-driven flow is imposed.

## 9. Numerics

### Why implicit (and why RK4 is only a check) — `numerics/cfl.py`

The explicit stability limit is set by pressure diffusion:
Δt ≤ ½Δx²/D_P with D_P ≈ ε_b³d_p² P / (150(1−ε_b)² ε_t μ). For the default grid
this is ~1e-7 s, i.e. ~1e11 RK4 steps for one 8000 s adsorption step.
`core/ode_solver.py` (RK4) refuses runs that would need more than its step budget;
it exists to cross-check the implicit solver on tiny problems.

### Time integration — `core/solver.py`

Own implicit one-step integrator, `solver.scheme`:

* `sdirk2` (default): Alexander's L-stable, stiffly accurate 2-stage SDIRK,
  γ = 1 − 1/√2; both stages damp stiff modes.
* `trbdf2`: TR-BDF2 (Bank et al. 1985; Hosea & Shampine 1996), γ = 2 − √2.

Both stages share the iteration matrix M = I − γh J. Local error: difference to a
3rd-order quadrature of the stage derivatives, filtered by M⁻¹ (Shampine) so
that stiff components do not cause spurious rejections; h_new = h·clip(0.9 err^(−1/3), 0.2, 4).

Step-size control also respects Newton: when a stage fails to converge **with a
fresh Jacobian**, h is cut ×0.3 and an upper bound h ≤ 0.5 h_fail is remembered,
relaxed by ×1.15 per accepted step. Without this the accuracy controller grows h
×4 right back into the same failure (each failure costs ~200 RHS evaluations).

**Conservation polish.** A step-size convergence test alone is not enough for exact
balances: in stiff directions (pressure) a state error far below tolerance still leaves
a large residual R, and Σ_cells V·R is exactly the stage's mismatch between inventory
change and boundary flows. After convergence the solver therefore checks, for every
species (gas + adsorbed) and for energy, |W R| ≤ `solver.conservation_tol`·W(atol + rtol|y|)
(W = inventory weights) and, if needed, does up to two extra accurate iterations. This is
what makes the per-step mass and energy balances close to ~10⁻⁷.

Accumulators (boundary flows, heat duties, pump/blower work) are advanced with
the same stage weights as the state, so the mass and energy balances close to
the Newton tolerance.

### Nonlinear solver

`solver.newton_mode`:

* `krylov` (default): inexact Newton. Direction from GMRES on (I − γh J(z)) Δz = −R,
  with J·v evaluated as a finite-difference directional derivative (the *true*
  current Jacobian), **right**-preconditioned by the stored sparse LU (so GMRES
  minimises the true linear residual), with an absolute floor at 0.1κ so it does not
  chase round-off. Globalisation: backtracking (λ = 1 … 1/8) on the **affine-invariant**
  simplified Newton correction ‖M⁻¹R/s‖ (Deuflhard), not on ‖R/s‖ — CO2 at 400 ppm is
  carried by the stiff bulk flow, scaled Jacobian entries reach ~10⁷, and the raw residual
  is dominated by finite-difference noise that corresponds to state errors far below
  tolerance. Positivity truncation of concentrations/loadings (a component may drop by at
  most 90 % per iteration); early exit when the observed contraction rate cannot reach the
  tolerance within the remaining iterations (Hairer–Wanner test); if GMRES stalls with a
  Jacobian from an earlier step, the Jacobian is re-evaluated at the current iterate
  inside the stage solve and kept for the following steps.
* `hybrid`: a few modified-Newton iterations with the stored LU first.
* `modified`: classical modified Newton only.

Convergence is measured in the scaled norm ‖Δz/(atol + rtol|z|)‖_RMS ≤ κ
(`solver.newton_kappa`), with per-variable atol (`solver.atol_*`).

### Jacobian and preconditioner — `core/jacobian.py`, `numerics/sparsity.py`

* The exact sparsity pattern is built from the stencils (5-point MUSCL in z,
  3-point in r and θ, local isotherm/energy coupling, wall coupling); columns are
  coloured greedily so the full Jacobian costs one RHS evaluation per colour
  (~80 for the 3D model, independent of grid size). Verified against a dense FD
  Jacobian in `tests/`.
* `solver.preconditioner: theta_block` (default in 3D) factorises a reduced
  matrix that keeps couplings inside each θ-sector plus all couplings among the
  innermost `solver.precond_theta_rings` rings (whose small arcs r·Δθ make
  θ-coupling stiff). Exact for axisymmetric states, 5–10× cheaper LU in 3D; the
  dropped couplings are recovered by GMRES. `full` uses the whole Jacobian.
* SuperLU with COLAMD ordering; re-factorised when h changes by more than
  `solver.refactor_ratio`; Jacobian refreshed by age (`solver.jacobian_max_age`),
  after Newton failures, and when GMRES needs more than `solver.gmres_refresh`
  iterations in a step.

## 10. Performance indicators — `postprocessing/reporting.py`

Per cycle (the reported values are those of the last cycle; at CSS they are the
process KPIs):

| KPI | definition |
|---|---|
| recovery_pct | 100 · n_CO2,product / n_CO2,fed |
| purity_dry_pct | 100 · n_CO2,product / (n_CO2,product + n_N2,product) (water condensed, Young 2021) |
| purity_wet_pct | 100 · n_CO2,product / Σ_k n_k,product |
| productivity | kg CO2 / (m³ bed · day) and mol / (kg sorbent · day) |
| working_capacity | n_CO2,product / sorbent mass [mol/kg] |
| specific_thermal_kJ_per_kg | jacket (+ exchanger) heat in steps whose jacket is hotter than ambient / kg CO2 |
| specific_electric_kJ_per_kg | vacuum pump + feed blower, isentropic work / efficiency (`energy.*`) |
| work_equivalent_kJ_per_kg | W_el + η_turb (1 − T_L/T_H) Q_th (Young 2021) |

Product = everything leaving through the step's product port (V3 by default).

## 11. Balances and cyclic steady state — `validation/`

* Mass: for each step and species, (inventory_end − inventory_start) − (in − out)
  = residual. Energy: ΔE_bed+wall − (H_in − H_out + Q_jacket + Q_hx + Q_ambient)
  = residual. Both are written per step and per cycle (`mass_balance.csv`,
  `energy_balance.csv`).
* CSS: (1) volume-weighted RMS change of q, T and c between successive cycle
  ends < `run.css_tolerance`, and (2) cycle CO2 balance |fed − product − vented|/fed
  < `run.css_mb_tolerance` (Young 2021 used 0.5 %), after at least `run.min_cycles`.

## 12. Known limitations (read before trusting a number)

* No water condensation: with the jacket at 288 K and humid desorption gas the real
  bed/wall could condense; the model caps RH at 0.97 in the isotherm and carries
  supersaturated vapour instead.
* LDF constants are lumped fits to breakthrough curves at adsorption conditions;
  desorption kinetics at 100 °C are extrapolated (set `kinetics.Ea_*` > 0 to test).
* The isotherm parameters are for one batch of Lewatit VP OC 1065; Young 2021
  report sizeable batch-to-batch and lab-to-lab differences.
* Rim (jacket-only) heating is conduction-limited. During desorption under vacuum the bed
  is a distributed heat sink s (water and CO2 keep desorbing as the pump removes them);
  steady conduction from the rim to the axis gives T_wall − T_axis ≈ s R²/(4λ) ≈ 30 K for
  R = 5 cm with λ ≈ 0.15 W/m/K. Hence the default d = 5 cm (≈ 7 K) and the case
  `young2021_flat_bed_rim_jacket.yaml` (d = 10 cm) for comparison. This is physics, not a
  numerical artefact.
* The jacket coefficient h = 26 W/m²/K (Galanti 2025) is low for a liquid jacket; in the
  default case the wall lags the jacket by tens of K early in heating, so h_jacket is one of
  the most influential parameters (it is in the sensitivity list).
* Valve flux splitting allows a small two-way exchange (≈ 1 % of the feed flow) while an
  open valve sits within ~1 Pa of its reservoir (e.g. the hold at the end of a
  pressurisation step); it affects the KPIs at the 0.1 % level.
* Wall resolved only in (θ, z); the jacket is an ideal fluid at a set temperature.

## 13. Errata found in the source literature (and how this code handles them)

* Zhang et al. 2025 (Energy 334:137733), Table A.1 lists ΔH₀,dry = −70 000 J/mol and
  q_m = 3.36 mol/kg; the values consistent with Young 2021 (from which they are
  taken) are ΔH₀ = −117 798 J/mol (Toth affinity) and q_m = 3.63 mol/kg (GAB). The
  −70 kJ/mol is Young's *isosteric heat* used in the energy balance, a different
  quantity. This code uses Young's values.
* Zhang 2025 eq. A.1 evaluates the wet-site affinity with b_dry; it must be b_wet.
* The project's equation-set notes copy the b_dry typo, and write the work-equivalent
  of heat as Q/η_turb; it is η_turb·(1 − T_L/T_H)·Q (a heat engine converts heat to
  work with efficiency < 1, it does not amplify it).
