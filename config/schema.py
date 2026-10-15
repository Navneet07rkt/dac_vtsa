"""Validation schema + loader for dac_vtsa configuration.

Every tunable number lives in ``defaults.yaml``.  This module

* declares, for every key, its type, unit, allowed range, meaning and literature
  source (``SPECS``) -- the single place that documents parameters;
* loads ``defaults.yaml``, deep-merges an optional user YAML and ``--set``
  overrides, resolves string references such as ``operation.T_heat``;
* validates everything (unknown keys are errors -> typos are caught before a
  long run silently produces garbage);
* returns an immutable-ish :class:`Config` with attribute access and a few
  derived quantities (bulk density, particle porosity, ...).

``python -m config.schema --markdown`` regenerates PARAMETERS.md.
"""
from __future__ import annotations

import copy
import math
import os
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULTS_PATH = os.path.join(HERE, "defaults.yaml")

Y21 = "Young et al. 2021, Energy Environ. Sci. 14:5377 (doi:10.1039/D1EE01272J), ESI Table S9"
Y21T2 = "Young et al. 2021 EES, Table 2 (Toth fit, Lewatit VP OC 1065)"
Y21T3 = "Young et al. 2021 EES, Table 3 (GAB fit)"
Y21T4 = "Young et al. 2021 EES, Table 4 (co-adsorption fits)"
G25 = "Galanti et al. 2025, RSC Energy Adv. (doi:10.1039/D5YA00336A), SI Tables S7-S8"
ER68 = "Edwards & Richardson 1968, Chem. Eng. Sci. 23:109"
ZS70 = "Zehner & Schlünder 1970, Chem.-Ing.-Tech. 42:933 (VDI Heat Atlas D6.3)"
MN93 = "Martin & Nilles 1993, Chem.-Ing.-Tech. 65:1468 (as used in G25 SI)"
ASSUME = "modelling assumption / user choice"
Z25 = "Zhang et al. 2025, Energy 334:137733, Table 1"


@dataclass(frozen=True)
class Spec:
    key: str
    type: str                      # float | int | str | bool | list | any
    unit: str = "-"
    lo: Optional[float] = None
    hi: Optional[float] = None
    choices: Optional[Tuple[Any, ...]] = None
    desc: str = ""
    source: str = ASSUME
    nullable: bool = False


def _s(key, typ, unit="-", lo=None, hi=None, desc="", source=ASSUME, choices=None, nullable=False):
    return Spec(key, typ, unit, lo, hi, tuple(choices) if choices else None, desc, source, nullable)


SPECS: List[Spec] = [
    # ---------------------------------------------------------------- run
    _s("run.name", "str", desc="run name; outputs go to <results_dir>/<name>/"),
    _s("run.results_dir", "str", desc="output root directory (relative to project root)"),
    _s("run.max_cycles", "int", "-", 1, 500, "maximum number of TVSA cycles"),
    _s("run.min_cycles", "int", "-", 1, 500, "minimum number of cycles before CSS may be declared"),
    _s("run.css_tolerance", "float", "-", 1e-8, 1.0,
       "cyclic-steady-state tolerance on normalised end-of-cycle state change"),
    _s("run.css_mb_tolerance", "float", "-", 1e-8, 1.0,
       "CSS tolerance on the cycle CO2 balance |in-out|/in", "Young 2021 used 0.5 %"),
    _s("run.handoff_via_csv", "bool", desc="each step reloads the state CSV written by the previous step"),
    _s("run.snapshots_per_step", "int", "-", 2, 10000, "number of field snapshots stored per step"),
    _s("run.write_fields_every_cycle", "bool", desc="write field CSVs for every cycle (else only the last)"),
    _s("run.initial_state", "str", choices=("clean", "feed_equilibrium", "csv"),
       desc="initial bed: clean sorbent in feed air, saturated with feed, or loaded from CSV"),
    _s("run.initial_state_csv", "str", nullable=True, desc="state CSV used when initial_state == csv"),
    _s("run.make_plots", "bool", desc="write PNG plots at the end of the run"),
    # ---------------------------------------------------------------- components
    _s("components.water", "bool", desc="include H2O (humid air, co-adsorption). false = dry CO2/N2 model"),
    # ---------------------------------------------------------------- geometry
    _s("geometry.bed_length", "float", "m", 1e-4, 50.0, "packed length L (flow direction z)", Y21),
    _s("geometry.bed_radius", "float", "m", 1e-4, 10.0, "inner radius R of the bed",
       "design choice: d = 5 cm; Young 2021 ESI S9 flat bed has d = 10 cm, which is conduction-limited "
       "with a rim jacket (cases/young2021_flat_bed_rim_jacket.yaml)"),
    _s("geometry.wall_thickness", "float", "m", 1e-5, 0.5, "metal wall thickness", ASSUME),
    # ---------------------------------------------------------------- grid
    _s("grid.Nr", "int", "-", 1, 200, "radial finite volumes"),
    _s("grid.Ntheta", "int", "-", 1, 360, "angular finite volumes (1 = axisymmetric 2D)"),
    _s("grid.Nz", "int", "-", 2, 2000, "axial finite volumes"),
    _s("grid.radial_spacing", "str", choices=("uniform", "equal_area"), desc="radial face placement"),
    # ---------------------------------------------------------------- sorbent
    _s("sorbent.name", "str", desc="label only"),
    _s("sorbent.particle_diameter", "float", "m", 1e-6, 0.05, "bead diameter d_p", Y21 + " (R_p = 0.26 mm)"),
    _s("sorbent.particle_density", "float", "kg/m3", 50, 5000, "pellet density rho_p", Y21),
    _s("sorbent.bed_voidage", "float", "-", 0.2, 0.95, "inter-particle voidage eps_b", Y21),
    _s("sorbent.total_voidage", "float", "-", 0.2, 0.99, "total voidage eps_t = eps_b + eps_p(1-eps_b)", Y21),
    _s("sorbent.cp_solid", "float", "J/kg/K", 100, 5000, "sorbent heat capacity", Y21),
    _s("sorbent.particle_conductivity", "float", "W/m/K", 0.01, 50, "solid (bead) thermal conductivity",
       Z25 + " (Galanti 2025 SI uses 0.43)"),
    _s("sorbent.max_temperature", "float", "K", 300, 700, "temperature above which a warning is logged",
       "Young 2021 (max 100 C regeneration for oxidative stability)"),
    # ---------------------------------------------------------------- isotherm CO2
    _s("isotherm.co2.q_inf0", "float", "mol/kg", 0.01, 50, "Toth saturation capacity at T0", Y21T2),
    _s("isotherm.co2.chi", "float", "-", -50, 50, "Toth capacity T-dependence", Y21T2),
    _s("isotherm.co2.b0", "float", "1/Pa", 0, 1e10, "Toth affinity pre-factor (b = b0 exp(-dH0/RT))", Y21T2),
    _s("isotherm.co2.dH0", "float", "J/mol", -5e5, 0, "Toth affinity energy (negative)", Y21T2),
    _s("isotherm.co2.tau0", "float", "-", 0.01, 5, "Toth heterogeneity at T0", Y21T2),
    _s("isotherm.co2.alpha", "float", "-", -10, 10, "Toth heterogeneity T-dependence", Y21T2),
    _s("isotherm.co2.T0", "float", "K", 200, 500, "Toth reference temperature", Y21T2),
    _s("isotherm.co2.p_regularization", "float", "Pa", 0.0, 10.0,
       "trace-level regularisation q* -> q* p/(p+p_reg); 0 disables (numerically harsh with WADST)"),
    # ---------------------------------------------------------------- isotherm H2O
    _s("isotherm.h2o.q_m", "float", "mol/kg", 0.01, 100, "GAB monolayer capacity", Y21T3),
    _s("isotherm.h2o.C", "float", "J/mol", 0, 2e5, "GAB: E1 = C - exp(D T)", Y21T3),
    _s("isotherm.h2o.D", "float", "1/K", 0, 1, "GAB: E1 = C - exp(D T)", Y21T3),
    _s("isotherm.h2o.F", "float", "J/mol", 0, 2e5, "GAB: E2-9 = F + G T", Y21T3),
    _s("isotherm.h2o.G", "float", "J/mol/K", -1000, 1000, "GAB: E2-9 = F + G T", Y21T3),
    _s("isotherm.h2o.E10_slope", "float", "J/mol/K", -1000, 1000, "E10+ = slope*T + intercept", "Young 2021 eq. 8"),
    _s("isotherm.h2o.E10_intercept", "float", "J/mol", 0, 2e5, "E10+ = slope*T + intercept", "Young 2021 eq. 8"),
    _s("isotherm.h2o.rh_cap", "float", "-", 0.5, 0.999,
       "numerical cap on relative humidity passed to GAB (no condensation model)"),
    # ---------------------------------------------------------------- co-adsorption
    _s("isotherm.coadsorption.model", "str", choices=("none", "wadst", "mechanistic", "stampi_bombelli"),
       desc="effect of adsorbed water on CO2 isotherm", source="Young 2021 EES"),
    _s("isotherm.coadsorption.water_loading_basis", "str", choices=("state", "equilibrium"),
       desc="use local adsorbed water q_H2O (state) or q*_H2O in the CO2 isotherm"),
    _s("isotherm.coadsorption.wadst.b0_wet", "float", "1/Pa", 0, 1e10, "WADST wet-site affinity pre-factor", Y21T4),
    _s("isotherm.coadsorption.wadst.q_inf0_wet", "float", "mol/kg", 0, 100, "WADST wet-site capacity", Y21T4),
    _s("isotherm.coadsorption.wadst.tau0_wet", "float", "-", 0.001, 5, "WADST wet-site Toth tau0", Y21T4),
    _s("isotherm.coadsorption.wadst.chi_wet", "float", "-", -50, 50, "WADST wet-site chi", Y21T4),
    _s("isotherm.coadsorption.wadst.alpha_wet", "float", "-", -10, 10, "WADST wet-site alpha", Y21T4),
    _s("isotherm.coadsorption.wadst.dH_wet", "float", "J/mol", -5e5, 0, "WADST wet-site affinity energy", Y21T4),
    _s("isotherm.coadsorption.wadst.A", "float", "mol/kg", 1e-6, 100, "WADST critical water loading", Y21T4),
    _s("isotherm.coadsorption.mechanistic.f_blocked_max", "float", "-", 0, 1, "max blocked amine fraction", Y21T4),
    _s("isotherm.coadsorption.mechanistic.k", "float", "kg/mol", 0, 100, "Avrami blocking constant", Y21T4),
    _s("isotherm.coadsorption.mechanistic.phi_dry", "float", "-", 0.01, 2, "dry amine efficiency", Y21T4),
    _s("isotherm.coadsorption.mechanistic.A", "float", "mol/kg", 1e-6, 100, "critical water loading", Y21T4),
    _s("isotherm.coadsorption.mechanistic.dH_wet", "float", "J/mol", -5e5, 0, "wet affinity energy", Y21T4),
    _s("isotherm.coadsorption.mechanistic.n", "float", "-", 0.1, 10, "Avrami exponent", Y21T4),
    _s("isotherm.coadsorption.stampi_bombelli.gamma", "float", "kg/mol", -10, 10, "capacity factor", Y21T4),
    _s("isotherm.coadsorption.stampi_bombelli.beta", "float", "kg/mol", -10, 100, "affinity factor", Y21T4),
    # ---------------------------------------------------------------- kinetics
    _s("kinetics.k_ldf_co2", "float", "1/s", 1e-7, 100, "LDF coefficient CO2 at T_ref", Y21 + " (breakthrough fit)"),
    _s("kinetics.k_ldf_h2o", "float", "1/s", 1e-7, 100, "LDF coefficient H2O at T_ref", Y21 + " (breakthrough fit)"),
    _s("kinetics.Ea_co2", "float", "J/mol", 0, 2e5, "Arrhenius activation energy of k_LDF,CO2 (0 = constant)",
       "Young 2021 / Galanti 2025 use constant k"),
    _s("kinetics.Ea_h2o", "float", "J/mol", 0, 2e5, "Arrhenius activation energy of k_LDF,H2O"),
    _s("kinetics.T_ref", "float", "K", 200, 500, "reference temperature of k_LDF"),
    # ---------------------------------------------------------------- thermo
    _s("thermo.dH_ads_co2", "float", "J/mol", -3e5, 0, "heat of adsorption CO2 (energy balance)", Y21),
    _s("thermo.dH_ads_h2o", "float", "J/mol", -3e5, 0, "heat of adsorption H2O (energy balance)", Y21),
    _s("thermo.cp_ads_co2", "float", "J/mol/K", 0, 500, "adsorbed-phase heat capacity CO2", ASSUME + " (gas-like)"),
    _s("thermo.cp_ads_h2o", "float", "J/mol/K", 0, 500, "adsorbed-phase heat capacity H2O",
       ASSUME + " (liquid water 75.3; gives Kirchhoff T-dependence of dH)"),
    _s("thermo.T_ref", "float", "K", 200, 500, "reference temperature of enthalpies / dH_ads"),
    # ---------------------------------------------------------------- properties
    _s("properties.backend", "str", choices=("ideal_gas", "coolprop"),
       desc="gas property backend (coolprop needs `pip install CoolProp`)"),
    _s("properties.cp_reference_T", "float", "K", 200, 600,
       "temperature at which constant gas cp are evaluated for the energy balance"),
    # ---------------------------------------------------------------- transport
    _s("transport.dispersion_model", "str", choices=("edwards_richardson", "young2021", "wakao_kagei"),
       desc="axial dispersion correlation", source=ER68),
    _s("transport.er_gamma1", "float", "-", 0, 5, "E-R tortuosity factor gamma1", ER68),
    _s("transport.er_gamma2", "float", "-", 0, 5, "E-R gamma2", ER68),
    _s("transport.er_beta", "float", "-", 0, 100, "E-R beta (beta*gamma1 = 9.49)", ER68),
    _s("transport.Pe_mass_radial_inf", "float", "-", 1, 100, "asymptotic transverse mass Peclet number",
       "Delgado 2006 review (8-12)"),
    _s("transport.Pe_heat_radial_inf", "float", "-", 1, 100, "asymptotic radial thermal Peclet number", G25),
    _s("transport.Pe_heat_axial_inf", "float", "-", 0.1, 100, "asymptotic axial thermal Peclet number",
       "Yagi & Kunii (0.5 Re Pr)"),
    _s("transport.bed_conductivity_model", "str", choices=("zehner_schlunder", "parallel", "constant"),
       desc="stagnant effective bed conductivity", source=ZS70),
    _s("transport.bed_conductivity_constant", "float", "W/m/K", 1e-3, 50, "used if model == constant"),
    _s("transport.zs_shape_factor", "float", "-", 0.5, 5, "Zehner-Schluender shape factor (1.25 spheres)", ZS70),
    _s("transport.wall_htc_model", "str", choices=("martin_nilles", "constant"),
       desc="bed-to-wall heat transfer coefficient", source=MN93),
    _s("transport.wall_htc_constant", "float", "W/m2/K", 0.1, 1e5, "used if wall_htc_model == constant"),
    _s("transport.ergun_A", "float", "-", 0, 1e4, "Ergun viscous constant"),
    _s("transport.ergun_B", "float", "-", 0, 100, "Ergun inertial constant"),
    _s("transport.heat_transfer_mode", "str", choices=("resolved", "lumped"),
       desc="resolved: conduction resistances in the bed; lumped: 1D-style overall U at the wall"),
    _s("transport.lumped_U", "float", "W/m2/K", 0.01, 1e5, "overall coefficient for lumped mode",
       "Young 2021 ESI Table S9 (h_j = 14)"),
    # ---------------------------------------------------------------- wall / jacket
    _s("wall.enabled", "bool", desc="resolve wall heat capacity and conduction"),
    _s("wall.density", "float", "kg/m3", 100, 30000, "wall density", G25),
    _s("wall.cp", "float", "J/kg/K", 100, 5000, "wall heat capacity", G25),
    _s("wall.conductivity", "float", "W/m/K", 0.01, 500, "wall conductivity", G25),
    _s("jacket.h_jacket", "float", "W/m2/K", 0, 1e5, "jacket fluid to wall heat transfer coefficient",
       G25 + " (h_ext = 26)"),
    _s("jacket.theta_coverage_deg", "list", "deg", desc="[start, end] angular extent of the jacket (0-360)"),
    _s("jacket.z_coverage", "list", "-", desc="[start, end] axial extent of the jacket as fractions of L"),
    _s("jacket.h_ambient", "float", "W/m2/K", 0, 1e4, "uncovered wall to ambient coefficient"),
    _s("jacket.ramp_time", "float", "s", 0, 1e5, "first-order lag of jacket fluid temperature (0 = step)"),
    _s("internal_hx.enabled", "bool", desc="embedded volumetric heat exchanger"),
    _s("internal_hx.hA_per_volume", "float", "W/m3/K", 0, 1e7, "volumetric hA of embedded exchanger"),
    _s("internal_hx.follows_jacket", "bool", desc="exchanger fluid temperature = jacket temperature"),
    # ---------------------------------------------------------------- feed / ambient
    _s("feed.T", "float", "K", 200, 400, "feed temperature", Y21),
    _s("feed.P", "float", "Pa", 1e4, 1e6, "feed pressure (also outlet pressure in adsorption)", Y21),
    _s("feed.RH", "float", "-", 0, 1, "feed relative humidity", Y21),
    _s("feed.y_co2_dry", "float", "-", 0, 1, "CO2 mole fraction in dry air", Y21),
    _s("feed.velocity", "float", "m/s", 1e-6, 10, "superficial feed velocity at feed T,P", Y21),
    _s("feed.profile", "str", choices=("uniform", "cosine_theta", "linear_r", "central_disc", "sector"),
       desc="inlet flow distribution (maldistribution / partial-inlet studies)"),
    _s("feed.profile_amplitude", "float", "-", 0, 0.99, "amplitude of the inlet maldistribution"),
    _s("feed.profile_theta0_deg", "float", "deg", 0, 360, "angle of maximum flow (cosine_theta, sector)"),
    _s("feed.profile_radius_fraction", "float", "-", 0.05, 1.0, "inlet radius / R for central_disc"),
    _s("feed.profile_sector_deg", "float", "deg", 1.0, 360.0, "angular width of the inlet for sector"),
    _s("ambient.T", "float", "K", 200, 400, "ambient temperature (vent back-flow, heat loss)"),
    _s("ambient.P", "float", "Pa", 1e4, 1e6, "ambient pressure (atmospheric vent, pump discharge)"),
    _s("ambient.RH", "float", "-", 0, 1, "ambient relative humidity (vent back-flow)"),
    _s("ambient.y_co2_dry", "float", "-", 0, 1, "ambient CO2 dry mole fraction"),
    # ---------------------------------------------------------------- operation
    _s("operation.P_vacuum", "float", "Pa", 100, 1e5, "vacuum level (blowdown / desorption)", Y21),
    _s("operation.T_heat", "float", "K", 273, 600, "jacket temperature in heating", Y21),
    _s("operation.T_cool", "float", "K", 200, 600, "jacket temperature in cooling", Y21),
    _s("operation.T_adsorption_jacket", "float", "K", 200, 600, "jacket temperature in adsorption", Y21),
    _s("operation.valve_time_constant", "float", "s", 1e-3, 1e4,
       "time constant of the exponential pressure ramp at an opened valve"),
    _s("operation.vent_model", "str", choices=("valve", "ramp"),
       desc="atmospheric vent / feed header: finite-conductance valve or imposed exponential P ramp"),
    _s("operation.vent_conductance", "float", "m/s", 1e-6, 1e4,
       "vent valve conductance per unit cross-section, N = k_v (P_res - P)/(R T)"),
    _s("operation.vacuum_model", "str", choices=("pump", "ramp"),
       desc="vacuum port: finite-speed pump (P relaxes to P_vacuum) or imposed exponential P ramp"),
    _s("operation.pump_speed_per_area", "float", "m/s", 1e-7, 100.0,
       "vacuum pump volumetric speed per unit bed cross-section (S/A); sets the blowdown time constant "
       "eps_t L/(S/A)"),
    _s("operation.vacuum_check_valve", "bool", desc="prevent back-flow from vacuum lines"),
    _s("ports.V1", "str", choices=("top", "bottom"), desc="end of the column where V1 (inlet) sits"),
    _s("ports.V2", "str", choices=("top", "bottom"), desc="end of V2 (vent / blowdown)"),
    _s("ports.V3", "str", choices=("top", "bottom"), desc="end of V3 (product)", source="Young 2021 (z = L)"),
    _s("cycle", "list", desc="ordered list of steps (see README 'Cycle definition')"),
    # ---------------------------------------------------------------- energy accounting
    _s("energy.eta_vacuum", "float", "-", 0.01, 1, "vacuum pump isentropic efficiency", Y21),
    _s("energy.vacuum_efficiency_model", "str", choices=("constant", "krishnamurthy"),
       desc="constant, or 0.8*19.55p/(1+19.55p) (p in bar) as in G25 / Krishnamurthy et al. 2014"),
    _s("energy.eta_blower", "float", "-", 0.01, 1, "feed blower efficiency", Y21),
    _s("energy.gamma_gas", "float", "-", 1.0, 1.7, "cp/cv for compression work", Y21),
    _s("energy.pump_handles_water", "bool",
       desc="include H2O in the vacuum-pump molar flow (false: water condensed before compression)",
       source="assumption; Young 2021 condense water 'in the compression process' without giving the work convention"),
    _s("energy.eta_turbine", "float", "-", 0, 1, "efficiency for heat -> work equivalent", Y21),
    _s("energy.T_low_equiv", "float", "K", 200, 400, "T_L in work equivalent", Y21),
    # ---------------------------------------------------------------- solver
    _s("solver.method", "str", choices=("trbdf2", "scipy_bdf", "rk4"),
       desc="time integrator: implicit Newton TR-BDF2 (default), scipy BDF cross-check, explicit RK4"),
    _s("solver.scheme", "str", choices=("sdirk2", "trbdf2"),
       desc="2-stage implicit RK used by the Newton integrator: L-stable SDIRK2 or TR-BDF2"),
    _s("solver.rtol", "float", "-", 1e-10, 1e-1, "relative tolerance"),
    _s("solver.atol_c_co2", "float", "mol/m3", 1e-14, 1, "absolute tolerance CO2 concentration"),
    _s("solver.atol_c_h2o", "float", "mol/m3", 1e-14, 1, "absolute tolerance H2O concentration"),
    _s("solver.atol_c_inert", "float", "mol/m3", 1e-14, 10, "absolute tolerance inert concentration"),
    _s("solver.atol_q", "float", "mol/kg", 1e-14, 1, "absolute tolerance loadings"),
    _s("solver.atol_T", "float", "K", 1e-8, 10, "absolute tolerance temperature"),
    _s("solver.h_initial", "float", "s", 1e-10, 100, "initial time step"),
    _s("solver.h_max_fraction", "float", "-", 1e-5, 1, "max time step as a fraction of step duration"),
    _s("solver.max_newton_iter", "int", "-", 2, 50, "Newton iterations per stage"),
    _s("solver.preconditioner", "str", choices=("theta_block", "full"),
       desc="LU used by Newton: full Jacobian, or theta-sector-decoupled (exact for axisymmetric states)"),
    _s("solver.precond_theta_rings", "int", "-", 0, 1000,
       "number of innermost rings whose theta couplings are kept in the theta_block preconditioner"),
    _s("solver.refactor_ratio", "float", "-", 1.05, 100.0, "h change factor that triggers re-factorisation"),
    _s("solver.gmres_refresh", "int", "-", 1, 100000, "GMRES iterations per step that trigger a new Jacobian"),
    _s("solver.newton_mode", "str", choices=("hybrid", "krylov", "modified"),
       desc="stage solver: modified Newton with stored LU, Newton-Krylov (GMRES, LU-preconditioned), or both"),
    _s("solver.newton_kappa", "float", "-", 1e-4, 1, "Newton convergence safety factor"),
    _s("solver.conservation_tol", "float", "-", 0.0, 1.0,
       "stage balance tolerance: |W R| <= tol * W (atol + rtol|y|) per species and energy (0 = off)"),
    _s("solver.jacobian_max_age", "int", "-", 1, 10000, "accepted steps before forced Jacobian refresh"),
    _s("solver.limiter", "str", choices=("van_albada", "van_albada_tvd", "van_leer", "minmod", "upwind"),
       desc="axial TVD flux limiter"),
    _s("solver.max_steps_per_step", "int", "-", 10, 10 ** 8, "safety cap on time steps per cycle step"),
    _s("solver.outlet_record_every", "int", "-", 1, 10 ** 6, "record outlet time series every n steps"),
    # ---------------------------------------------------------------- sensitivity
    _s("sensitivity.relative_step", "float", "-", 1e-4, 0.9, "relative perturbation for sensitivities"),
    _s("sensitivity.cycles", "int", "-", 1, 200, "cycles per perturbed run"),
    _s("sensitivity.grid", "any", desc="grid override for sensitivity runs {Nr, Ntheta, Nz}"),
    _s("sensitivity.processes", "int", "-", 1, 256, "parallel worker processes"),
    _s("sensitivity.parameters", "list", desc="dotted parameter keys to perturb"),
    _s("sensitivity.kpis", "list", desc="KPI names to report"),
]

SPEC_BY_KEY: Dict[str, Spec] = {s.key: s for s in SPECS}

STEP_TYPES = ("adsorption", "blowdown", "heating", "cooling", "pressurization")
VALVES = ("V1", "V2", "V3")
VALVE_STATES = ("CLOSED", "OPEN", "OPEN_ATMOSPHERIC", "OPEN_VACUUM")


class ConfigError(ValueError):
    pass


class Config:
    """Nested, attribute-accessible view of a validated configuration dict."""

    def __init__(self, data: Dict[str, Any], _root: Optional["Config"] = None):
        object.__setattr__(self, "_data", data)
        object.__setattr__(self, "_root", _root if _root is not None else self)

    def __getattr__(self, name):
        try:
            v = self._data[name]
        except KeyError as exc:
            raise AttributeError(name) from exc
        if isinstance(v, dict):
            return Config(v, self._root)
        return v

    def __setattr__(self, key, value):
        raise AttributeError("Config is read-only; use parameters.load_config(overrides=...)")

    def __getitem__(self, name):
        return self.__getattr__(name)

    def get(self, dotted: str, default=None):
        node: Any = self._data
        for part in dotted.split("."):
            if isinstance(node, dict) and part in node:
                node = node[part]
            else:
                return default
        return node

    def to_dict(self) -> Dict[str, Any]:
        return copy.deepcopy(self._data)

    def __repr__(self):
        return f"Config({list(self._data.keys())})"


# ----------------------------------------------------------------------------- helpers

def deep_merge(base: Dict, upd: Dict) -> Dict:
    out = copy.deepcopy(base)
    for k, v in (upd or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def set_dotted(d: Dict, dotted: str, value: Any) -> None:
    parts = dotted.split(".")
    node = d
    for p in parts[:-1]:
        if p not in node or not isinstance(node[p], dict):
            node[p] = {}
        node = node[p]
    node[parts[-1]] = value


def get_dotted(d: Dict, dotted: str):
    node = d
    for p in dotted.split("."):
        node = node[p]
    return node


def _flatten(d: Dict, prefix: str = "") -> Dict[str, Any]:
    out = {}
    for k, v in d.items():
        key = f"{prefix}.{k}" if prefix else str(k)
        if isinstance(v, dict) and key not in SPEC_BY_KEY:
            out.update(_flatten(v, key))
        else:
            out[key] = v
    return out


def parse_override(text: str) -> Tuple[str, Any]:
    """'a.b.c=1.5' -> ('a.b.c', 1.5) using YAML scalar parsing."""
    if "=" not in text:
        raise ConfigError(f"override '{text}' must look like key=value")
    k, v = text.split("=", 1)
    return k.strip(), yaml.safe_load(v)


def _check_value(spec: Spec, v: Any) -> Any:
    if v is None:
        if spec.nullable:
            return None
        raise ConfigError(f"{spec.key}: value is required")
    t = spec.type
    if t in ("float", "int") and isinstance(v, str):
        # PyYAML (YAML 1.1) reads e.g. 2.0e4 as a string -> accept numeric strings
        try:
            v = float(v)
        except ValueError:
            raise ConfigError(f"{spec.key}: expected a number, got {v!r}")
    if t == "float":
        if isinstance(v, bool) or not isinstance(v, (int, float)):
            raise ConfigError(f"{spec.key}: expected a number, got {v!r}")
        v = float(v)
        if not math.isfinite(v):
            raise ConfigError(f"{spec.key}: must be finite")
    elif t == "int":
        if isinstance(v, bool) or not isinstance(v, (int, float)) or int(v) != v:
            raise ConfigError(f"{spec.key}: expected an integer, got {v!r}")
        v = int(v)
    elif t == "bool":
        if not isinstance(v, bool):
            raise ConfigError(f"{spec.key}: expected true/false, got {v!r}")
    elif t == "str":
        if not isinstance(v, str):
            raise ConfigError(f"{spec.key}: expected a string, got {v!r}")
    elif t == "list":
        if not isinstance(v, (list, tuple)):
            raise ConfigError(f"{spec.key}: expected a list, got {v!r}")
        v = list(v)
    if spec.choices is not None and v not in spec.choices:
        raise ConfigError(f"{spec.key}: {v!r} not in {spec.choices}")
    if t in ("float", "int"):
        if spec.lo is not None and v < spec.lo:
            raise ConfigError(f"{spec.key} = {v} below allowed minimum {spec.lo} [{spec.unit}]")
        if spec.hi is not None and v > spec.hi:
            raise ConfigError(f"{spec.key} = {v} above allowed maximum {spec.hi} [{spec.unit}]")
    return v


def _is_number(x: Any) -> bool:
    try:
        float(x)
        return True
    except (TypeError, ValueError):
        return False


def _resolve_ref(raw: Dict, value: Any) -> Any:
    """Allow step entries like  T: operation.T_heat  (string reference)."""
    if isinstance(value, str) and "." in value:
        try:
            return float(get_dotted(raw, value))
        except (KeyError, TypeError):
            raise ConfigError(f"cannot resolve reference '{value}'")
    return value


def _validate_cycle(raw: Dict) -> List[Dict[str, Any]]:
    steps = raw.get("cycle")
    if not isinstance(steps, list) or not steps:
        raise ConfigError("cycle: must be a non-empty list of steps")
    out = []
    names = set()
    for i, st in enumerate(steps):
        if not isinstance(st, dict):
            raise ConfigError(f"cycle[{i}] must be a mapping")
        st = copy.deepcopy(st)
        allowed = {"name", "type", "duration", "valves", "jacket", "stop_on", "product_port", "notes"}
        extra = set(st) - allowed
        if extra:
            raise ConfigError(f"cycle[{i}]: unknown keys {sorted(extra)}")
        name = st.get("name") or st.get("type")
        if name in names:
            raise ConfigError(f"cycle: duplicate step name '{name}'")
        names.add(name)
        st["name"] = name
        if st.get("type") not in STEP_TYPES:
            raise ConfigError(f"cycle[{i}].type must be one of {STEP_TYPES}")
        dur = _resolve_ref(raw, st.get("duration"))
        if not isinstance(dur, (int, float)) or dur <= 0:
            raise ConfigError(f"cycle[{i}].duration must be > 0 s")
        st["duration"] = float(dur)
        valves = st.get("valves") or {}
        vv = {}
        for v in VALVES:
            s = valves.get(v, "CLOSED")
            if s not in VALVE_STATES:
                raise ConfigError(f"cycle[{i}].valves.{v} = {s!r} not in {VALVE_STATES}")
            vv[v] = s
        extra_v = set(valves) - set(VALVES)
        if extra_v:
            raise ConfigError(f"cycle[{i}].valves: unknown valves {sorted(extra_v)}")
        st["valves"] = vv
        jk = st.get("jacket") or {"mode": "adiabatic"}
        mode = jk.get("mode", "temperature")
        if mode not in ("temperature", "adiabatic"):
            raise ConfigError(f"cycle[{i}].jacket.mode must be temperature|adiabatic")
        jout = {"mode": mode}
        if mode == "temperature":
            # a reference such as 'operation.T_heat' is kept next to its value (T_source), so that a
            # re-validated config (config_used.yaml, sensitivity runs with operation.T_heat
            # perturbed) re-resolves it instead of freezing the old number
            src = jk.get("T_source", jk.get("T"))
            T = _resolve_ref(raw, src)
            if not isinstance(T, (int, float)) or not (150 < T < 900):
                raise ConfigError(f"cycle[{i}].jacket.T invalid ({T!r})")
            jout["T"] = float(T)
            if isinstance(src, str) and "." in src and not _is_number(src):
                jout["T_source"] = src
        st["jacket"] = jout
        stop = st.get("stop_on")
        if stop is not None:
            crit = stop.get("criterion")
            if crit not in ("outlet_ratio", "end_saturation", "min_temperature", "max_temperature"):
                raise ConfigError(f"cycle[{i}].stop_on.criterion invalid: {crit!r}")
            vsrc = stop.get("value_source", stop.get("value"))
            stop_new = {"criterion": crit, "value": float(_resolve_ref(raw, vsrc))}
            if isinstance(vsrc, str) and "." in vsrc and not _is_number(vsrc):
                stop_new["value_source"] = vsrc
            stop = stop_new
        st["stop_on"] = stop
        pp = st.get("product_port")
        if pp is None:
            pp = "V3"
        if pp not in VALVES + ("none",):
            raise ConfigError(f"cycle[{i}].product_port must be V1/V2/V3/none")
        st["product_port"] = pp
        out.append(st)
    return out


def _validate_ports_vs_cycle(raw: Dict, steps: List[Dict]) -> None:
    ports = raw["ports"]
    for st in steps:
        open_by_end: Dict[str, List[str]] = {"top": [], "bottom": []}
        for v, s in st["valves"].items():
            if s != "CLOSED":
                open_by_end[ports[v]].append(v)
        for end, lst in open_by_end.items():
            if len(lst) > 1:
                raise ConfigError(f"step '{st['name']}': valves {lst} are open at the same end ({end})")
        for v, s in st["valves"].items():
            if s == "OPEN" and v != "V1":
                raise ConfigError(f"step '{st['name']}': state OPEN (feed blower) is only valid for V1")


def validate(raw: Dict[str, Any]) -> Config:
    raw = copy.deepcopy(raw)
    flat = _flatten({k: v for k, v in raw.items() if k != "cycle"})
    unknown = [k for k in flat if k not in SPEC_BY_KEY]
    if unknown:
        raise ConfigError(f"unknown configuration keys (typo?): {unknown}")
    for spec in SPECS:
        if spec.key == "cycle":
            continue
        try:
            v = get_dotted(raw, spec.key)
        except (KeyError, TypeError):
            raise ConfigError(f"missing configuration key: {spec.key}")
        set_dotted(raw, spec.key, _check_value(spec, v))
    # ---------------- cross-field checks
    s = raw["sorbent"]
    if s["total_voidage"] < s["bed_voidage"]:
        raise ConfigError("sorbent.total_voidage must be >= sorbent.bed_voidage")
    if raw["run"]["min_cycles"] > raw["run"]["max_cycles"]:
        raise ConfigError("run.min_cycles > run.max_cycles")
    th = raw["jacket"]["theta_coverage_deg"]
    zc = raw["jacket"]["z_coverage"]
    if len(th) != 2 or not (0 <= th[0] < th[1] <= 360):
        raise ConfigError("jacket.theta_coverage_deg must be [start, end] with 0 <= start < end <= 360")
    if len(zc) != 2 or not (0 <= zc[0] < zc[1] <= 1):
        raise ConfigError("jacket.z_coverage must be [start, end] fractions with 0 <= start < end <= 1")
    if raw["operation"]["P_vacuum"] >= raw["ambient"]["P"]:
        raise ConfigError("operation.P_vacuum must be below ambient.P")
    if raw["grid"]["Ntheta"] == 2:
        raise ConfigError("grid.Ntheta = 2 is degenerate; use 1 (axisymmetric) or >= 3")
    steps = _validate_cycle(raw)
    _validate_ports_vs_cycle(raw, steps)
    raw["cycle"] = steps
    # ---------------- derived quantities
    eb, et = s["bed_voidage"], s["total_voidage"]
    raw["derived"] = {
        "particle_porosity": (et - eb) / (1.0 - eb),
        "bulk_density": s["particle_density"] * (1.0 - eb),
        "bed_volume": math.pi * raw["geometry"]["bed_radius"] ** 2 * raw["geometry"]["bed_length"],
        "cross_section": math.pi * raw["geometry"]["bed_radius"] ** 2,
    }
    raw["derived"]["sorbent_mass"] = raw["derived"]["bulk_density"] * raw["derived"]["bed_volume"]
    return Config(raw)


def load_raw_defaults() -> Dict[str, Any]:
    with open(DEFAULTS_PATH, "r") as fh:
        return yaml.safe_load(fh)


def load(user_yaml: Optional[str] = None, overrides: Optional[Sequence[str]] = None,
         override_dict: Optional[Dict[str, Any]] = None) -> Config:
    raw = load_raw_defaults()
    if user_yaml:
        with open(user_yaml, "r") as fh:
            user = yaml.safe_load(fh) or {}
        user.pop("derived", None)      # derived quantities are always recomputed (config_used.yaml has them)
        if "cycle" in user:            # a user cycle replaces the default cycle entirely
            raw["cycle"] = user.pop("cycle")
        raw = deep_merge(raw, user)
    for k, v in (override_dict or {}).items():
        set_dotted(raw, k, v)
    for ov in overrides or []:
        k, v = parse_override(ov)
        set_dotted(raw, k, v)
    return validate(raw)


# ----------------------------------------------------------------------------- docs

def markdown_table() -> str:
    raw = load_raw_defaults()
    lines = ["| key | default | unit | allowed | meaning | source |", "|---|---|---|---|---|---|"]
    for sp in SPECS:
        if sp.key == "cycle":
            continue
        try:
            dv = get_dotted(raw, sp.key)
        except (KeyError, TypeError):
            dv = ""
        if sp.choices:
            allowed = " / ".join(str(c) for c in sp.choices)
        elif sp.lo is not None or sp.hi is not None:
            allowed = f"[{sp.lo}, {sp.hi}]"
        else:
            allowed = sp.type
        esc = lambda x: str(x).replace("|", "\\|")      # noqa: E731  (pipes break markdown tables)
        lines.append(f"| `{sp.key}` | `{dv}` | {esc(sp.unit)} | {esc(allowed)} | {esc(sp.desc)} | {esc(sp.source)} |")
    return "\n".join(lines)


if __name__ == "__main__":
    import sys
    if "--markdown" in sys.argv:
        print(markdown_table())
    else:
        cfg = load()
        print("defaults.yaml is valid. derived:", cfg.derived.to_dict() if hasattr(cfg.derived, 'to_dict') else cfg.derived)
