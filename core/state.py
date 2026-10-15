"""CycleState: the container handed from step to step (C, q, T, P, t).

Fields are cell averages on the (Nr, Ntheta, Nz) grid:
    c  (ng, Nr, Nth, Nz)  gas concentrations of [CO2, (H2O), N2]   [mol/m3]
    q  (na, Nr, Nth, Nz)  adsorbed loadings of [CO2, (H2O)]        [mol/kg]
    T  (Nr, Nth, Nz)      bed temperature                           [K]
    Tw (Nth, Nz)          wall temperature                          [K]
    P  derived: sum(c) R T                                          [Pa]
    t  cumulative process time                                      [s]
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple

import numpy as np

from properties import R_GAS


@dataclass
class CycleState:
    c: np.ndarray
    q: np.ndarray
    T: np.ndarray
    Tw: np.ndarray
    t: float = 0.0
    cycle: int = 1
    step: str = "initial"
    species: Tuple[str, ...] = ("CO2", "H2O", "N2")
    adsorbed: Tuple[str, ...] = ("CO2", "H2O")
    meta: Dict = field(default_factory=dict)

    @property
    def P(self) -> np.ndarray:
        return self.c.sum(axis=0) * R_GAS * self.T

    @property
    def y(self) -> np.ndarray:
        return self.c / np.maximum(self.c.sum(axis=0), 1e-30)

    def gas(self, name: str) -> np.ndarray:
        return self.c[self.species.index(name)]

    def loading(self, name: str) -> np.ndarray:
        return self.q[self.adsorbed.index(name)]

    def copy(self) -> "CycleState":
        return copy.deepcopy(self)


def _air_state(cfg, grid, props, species, T, P, y):
    ng = len(species)
    C = P / (R_GAS * T)
    c = np.empty((ng,) + grid.shape)
    for k in range(ng):
        c[k] = C * y[k]
    return c


def initial_state(cfg, grid, props, eq_model=None) -> CycleState:
    """Initial bed state according to run.initial_state."""
    from valves.valve_logic import humid_air_composition
    species = tuple(props.species)
    adsorbed = ("CO2", "H2O") if cfg.components.water else ("CO2",)
    mode = cfg.run.initial_state
    if mode == "csv":
        from simio.checkpoint import load_state_csv
        st = load_state_csv(cfg.run.initial_state_csv, grid, species, adsorbed)
        st.step, st.cycle = "initial", 1
        return st
    T0, P0 = cfg.feed.T, cfg.feed.P
    y = humid_air_composition(species, P0, T0, cfg.feed.RH if cfg.components.water else 0.0,
                              cfg.feed.y_co2_dry, float(props.p_sat(T0)))
    c = _air_state(cfg, grid, props, species, T0, P0, y)
    T = np.full(grid.shape, T0)
    q = np.zeros((len(adsorbed),) + grid.shape)
    if mode == "feed_equilibrium":
        if eq_model is None:
            from physics.isotherms import EquilibriumModel
            eq_model = EquilibriumModel(cfg)
        p = c * R_GAS * T
        iw = species.index("H2O") if "H2O" in species else None
        psat = props.p_sat(T)
        qw = None
        if iw is not None:
            qw = eq_model.q_h2o(p[iw], T, psat)
        qs = eq_model.equilibrium(p[0], p[iw] if iw is not None else None, T, qw, psat)
        q = np.array(qs)
    Tw = np.full(grid.wall_shape, T0)
    return CycleState(c=c, q=q, T=T, Tw=Tw, t=0.0, cycle=1, step="initial",
                      species=species, adsorbed=adsorbed)
