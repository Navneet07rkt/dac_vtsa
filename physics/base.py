"""Shared interface for physics modules.

Every physics module is a small object constructed once per step from
``(cfg, grid, props)`` and exposing pure functions of the local field.  The
semi-discrete model (core/residual.py) evaluates the shared intermediate
quantities exactly once per right-hand-side call and stores them in a
:class:`FieldContext`, which is then handed to each module.  This keeps
"how much CO2 the bed wants to hold" (isotherms), "how fast it gets there"
(kinetics), transport (mass/energy/pressure) and correlations (dispersion,
heat transfer) independent and swappable.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Optional

import numpy as np

R_GAS = 8.314462618


class PhysicsModule:
    """Base class: stores config/grid/properties; subclasses add methods."""

    name = "physics"

    def __init__(self, cfg, grid=None, props=None):
        self.cfg = cfg
        self.grid = grid
        self.props = props

    def describe(self) -> str:   # pragma: no cover - cosmetic
        return self.name


@dataclass
class FieldContext:
    """Intermediate quantities of one RHS evaluation (all cell arrays (Nr,Nth,Nz))."""
    t: float
    c: np.ndarray            # (ng, Nr, Nth, Nz) gas concentrations [mol/m3]
    q: np.ndarray            # (na, ...) loadings [mol/kg]
    e: np.ndarray            # energy density [J/m3 bed]
    Tw: Optional[np.ndarray]  # wall temperature (Nth, Nz)
    T: np.ndarray = None
    Ctot: np.ndarray = None
    P: np.ndarray = None
    y: np.ndarray = None
    mu: np.ndarray = None
    kg: np.ndarray = None
    rho: np.ndarray = None
    cp_mix: np.ndarray = None
    Dm: np.ndarray = None
    lam0: np.ndarray = None
    qstar: np.ndarray = None
    rate: np.ndarray = None   # dq/dt (na, ...)
    extra: Dict[str, Any] = field(default_factory=dict)
