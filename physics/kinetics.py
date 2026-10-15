"""Linear-driving-force uptake rates (rate ONLY -- equilibrium lives in isotherms.py).

    dq_i/dt = k_i(T) (q*_i - q_i),    k_i(T) = k_i,ref exp(-Ea_i/R (1/T - 1/T_ref))

Young et al. 2021 fitted temperature-independent k (Ea = 0) to breakthrough
curves: k_CO2 = 0.003 1/s, k_H2O = 0.0086 1/s (ESI Table S9).  The Arrhenius
form is provided because amine sorbents are known to speed up at desorption
temperatures (k(T) measurements for Lewatit at 20-70 C exist, e.g. the 2024
ChemRxiv study of Lewatit/A110/TIFSIX kinetics); set Ea > 0 to explore it.
"""
from __future__ import annotations

import numpy as np

from .base import PhysicsModule, R_GAS


class LDFKinetics(PhysicsModule):
    name = "kinetics"

    def __init__(self, cfg, grid=None, props=None):
        super().__init__(cfg, grid, props)
        k = cfg.kinetics
        self.water = bool(cfg.components.water)
        self.k_ref = np.array([k.k_ldf_co2, k.k_ldf_h2o]) if self.water else np.array([k.k_ldf_co2])
        self.Ea = np.array([k.Ea_co2, k.Ea_h2o]) if self.water else np.array([k.Ea_co2])
        self.T_ref = k.T_ref

    def k_ldf(self, T):
        """(na, ...) LDF coefficients at temperature T."""
        shape = (len(self.k_ref),) + (1,) * np.ndim(T)
        kr = self.k_ref.reshape(shape)
        if np.all(self.Ea == 0.0):
            return kr * np.ones_like(T)
        Ea = self.Ea.reshape(shape)
        return kr * np.exp(-Ea / R_GAS * (1.0 / T - 1.0 / self.T_ref))

    def rates(self, qstar, q, T):
        return self.k_ldf(T) * (qstar - q)
