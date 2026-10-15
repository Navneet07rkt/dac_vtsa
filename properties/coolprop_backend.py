"""Real-gas pure-component properties from CoolProp (optional dependency).

CoolProp is far too slow to call per cell per right-hand-side evaluation, so
viscosity and conductivity of each pure fluid are tabulated once on a fine
temperature grid (at 1 atm for CO2/Air; at a sub-saturation pressure for water
vapour so that the vapour branch is used) and linearly interpolated.  Mixing
rules are inherited from :class:`PropertyBackend`.  The constant heat capacity
used in the energy balance is CoolProp's ideal-gas cp0 at cp_reference_T.
Water saturation pressure comes from IAPWS-95 (CoolProp) above 273.16 K.
Binary diffusivities are not available in CoolProp -> Fuller (same as ideal).
"""
from __future__ import annotations

import numpy as np

from .ideal_gas_backend import PropertyBackend, arden_buck

FLUID = {"CO2": "CO2", "H2O": "Water", "N2": "Air"}


class CoolPropBackend(PropertyBackend):
    name = "coolprop"

    def __init__(self, species, cp_reference_T=330.0, T_min=250.0, T_max=520.0, n=271):
        try:
            import CoolProp.CoolProp as CP   # noqa: F401
        except ImportError as exc:   # pragma: no cover
            raise ImportError("properties.backend = coolprop needs `pip install CoolProp`") from exc
        self._CP = CP
        self._Tg = np.linspace(T_min, T_max, n)
        super().__init__(species, cp_reference_T)
        self._mu_tab = {}
        self._k_tab = {}
        for s in self.species:
            mus, ks = [], []
            for T in self._Tg:
                if s == "H2O":
                    T = max(T, 274.0)      # IAPWS-95 not defined below the triple point
                P = self._p_for(s, T)
                mus.append(CP.PropsSI("V", "T", T, "P", P, FLUID[s]))
                ks.append(CP.PropsSI("L", "T", T, "P", P, FLUID[s]))
            self._mu_tab[s] = np.array(mus)
            self._k_tab[s] = np.array(ks)
        Ts = self._Tg[self._Tg > 273.2]
        self._psat_T = Ts
        self._psat = np.array([CP.PropsSI("P", "T", T, "Q", 0, "Water") for T in Ts])

    def _p_for(self, s, T):
        if s == "H2O":
            return max(50.0, 0.3 * float(arden_buck(min(T, 640.0)))) if T < 640 else 101325.0
        return 101325.0

    def _cp_pure(self, s, T):
        return float(self._CP.PropsSI("Cp0molar", "T", T, "P", 101325.0, FLUID[s]))

    def mu_pure(self, T):
        T = np.asarray(T, dtype=float)
        return np.stack([np.interp(T, self._Tg, self._mu_tab[s]) for s in self.species])

    def k_pure(self, T):
        T = np.asarray(T, dtype=float)
        return np.stack([np.interp(T, self._Tg, self._k_tab[s]) for s in self.species])

    def p_sat(self, T):
        T = np.asarray(T, dtype=float)
        hi = np.interp(T, self._psat_T, self._psat)
        return np.where(T > 273.2, hi, arden_buck(T))
