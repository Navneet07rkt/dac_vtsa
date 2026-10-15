"""Ideal-gas property backend (fast, dependency-free).

Pure-component data
-------------------
* heat capacity  cp_i(T) = A T^2 + B T + C       Young et al. 2021 ESI Table S6
  (evaluated once at ``properties.cp_reference_T``; the gas holds < 0.1 % of the
  bed heat capacity, so a constant cp is used in the conservative energy balance)
* viscosity      mu_i(T) = D T + E               Young et al. 2021 ESI Table S7
* conductivity   k_i(T)  = k_300 (T/300)^n        fitted to Incropera & DeWitt Table A.4
                                                  (H2O: IAPWS dilute-gas values)
* water vapour pressure: Arden Buck equation (as in Young et al. 2021 ESI)

The inert species is called ``N2`` in the state vector but represents dry air
without CO2 (0.79 N2 + 0.21 O2); its properties are mole-averaged accordingly.

Mixture rules: Wilke (viscosity), Mason-Saxena (conductivity), molar average (cp).
"""
from __future__ import annotations

import numpy as np

from utils.diffusion import fuller

R_GAS = 8.314462618

CP_POLY = {  # A, B, C  (J/mol/K)  Young 2021 Table S6
    "N2": (1.06569e-5, -6.4057779e-3, 30.1318),
    "CO2": (-3.23753e-5, 6.17517420e-2, 21.9165),
    "H2O": (1.72619108e-3, -1.22971184026, 250.9868),
    "O2": (9.08259e-6, 1.00528795e-3, 28.3000),
}
MU_LIN = {  # D, E  (Pa s)  Young 2021 Table S7
    "N2": (3.96929e-8, 6.15036e-6),
    "CO2": (4.41697e-8, 1.90440e-6),
    "H2O": (3.55595e-8, -9.930767e-7),
    "O2": (4.75391e-8, 6.53871e-6),
}
K_POW = {  # k_300 (W/m/K), exponent
    "N2": (0.02598, 0.81),
    "O2": (0.02676, 0.84),
    "CO2": (0.01655, 1.335),
    "H2O": (0.0186, 1.18),
}
MOLAR_MASS = {"CO2": 0.04401, "H2O": 0.018015, "N2": 0.028960}   # inert = air (no CO2)
AIR_FRAC = {"N2": 0.79, "O2": 0.21}


def arden_buck(T):
    """Saturation pressure of water over liquid [Pa]."""
    Tc = np.asarray(T, dtype=float) - 273.15
    return 611.21 * np.exp((18.678 - Tc / 234.5) * (Tc / (257.14 + Tc)))


class PropertyBackend:
    """Common mixing rules; subclasses supply pure-component functions."""

    name = "base"

    def __init__(self, species, cp_reference_T=330.0):
        self.species = list(species)
        self.ng = len(self.species)
        self.M = np.array([MOLAR_MASS[s] for s in self.species])
        self.cp = np.array([self._cp_pure(s, cp_reference_T) for s in self.species])
        self.i_co2 = self.species.index("CO2")
        self.i_h2o = self.species.index("H2O") if "H2O" in self.species else None
        self.i_inert = self.species.index("N2")
        # Wilke Phi_ij needs mass ratios only once
        Mi = self.M[:, None]
        Mj = self.M[None, :]
        self._wilke_den = np.sqrt(8.0 * (1.0 + Mi / Mj))
        self._mass_q = (Mj / Mi) ** 0.25

    # -------- pure-component hooks
    def _cp_pure(self, s, T):
        raise NotImplementedError

    def mu_pure(self, T):
        raise NotImplementedError

    def k_pure(self, T):
        raise NotImplementedError

    def p_sat(self, T):
        return arden_buck(T)

    # -------- mixtures
    def _phi(self, mu):
        # mu: (ng, ...) -> Phi (ng, ng, ...)
        ratio = np.sqrt(mu[:, None] / mu[None, :])
        mq = self._mass_q.reshape(self.ng, self.ng, *([1] * (mu.ndim - 1)))
        den = self._wilke_den.reshape(self.ng, self.ng, *([1] * (mu.ndim - 1)))
        return (1.0 + ratio * mq) ** 2 / den

    def mu_k_mix(self, T, y):
        """Wilke viscosity and Mason-Saxena conductivity of the mixture."""
        mu = self.mu_pure(T)
        k = self.k_pure(T)
        ng = self.ng
        sq = np.sqrt(mu)
        mu_m = 0.0
        k_m = 0.0
        for i in range(ng):
            den = 0.0
            for j in range(ng):
                if i == j:
                    den = den + y[j]
                else:
                    phi = (1.0 + (sq[i] / sq[j]) * self._mass_q[i, j]) ** 2 / self._wilke_den[i, j]
                    den = den + y[j] * phi
            den = np.maximum(den, 1e-30)
            mu_m = mu_m + y[i] * mu[i] / den
            k_m = k_m + y[i] * k[i] / den
        return mu_m, k_m

    def cp_mix(self, y):
        out = self.cp[0] * y[0]
        for k in range(1, self.ng):
            out = out + self.cp[k] * y[k]
        return out

    def molar_mass_mix(self, y):
        out = self.M[0] * y[0]
        for k in range(1, self.ng):
            out = out + self.M[k] * y[k]
        return out

    def D_co2_mix(self, T, P, y):
        """CO2 diffusivity in the local mixture (Blanc's law over air and H2O)."""
        D_ca = fuller(T, P, "CO2", "AIR")
        yi = y[self.i_inert]
        if self.i_h2o is None:
            return D_ca
        D_cw = fuller(T, P, "CO2", "H2O")
        yw = y[self.i_h2o]
        num = yi + yw + 1e-9
        den = yi / D_ca + yw / D_cw + 1e-9 / D_ca
        return num / den


class IdealGasBackend(PropertyBackend):
    name = "ideal_gas"

    def _cp_pure(self, s, T):
        if s == "N2":   # inert air
            return sum(AIR_FRAC[c] * np.polyval(CP_POLY[c], T) for c in AIR_FRAC)
        return float(np.polyval(CP_POLY[s], T))

    def _mu1(self, s, T):
        if s == "N2":
            return sum(AIR_FRAC[c] * (MU_LIN[c][0] * T + MU_LIN[c][1]) for c in AIR_FRAC)
        a, b = MU_LIN[s]
        return a * T + b

    def _k1(self, s, T):
        if s == "N2":
            return sum(AIR_FRAC[c] * K_POW[c][0] * (T / 300.0) ** K_POW[c][1] for c in AIR_FRAC)
        k0, n = K_POW[s]
        return k0 * (T / 300.0) ** n

    def mu_pure(self, T):
        T = np.asarray(T, dtype=float)
        return np.stack([self._mu1(s, T) for s in self.species])

    def k_pure(self, T):
        T = np.asarray(T, dtype=float)
        return np.stack([self._k1(s, T) for s in self.species])
