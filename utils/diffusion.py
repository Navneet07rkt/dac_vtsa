"""Molecular diffusivity (raw diffusion only; dispersion lives in physics/dispersion.py).

Fuller, Schettler & Giddings (1966) correlation as given in Poling, Prausnitz &
O'Connell, *The Properties of Gases and Liquids*, 5th ed., eq. 11-4.4:

    D_AB [cm^2/s] = 0.00143 T^1.75 / ( P[bar] * M_AB^0.5 * (Sv_A^(1/3) + Sv_B^(1/3))^2 )
    M_AB = 2 / (1/M_A + 1/M_B)   [g/mol]

Diffusion volumes (Poling Table 11-1).
"""
from __future__ import annotations

import numpy as np

DIFFUSION_VOLUME = {"CO2": 26.7, "H2O": 13.1, "N2": 18.5, "O2": 16.3, "AIR": 19.7, "Ar": 16.2}
MOLAR_MASS_G = {"CO2": 44.01, "H2O": 18.015, "N2": 28.013, "O2": 31.999, "AIR": 28.96, "Ar": 39.948}


def fuller(T, P, a: str, b: str):
    """Binary diffusivity D_ab [m^2/s] at temperature T [K] and pressure P [Pa]."""
    Ma, Mb = MOLAR_MASS_G[a], MOLAR_MASS_G[b]
    Mab = 2.0 / (1.0 / Ma + 1.0 / Mb)
    s = DIFFUSION_VOLUME[a] ** (1 / 3) + DIFFUSION_VOLUME[b] ** (1 / 3)
    T = np.asarray(T, dtype=float)
    P = np.maximum(np.asarray(P, dtype=float), 1.0)
    return 1.43e-7 * T ** 1.75 / ((P / 1e5) * np.sqrt(Mab) * s * s)


def blanc_mixture(y_other, D_pairs, floor: float = 1e-12):
    """Blanc's law: effective diffusivity of species A in a mixture.

    y_other : list of mole fractions of the *other* species
    D_pairs : list of binary diffusivities D_A,j (same order)
    1/D_A,mix = sum_j (y_j / D_Aj) / sum_j y_j
    """
    num = sum(y_other) + floor
    den = sum(y / D for y, D in zip(y_other, D_pairs)) + floor / D_pairs[0]
    return num / den
