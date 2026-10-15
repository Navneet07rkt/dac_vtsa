"""Momentum balance (Ergun) and valve pressure ramps.

Ergun, written per face in the face-normal direction n (superficial velocity u):

    -dP/dn = a u + b |u| u,   a = A mu (1-eps)^2 / (eps^3 d_p^2),   b = B rho (1-eps) / (eps^3 d_p)

Inverting for u (smooth, no division by b):

    u = -2 g / (a + sqrt(a^2 + 4 b |g|)),   g = dP/dn

Valve opening ("blowdown pressure ramp"): when a valve opens onto a reservoir
at P_target, the face pressure follows

    P_face(t) = P_target + (P_start - P_target) exp(-t / tau_valve)

which models the finite opening/pump-down speed (cf. the exponential pressure
histories of Haghpanah et al. 2013) and keeps the problem well posed.
"""
from __future__ import annotations

import numpy as np

from .base import PhysicsModule


class Ergun(PhysicsModule):
    name = "pressure"

    def __init__(self, cfg, grid=None, props=None):
        super().__init__(cfg, grid, props)
        eps = cfg.sorbent.bed_voidage
        dp = cfg.sorbent.particle_diameter
        self.ka = cfg.transport.ergun_A * (1.0 - eps) ** 2 / (eps ** 3 * dp ** 2)   # * mu
        self.kb = cfg.transport.ergun_B * (1.0 - eps) / (eps ** 3 * dp)            # * rho

    def coefficients(self, mu, rho):
        return self.ka * mu, self.kb * rho

    @staticmethod
    def velocity(g, a, b):
        """Superficial velocity from pressure gradient g = dP/dn (Pa/m)."""
        return -2.0 * g / (a + np.sqrt(a * a + 4.0 * b * np.abs(g)))

    @staticmethod
    def dp_dn(u, a, b):
        """Pressure gradient that drives velocity u (inverse of velocity())."""
        return -(a * u + b * np.abs(u) * u)


def ramp(t, P_start, P_target, tau):
    return P_target + (P_start - P_target) * np.exp(-t / tau)


def smooth_positive(u, delta=1e-7):
    """max(u, 0) smoothed over |u| ~ delta (for check valves)."""
    return 0.5 * (u + np.sqrt(u * u + delta * delta)) - 0.5 * delta
