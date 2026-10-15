"""Axial / transverse dispersion correlations (molecular D_m comes from utils/diffusion).

All coefficients are *interstitial* (the flux is  J = -eps_b * D * C * grad(y)).

* Edwards & Richardson (1968):
      D_L = g1 D_m + g2 d_p v / (1 + beta g1 D_m / (d_p v)),  g1 = 0.73, g2 = 0.5, beta = 13
* young2021 (Young et al. 2021 ESI):  D_L = 0.7 D_m + 0.5 d_p v
* Wakao & Kagei (1982):               D_L = D_m/eps_b (20 + 0.5 Re Sc)   (used by Galanti 2025)

Transverse (r, theta):  D_T = g1 D_m + d_p v / Pe_inf   with Pe_inf ~ 8-12
(Delgado 2006 review).  v = |u_s| / eps_b is the interstitial velocity.
"""
from __future__ import annotations

import numpy as np

from .base import PhysicsModule


class Dispersion(PhysicsModule):
    name = "dispersion"

    def __init__(self, cfg, grid=None, props=None):
        super().__init__(cfg, grid, props)
        tr = cfg.transport
        self.model = tr.dispersion_model
        self.g1, self.g2, self.beta = tr.er_gamma1, tr.er_gamma2, tr.er_beta
        self.Pe_T = tr.Pe_mass_radial_inf
        self.dp = cfg.sorbent.particle_diameter
        self.eps = cfg.sorbent.bed_voidage

    def axial(self, Dm, u_s, rho=None, mu=None):
        v = np.abs(u_s) / self.eps
        dp = self.dp
        if self.model == "edwards_richardson":
            return self.g1 * Dm + self.g2 * dp * v / (1.0 + self.beta * self.g1 * Dm / (dp * v + 1e-30))
        if self.model == "young2021":
            return 0.7 * Dm + 0.5 * dp * v
        if self.model == "wakao_kagei":
            Re = rho * np.abs(u_s) * dp / mu
            Sc = mu / (rho * Dm)
            return Dm / self.eps * (20.0 + 0.5 * Re * Sc)
        raise ValueError(self.model)

    def transverse(self, Dm, u_s_nominal):
        v = np.abs(u_s_nominal) / self.eps
        return self.g1 * Dm + self.dp * v / self.Pe_T
