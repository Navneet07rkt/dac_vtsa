"""Effective bed conductivity and bed-to-wall heat transfer (suggested addition).

Stagnant effective conductivity -- Zehner & Schlünder (1970), VDI Heat Atlas D6.3:

    lam0/lam_g = 1 - sqrt(1-eps) + 2 sqrt(1-eps)/(1 - B/kap) *
                 [ (1 - 1/kap) B/(1 - B/kap)^2 ln(kap/B) - (B+1)/2 - (B-1)/(1 - B/kap) ]
    B = Cf ((1-eps)/eps)^(10/9),  kap = lam_s/lam_g,  Cf = 1.25 (spheres)

Dynamic (flow) contributions:  lam_r = lam0 + C cp u d_p / Pe_h,r ;  lam_z = lam0 + C cp |u| d_p / Pe_h,z

Bed-to-wall coefficient -- Martin & Nilles (1993), as used in Galanti et al. 2025 SI:

    Nu_w = (1.3 + 5/N) lam0/lam_g + 1 / (1/(0.3 Pr^(1/3) Re^0.75) + 1/(0.054 Pr Re)),
    N = D_bed / d_p,  alpha_w = Nu_w lam_g / d_p
"""
from __future__ import annotations

import numpy as np

from .base import PhysicsModule


def zehner_schlunder(lam_g, lam_s, eps, Cf=1.25):
    kap = lam_s / lam_g
    B = Cf * ((1.0 - eps) / eps) ** (10.0 / 9.0)
    s = np.sqrt(1.0 - eps)
    one_m = 1.0 - B / kap
    kc = (2.0 / one_m) * ((1.0 - 1.0 / kap) * B / one_m ** 2 * np.log(kap / B)
                          - 0.5 * (B + 1.0) - (B - 1.0) / one_m)
    return lam_g * (1.0 - s + s * kc)


class HeatTransfer(PhysicsModule):
    name = "heat_transfer"

    def __init__(self, cfg, grid=None, props=None):
        super().__init__(cfg, grid, props)
        tr = cfg.transport
        self.model_bed = tr.bed_conductivity_model
        self.lam_const = tr.bed_conductivity_constant
        self.Cf = tr.zs_shape_factor
        self.lam_s = cfg.sorbent.particle_conductivity
        self.eps = cfg.sorbent.bed_voidage
        self.dp = cfg.sorbent.particle_diameter
        self.Pe_r = tr.Pe_heat_radial_inf
        self.Pe_z = tr.Pe_heat_axial_inf
        self.wall_model = tr.wall_htc_model
        self.alpha_const = tr.wall_htc_constant
        self.N = 2.0 * cfg.geometry.bed_radius / self.dp

    def stagnant(self, lam_g):
        if self.model_bed == "constant":
            return self.lam_const * np.ones_like(lam_g)
        if self.model_bed == "parallel":     # Zhang et al. 2025 eq. 7 (upper Wiener bound)
            return self.eps * lam_g + (1.0 - self.eps) * self.lam_s
        return zehner_schlunder(lam_g, self.lam_s, self.eps, self.Cf)

    def radial(self, lam0, C, cp, u_nom):
        return lam0 + C * cp * np.abs(u_nom) * self.dp / self.Pe_r

    def axial(self, lam0, C, cp, u_s):
        return lam0 + C * cp * np.abs(u_s) * self.dp / self.Pe_z

    def wall_coefficient(self, lam0, lam_g, rho, mu, cp_mass, u_nom):
        if self.wall_model == "constant":
            return self.alpha_const * np.ones_like(lam0)
        Re = rho * np.abs(u_nom) * self.dp / mu
        Pr = cp_mass * mu / lam_g
        nu0 = (1.3 + 5.0 / self.N) * lam0 / lam_g
        if np.all(Re <= 0):
            dyn = 0.0
        else:
            a = 0.3 * Pr ** (1.0 / 3.0) * np.maximum(Re, 1e-12) ** 0.75
            b = 0.054 * Pr * np.maximum(Re, 1e-12)
            dyn = 1.0 / (1.0 / a + 1.0 / b)
        return (nu0 + dyn) * lam_g / self.dp
