"""Pseudo-homogeneous energy balance in CONSERVATIVE form (bed) + wall.

State variable: e = internal energy per unit bed volume [J/m3]

    e = eps_t sum_k c_k (cp_k (T - T_ref) - R T)                 gas  (u = h - RT)
      + rho_b [ cp_s (T - T_ref) + sum_j q_j (cp_a,j (T - T_ref) + dH_j) ]   solid + adsorbed

    d(e V)/dt = - sum_faces [ sum_k h_k(T_f) N_k  - lam grad T ] . n A  + Q_wall + Q_hx

with h_k(T) = cp_k (T - T_ref).  The heat of adsorption is NOT a source term:
it is released automatically because an adsorbed mole carries enthalpy
h_k + dH_j < h_k (dH_j < 0).  Pressure work (the eps dP/dt term of temperature
-form models) and the (-dH - RT) constant-volume adsorption energy are likewise
contained in the definition of e.  Consequences:

* T is an algebraic, linear function of (e, c, q) -> :meth:`Thermo.temperature`;
* bed + wall energy is conserved to solver tolerance, so the energy balance
  closure reported in results/ is a real check, not an approximation.

Using cp_a,H2O = 75.3 J/mol/K (liquid-like adsorbed water) makes
dH_H2O(T) = dH_H2O(T_ref) + (cp_a - cp_g)(T - T_ref) follow Kirchhoff's law
(the heat of vaporisation of water falls by ~3.5 kJ/mol between 15 and 100 C).
"""
from __future__ import annotations

import numpy as np

from .base import PhysicsModule, R_GAS


class Thermo(PhysicsModule):
    name = "energy"

    def __init__(self, cfg, grid=None, props=None):
        super().__init__(cfg, grid, props)
        self.eps_t = cfg.sorbent.total_voidage
        self.rho_b = cfg.derived.bulk_density
        self.cp_s = cfg.sorbent.cp_solid
        self.cpg = np.asarray(props.cp, dtype=float)
        water = bool(cfg.components.water)
        th = cfg.thermo
        self.cpa = np.array([th.cp_ads_co2, th.cp_ads_h2o] if water else [th.cp_ads_co2])
        self.dH = np.array([th.dH_ads_co2, th.dH_ads_h2o] if water else [th.dH_ads_co2])
        self.T_ref = th.T_ref

    def _terms(self, c, q):
        Cg = c[0] + c[1] if len(c) == 2 else c[0] + c[1] + c[2]
        gas_cp = _wsum(self.cpg, c)
        A = self.eps_t * gas_cp + self.rho_b * (self.cp_s + _wsum(self.cpa, q))
        B = self.eps_t * R_GAS * Cg
        Hq = self.rho_b * _wsum(self.dH, q)
        return A, B, Hq

    def temperature(self, e, c, q):
        A, B, Hq = self._terms(c, q)
        return (e - Hq + A * self.T_ref) / (A - B)

    def energy(self, T, c, q):
        A, B, Hq = self._terms(c, q)
        return A * (T - self.T_ref) - B * T + Hq

    def heat_capacity(self, c, q):
        """Volumetric heat capacity at constant composition, de/dT [J/m3/K]."""
        A, B, _ = self._terms(c, q)
        return A - B

    def gas_enthalpy(self, T):
        """h_k(T) per mole, shape (ng, ...)."""
        shape = (len(self.cpg),) + (1,) * np.ndim(T)
        return self.cpg.reshape(shape) * (T - self.T_ref)


def _wsum(w, a):
    """sum_k w[k] a[k] for a small leading dimension (faster than tensordot)."""
    out = w[0] * a[0]
    for k in range(1, len(w)):
        out = out + w[k] * a[k]
    return out


def enthalpy_flux(F, T_face, cpg, T_ref):
    """sum_k F_k cp_k (T_face - T_ref) for species fluxes F (ng, ...)."""
    return _wsum(cpg, F) * (T_face - T_ref)


class Wall(PhysicsModule):
    """Thin metal wall ring, lumped in r, resolved in (theta, z)."""
    name = "wall"

    def __init__(self, cfg, grid, props=None):
        super().__init__(cfg, grid, props)
        w = cfg.wall
        self.enabled = bool(w.enabled)
        self.rho, self.cp, self.lam = w.density, w.cp, w.conductivity
        self.C = self.rho * self.cp * grid.V_wall          # J/K per wall cell
        self.T_ref = cfg.thermo.T_ref

    def conduction(self, Tw):
        """Net conductive heat INTO each wall cell from wall neighbours [W], (Nth, Nz)."""
        g = self.grid
        out = np.zeros_like(Tw)
        if g.Nz > 1:
            kz = self.lam * g.A_wall_axial / g.dz
            fz = kz * np.diff(Tw, axis=1)            # flux from k to k+1 is -fz
            out[:, :-1] += fz
            out[:, 1:] -= fz
        if g.Nth > 1:
            kt = self.lam * g.wall_thickness * g.dz / (g.R_wall_mid * g.dth)
            n = Tw.shape[0]
            ft = kt * (np.take(Tw, np.r_[1:n, 0], axis=0) - Tw)
            out += ft
            out -= np.take(ft, np.r_[n - 1, 0:n - 1], axis=0)
        return out

    def energy(self, Tw):
        return float(np.sum(self.C * (Tw - self.T_ref)))
