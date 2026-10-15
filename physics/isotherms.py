"""Equilibrium loadings q* (equilibrium ONLY -- no rates here).

CO2 on Lewatit VP OC 1065: temperature-dependent Toth (Young et al. 2021, eqs 1-4)

    q*  = q_inf(T) b(T) p / (1 + (b(T) p)^tau(T))^(1/tau(T))
    q_inf(T) = q_inf0 exp(chi (1 - T/T0))
    b(T)     = b0 exp(-dH0 / (R T))            (van 't Hoff affinity, no T_ref)
    tau(T)   = tau0 + alpha (1 - T0/T)

H2O: GAB with Anderson temperature dependence (Young 2021 eqs 5-8)

    q*_w = q_m k c x / ((1 - k x)(1 + (c - 1) k x)),   x = p_w / p_sat(T)
    c = exp((E1 - E10)/RT),  k = exp((E2_9 - E10)/RT)
    E1 = C - exp(D T),  E2_9 = F + G T,  E10 = slope T + intercept

CO2-H2O co-adsorption (q_w = adsorbed water, mol/kg):

* ``wadst``  q = (1 - w) Toth_dry + w Toth_wet,  w = exp(-A/q_w)        (Young eq. 14)
* ``mechanistic`` q = (phi/phi_dry) Toth(dH0 -> dH_ave)                  (Young eqs 9-13)
      f_blocked = f_bmax (1 - exp(-(k q_w)^n)); phi_av = 1 - f_blocked
      phi = phi_dry + (phi_av - phi_dry) w;  dH_ave = (1-w) dH_dry + w dH_wet
  (Young et al. state that dH_ave enters 'the affinity'; we evaluate b with it.)

Trace-level regularisation: q*_CO2 is multiplied by p/(p + p_reg) with
p_reg = isotherm.co2.p_regularization (default 0.01 Pa = 0.025 % of the feed CO2
partial pressure).  The WADST wet site (b0_wet = 1.2e-18 1/Pa, dH = -204 kJ/mol,
tau0 = 0.053) loads ~0.6 mol/kg already at 1e-3 Pa (25 ppb), i.e. dq*/dp ~ 1e20
mol/kg/Pa at p -> 0.  That is an artefact of the fit far outside its data range
and makes Newton iterations at the leading edge of the CO2 front ill-posed.
The factor caps the slope, makes q*(p) C1 at p = 0, and changes q* by < 0.03 %
at the feed partial pressure (< 1 % for p > 1 Pa).
* ``stampi_bombelli`` q_inf -> q_inf/(1 - gamma q_w), b -> b (1 + beta q_w)  (Young eqs 15-16)
* ``none``   dry Toth only.

All functions are vectorised and safe for slightly negative inputs produced by
Newton iterations (pressures and loadings are clipped at zero).
"""
from __future__ import annotations

import numpy as np

from .base import PhysicsModule, R_GAS


def vant_hoff(b0, dH, T):
    """Affinity b(T) = b0 exp(-dH/(R T)) with dH < 0 for exothermic adsorption."""
    return b0 * np.exp(-dH / (R_GAS * T))


def toth(p, T, q_inf0, chi, b0, dH0, tau0, alpha, T0):
    p = np.maximum(p, 0.0)
    qinf = q_inf0 * np.exp(chi * (1.0 - T / T0))
    b = vant_hoff(b0, dH0, T)
    tau = tau0 + alpha * (1.0 - T0 / T)
    bp = b * p
    return qinf * bp / (1.0 + bp ** tau) ** (1.0 / tau)


def gab(x, T, q_m, C, D, F, G, E10_slope, E10_intercept):
    E1 = C - np.exp(D * T)
    E29 = F + G * T
    E10 = E10_slope * T + E10_intercept
    RT = R_GAS * T
    c = np.exp((E1 - E10) / RT)
    k = np.exp((E29 - E10) / RT)
    kx = k * x
    return q_m * c * kx / ((1.0 - kx) * (1.0 + (c - 1.0) * kx))


def smooth_cap(x, cap, width=0.02):
    """min(x, cap) smoothed over `width` (C1), identity well below the cap.

    Slightly negative x (Newton iterates) pass through unchanged -- GAB is a
    smooth function there -- so no kink is introduced at x = 0."""
    x = np.maximum(x, -0.2)
    lo = cap - width
    over = np.maximum(x - lo, 0.0)
    # quadratic blend on [lo, lo + 2 width] that ends with slope 0 at cap
    blend = lo + over - over ** 2 / (4.0 * width)
    return np.where(x <= lo, x, np.where(over >= 2.0 * width, cap, blend))


class EquilibriumModel(PhysicsModule):
    name = "isotherms"

    def __init__(self, cfg, grid=None, props=None):
        super().__init__(cfg, grid, props)
        iso = cfg.isotherm
        c = iso.co2
        self.co2 = dict(q_inf0=c.q_inf0, chi=c.chi, b0=c.b0, dH0=c.dH0, tau0=c.tau0, alpha=c.alpha, T0=c.T0)
        h = iso.h2o
        self.h2o = dict(q_m=h.q_m, C=h.C, D=h.D, F=h.F, G=h.G, E10_slope=h.E10_slope,
                        E10_intercept=h.E10_intercept)
        self.rh_cap = h.rh_cap
        self.p_reg = c.p_regularization
        self.water = bool(cfg.components.water)
        self.model = iso.coadsorption.model if self.water else "none"
        self.basis = iso.coadsorption.water_loading_basis
        co = iso.coadsorption
        self.wadst = co.wadst.to_dict() if hasattr(co.wadst, "to_dict") else dict(co.wadst._data)
        self.mech = dict(co.mechanistic._data)
        self.sb = dict(co.stampi_bombelli._data)

    # ------------------------------------------------------------------ pieces
    def q_co2_dry(self, p, T):
        return toth(p, T, **self.co2)

    def q_h2o(self, p_w, T, p_sat):
        x = smooth_cap(p_w / p_sat, self.rh_cap)
        return gab(x, T, **self.h2o)

    def q_co2(self, p, T, qw):
        """CO2 equilibrium loading incl. co-adsorption and trace-level regularisation."""
        pp = np.maximum(p, 0.0)
        return self._q_co2_raw(pp, T, qw) * (pp / (pp + self.p_reg))

    def _q_co2_raw(self, p, T, qw):
        qw = np.maximum(qw, 1e-12)
        if self.model == "none":
            return self.q_co2_dry(p, T)
        w = np.exp(-self.wadst["A"] / qw) if self.model == "wadst" else None
        if self.model == "wadst":
            wd = self.wadst
            q_wet = toth(p, T, wd["q_inf0_wet"], wd["chi_wet"], wd["b0_wet"], wd["dH_wet"],
                         wd["tau0_wet"], wd["alpha_wet"], self.co2["T0"])
            return (1.0 - w) * self.q_co2_dry(p, T) + w * q_wet
        if self.model == "mechanistic":
            m = self.mech
            w = np.exp(-m["A"] / qw)
            f_bl = m["f_blocked_max"] * (1.0 - np.exp(-(m["k"] * qw) ** m["n"]))
            phi_av = 1.0 - f_bl
            phi = m["phi_dry"] + (phi_av - m["phi_dry"]) * w
            dH_ave = (1.0 - w) * self.co2["dH0"] + w * m["dH_wet"]
            pars = dict(self.co2)
            pars["dH0"] = dH_ave
            return (phi / m["phi_dry"]) * toth(p, T, **pars)
        if self.model == "stampi_bombelli":
            s = self.sb
            pars = dict(self.co2)
            pars["q_inf0"] = self.co2["q_inf0"] / (1.0 - s["gamma"] * qw)
            pars["b0"] = self.co2["b0"] * (1.0 + s["beta"] * qw)
            return toth(p, T, **pars)
        raise ValueError(self.model)

    # ------------------------------------------------------------------ main entry
    def equilibrium(self, p_co2, p_h2o, T, q_h2o_state, p_sat):
        """Return q* array (na, ...) ordered [CO2, (H2O)]."""
        if not self.water:
            return self.q_co2(p_co2, T, None)[None]
        qw_eq = self.q_h2o(p_h2o, T, p_sat)
        qw_basis = q_h2o_state if self.basis == "state" else qw_eq
        qc = self.q_co2(p_co2, T, qw_basis)
        return np.stack([qc, qw_eq])
