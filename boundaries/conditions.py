"""Boundary-condition mathematics for the finite-volume model.

Finite volumes need boundary *fluxes*, not ghost-node formulas.  This module
turns the valve state reported by valves/valve_logic.py into face fluxes:

* closed end  (Neumann / zero flux):           N_k = 0, conductive flux = 0
* feed inlet  (Danckwerts, fixed molar flow):   N_k = F y_feed,k,  H = F sum y_k h_k(T_feed)
* vacuum pump port (default for OPEN_VACUUM) -- finite pumping speed S per unit
  bed cross-section k_p = S/A [m/s], ultimate (set) pressure P_vac:
      N_out = k_p (P_cell - P_vac)/(R T_cell)            (check valve: no back-flow)
  The bed pressure then relaxes exponentially towards P_vac with time constant
  ~ eps_t L / k_p: this IS the 'blowdown pressure ramp', derived from the pump
  rather than imposed.  It also keeps the boundary non-stiff, which matters for
  exact mass/energy balance closure (a Dirichlet-P vacuum port is ~1e6 1/s stiff).
* atmospheric vent / feed header (OPEN_ATMOSPHERIC, default) -- linear valve of
  conductance k_v [m/s] to the reservoir at P_atm (or P_feed):
      N_in = k_v (P_res - P_cell)/(R T_upwind)    (+ = into the bed; back-flow allowed)
  k_v = 20 m/s passes the default feed flow with ~0.3 % pressure rise.  A valve
  (instead of fixing P on the face) keeps the Jacobian ~1e3x better conditioned.
* optional Dirichlet pressure port (vent_model / vacuum_model = ramp):
      P_face(t) = P_target + (P_start - P_target) exp(-t/tau)
      u_out     = Ergun^-1( (P_face - P_cell) / (dz/2) )      (+ = leaving the bed)
      N_k       = u_out A c_k,cell          if u_out >= 0   (zero-gradient outlet)
                = u_out A y_ext,k P_face/(R T_ext)  if u_out < 0  (back-flow of external gas)
      optional check valve: u_out -> max(u_out, 0) (smoothed)
* outer wall r = R (Robin):  q = (T_w - T_N) / (1/alpha_w + (R - r_N)/lam_r)
* axis r = 0 (symmetry/regularity): handled geometrically -- zero face area.
* theta (periodic): handled by np.roll in physics/mass.py & energy.py.

Scalar helper functions for the classical BC types are provided for unit tests
and for anyone extending the code with a finite-difference discretisation.
"""
from __future__ import annotations

from typing import Dict

import numpy as np

from numerics.stencils import upwind_weight
from physics.pressure import Ergun, ramp, smooth_positive
from properties import R_GAS

VALVE_DP_SMOOTH = 0.2      # Pa: smoothing of the valve flux splitting (see 'valve')


def softplus(x, width):
    """Smooth max(x, 0) that stays strictly positive: width * log(1 + exp(x/width))."""
    return width * np.logaddexp(0.0, x / width)


# ---------------------------------------------------------------- classical BC helpers

def dirichlet_ghost(value_boundary, value_interior):
    """Ghost value that puts `value_boundary` exactly on a cell face (2nd order)."""
    return 2.0 * value_boundary - value_interior


def neumann_ghost(value_interior, grad, h):
    """Ghost value for prescribed outward gradient `grad` at a face at distance h/2."""
    return value_interior + grad * h


def robin_flux(T_ext, T_int, h, half_dx, lam):
    """Heat flux INTO the domain through a Robin boundary with resistance 1/h + half_dx/lam."""
    return (T_ext - T_int) / (1.0 / h + half_dx / lam)


def periodic_neighbours(a, axis):
    """(previous, next) neighbours along a periodic axis."""
    return np.roll(a, 1, axis=axis), np.roll(a, -1, axis=axis)


# ---------------------------------------------------------------- end-face fluxes

def end_face_fluxes(cond, t: float, P_start: float, layer: Dict[str, np.ndarray], Az: np.ndarray,
                    half_dz: float, ergun: Ergun, cpg: np.ndarray, T_ref: float, feed_flux: np.ndarray,
                    tau_valve: float) -> Dict[str, np.ndarray]:
    """Fluxes through one end of the column.

    layer : c (ng,Nr,Nth), T, P, mu, rho (Nr,Nth) of the boundary cell layer
    Az    : (Nr,1) axial face areas;  feed_flux: (Nr,Nth) total molar feed per face [mol/s]
    returns F_in (ng,Nr,Nth) [mol/s INTO bed], H_in (Nr,Nth) [W INTO bed],
            u_out (Nr,Nth) superficial velocity (+ = leaving), P_face (Nr,Nth)
    """
    c, T, P = layer["c"], layer["T"], layer["P"]
    ng = c.shape[0]
    if cond.kind == "closed":
        z = np.zeros_like(T)
        return {"F_in": np.zeros_like(c), "H_in": z, "u_out": z.copy(), "P_face": P.copy()}
    a, b = ergun.coefficients(layer["mu"], layer["rho"])
    if cond.kind == "feed":
        y = cond.y_ext.reshape(ng, 1, 1)
        F_in = feed_flux[None] * y
        H_in = np.tensordot(cpg, F_in, axes=(0, 0)) * (cond.T_ext - T_ref)
        C_cell = np.maximum(c.sum(axis=0), 1e-12)
        u_in = feed_flux / (Az * C_cell)
        P_face = P + (a * u_in + b * u_in * u_in) * half_dz
        return {"F_in": F_in, "H_in": H_in, "u_out": -u_in, "P_face": P_face}
    if cond.kind == "pump":
        dp = P - cond.P_target
        C_cell = np.maximum(c.sum(axis=0), 1e-12)
        y_cell = c / C_cell
        if cond.check_valve:
            # check valve: flow only OUT of the bed.  softplus keeps N_out > 0 and smooth (width
            # 1 Pa), so the stream always carries the bed composition -- no blending with the gas
            # of the vacuum line (a velocity-based upwind switch would blend for most of the
            # desorption step: the pump velocity k_p dP/P is only ~1e-5..1e-3 m/s).
            dp = softplus(dp, 1.0)
            w_out = np.ones_like(dp)
        else:
            w_out = 0.5 * (1.0 + np.tanh(dp / VALVE_DP_SMOOTH))
        N_out = cond.speed * dp / (R_GAS * T)                  # mol/m2/s (+ = leaving)
        y_use = w_out[None] * y_cell + (1.0 - w_out[None]) * cond.y_ext.reshape(ng, 1, 1)
        F_out = (N_out * Az)[None] * y_use
        T_up = w_out * T + (1.0 - w_out) * cond.T_ext
        H_out = np.tensordot(cpg, F_out, axes=(0, 0)) * (T_up - T_ref)
        u_out = N_out / C_cell
        P_face = P - (a * u_out + b * np.abs(u_out) * u_out) * half_dz
        return {"F_in": -F_out, "H_in": -H_out, "u_out": u_out, "P_face": P_face}
    if cond.kind == "valve":
        # linear valve of conductance k_v to a reservoir at P_target (+ = into the bed); the
        # valve opens over ~tau_valve: the reservoir pressure seen by the bed ramps from the
        # bed pressure at step start (no impulsive pressure wave when a vent opens on a bed
        # under vacuum -- physically a valve opening time, numerically a smooth transient)
        P_res = ramp(t, P_start, cond.P_target, tau_valve) if tau_valve > 0 else cond.P_target
        dp = P_res - P
        C_cell = np.maximum(c.sum(axis=0), 1e-12)
        # FLUX SPLITTING.  N = N_in - N_out with N_in = k_v softplus(dp)/(R T_ext) carrying reservoir
        # gas at T_ext and N_out = k_v softplus(-dp)/(R T) carrying bed gas at the bed T.  Each part
        # is upwinded exactly, and the split is smooth for Newton.  A blended composition/temperature
        # (one stream with T_up = w T_ext + (1-w) T) is NOT used: near dp = 0 it lets out-flowing gas
        # carry less enthalpy than h(T_bed), so a warmer cell pushes out gas and heats up further --
        # a spurious instability (+15 1/s) found at the end of pressurisation.
        p_in = softplus(dp, VALVE_DP_SMOOTH)
        p_out = softplus(-dp, VALVE_DP_SMOOTH)               # p_in - p_out == dp exactly
        N_in = cond.speed * p_in / (R_GAS * cond.T_ext)        # mol/m2/s, >= 0
        N_out = cond.speed * p_out / (R_GAS * T)               # mol/m2/s, >= 0
        y_ext = cond.y_ext.reshape(ng, 1, 1)
        y_cell = c / C_cell
        F_inpart = (N_in * Az)[None] * y_ext
        F_outpart = (N_out * Az)[None] * y_cell
        H_inpart = np.tensordot(cpg, F_inpart, axes=(0, 0)) * (cond.T_ext - T_ref)
        H_outpart = np.tensordot(cpg, F_outpart, axes=(0, 0)) * (T - T_ref)
        N_net = N_in - N_out
        u_out = -N_net / C_cell
        P_face = P - (a * u_out + b * np.abs(u_out) * u_out) * half_dz
        return {"F_in": F_inpart - F_outpart, "H_in": H_inpart - H_outpart, "u_out": u_out, "P_face": P_face,
                "F_inpart": F_inpart, "F_outpart": F_outpart, "H_inpart": H_inpart, "H_outpart": H_outpart}
    # pressure-controlled port (Dirichlet on P with an exponential ramp)
    P_face = ramp(t, P_start, cond.P_target, tau_valve) * np.ones_like(P)
    g_out = (P_face - P) / half_dz
    u_out = Ergun.velocity(g_out, a, b)
    if cond.check_valve:
        u_out = smooth_positive(u_out)
    c_ext = cond.y_ext.reshape(ng, 1, 1) * P_face[None] / (R_GAS * cond.T_ext)
    w_out = upwind_weight(u_out)
    F_out = u_out[None] * Az[None] * (w_out[None] * c + (1.0 - w_out[None]) * c_ext)
    T_up = w_out * T + (1.0 - w_out) * cond.T_ext
    H_out = np.tensordot(cpg, F_out, axes=(0, 0)) * (T_up - T_ref)
    return {"F_in": -F_out, "H_in": -H_out, "u_out": u_out, "P_face": P_face}


def wall_to_bed(Tw, T_outer, A, R_bed):
    """Heat flow [W] from the wall ring into the outermost bed cells (Robin)."""
    return A * (Tw - T_outer) / R_bed
