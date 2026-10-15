"""Valve states and their mapping to end-of-column boundary conditions.

Valves: V1 = inlet (top), V2 = bottom outlet (vent / blowdown), V3 = product
outlet (end configurable).  Each valve carries one of FOUR states (the
roadmap's three outlet states plus the feed blower):

=================  ==========================================  =============================================
state              meaning                                     boundary condition produced here
=================  ==========================================  =============================================
CLOSED             sealed line, no flow                        zero total flux (u = 0), zero dispersive flux,
                                                               adiabatic end
OPEN               (V1 only) feed blower, fixed molar flow     Danckwerts inlet: N_k = F_feed y_feed,k
OPEN_ATMOSPHERIC   passive vent / header at P_atm              linear valve N = k_v (P_atm - P)/(R T) (default)
                                                               or Dirichlet P ramp; outflow carries bed gas,
                                                               back-flow carries ambient (feed, for V1) gas
OPEN_VACUUM        pump line at P_vacuum                       vacuum pump of finite speed S (default):
                                                               N_out = (S/A)(P - P_vac)/(R T), convective
                                                               outflow, check valve; or Dirichlet P ramp
=================  ==========================================  =============================================

NOTE (deliberate deviation from the roadmap): the roadmap maps OPEN_VACUUM to
"Dirichlet C = 0 at the outlet".  That is not what a vacuum pump does -- it
fixes the TOTAL pressure.  During heating the gas leaving V3 is ~pure CO2+H2O,
so imposing C_CO2 = 0 would create an artificial diffusive sink (CO2 'leaving'
by dispersion into a zero-concentration ghost), break the CO2 balance and
over-predict desorption.  Here the composition at an outlet is convected out
(zero-gradient / Danckwerts outlet) and only P is imposed.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Dict, Optional

import numpy as np

from properties import R_GAS


class ValveState(str, Enum):
    CLOSED = "CLOSED"
    OPEN = "OPEN"
    OPEN_ATMOSPHERIC = "OPEN_ATMOSPHERIC"
    OPEN_VACUUM = "OPEN_VACUUM"


@dataclass
class EndCondition:
    end: str                      # 'top' (z = 0) or 'bottom' (z = L)
    valve: Optional[str]          # valve open at this end (None if all closed)
    state: ValveState
    kind: str                     # 'closed' | 'feed' | 'pressure'
    P_target: float = 101325.0
    y_ext: Optional[np.ndarray] = None   # composition of gas entering from outside
    T_ext: float = 298.15
    check_valve: bool = False
    is_product: bool = False
    is_vacuum: bool = False
    is_feed: bool = False
    speed: float = 0.0             # vacuum pump speed per unit cross-section [m/s] (kind == 'pump')


def humid_air_composition(species, P, T, RH, y_co2_dry, p_sat):
    """Mole fractions [CO2, (H2O), N2] of air at P, T, RH (y_co2_dry on dry basis)."""
    y_w = RH * p_sat / P if "H2O" in species else 0.0
    y = {"H2O": y_w, "CO2": y_co2_dry * (1.0 - y_w)}
    y["N2"] = 1.0 - y["CO2"] - y_w
    return np.array([y[s] for s in species])


def end_conditions(cfg, step: Dict, species, props) -> Dict[str, EndCondition]:
    """Boundary condition at the top and bottom ends for one cycle step."""
    feed = cfg.feed
    amb = cfg.ambient
    y_feed = humid_air_composition(species, feed.P, feed.T, feed.RH if "H2O" in species else 0.0,
                                   feed.y_co2_dry, float(props.p_sat(feed.T)))
    y_amb = humid_air_composition(species, amb.P, amb.T, amb.RH if "H2O" in species else 0.0,
                                  amb.y_co2_dry, float(props.p_sat(amb.T)))
    ports = cfg.ports
    out = {}
    for end in ("top", "bottom"):
        valve = None
        for v in ("V1", "V2", "V3"):
            if ports[v] == end and step["valves"][v] != "CLOSED":
                valve = v
        if valve is None:
            out[end] = EndCondition(end, None, ValveState.CLOSED, "closed")
            continue
        st = ValveState(step["valves"][valve])
        is_prod = (step.get("product_port", "V3") == valve)
        if st == ValveState.OPEN:
            out[end] = EndCondition(end, valve, st, "feed", feed.P, y_feed, feed.T,
                                    is_product=False, is_feed=True)
        elif st == ValveState.OPEN_ATMOSPHERIC:
            kind = "valve" if cfg.operation.vent_model == "valve" else "pressure"
            kv = cfg.operation.vent_conductance
            if valve == "V1":    # feed header at feed pressure (pressurisation with feed gas)
                out[end] = EndCondition(end, valve, st, kind, feed.P, y_feed, feed.T,
                                        is_product=is_prod, is_feed=True, speed=kv)
            else:
                out[end] = EndCondition(end, valve, st, kind, amb.P, y_amb, amb.T,
                                        is_product=is_prod, speed=kv)
        else:   # OPEN_VACUUM
            kind = "pump" if cfg.operation.vacuum_model == "pump" else "pressure"
            out[end] = EndCondition(end, valve, st, kind, cfg.operation.P_vacuum, y_amb, amb.T,
                                    check_valve=bool(cfg.operation.vacuum_check_valve),
                                    is_product=is_prod, is_vacuum=True,
                                    speed=cfg.operation.pump_speed_per_area)
    return out


def feed_molar_flux(cfg):
    """Feed molar flux per unit cross-section [mol/m2/s] at feed T, P."""
    return cfg.feed.velocity * cfg.feed.P / (R_GAS * cfg.feed.T)


def feed_profile(cfg, grid) -> np.ndarray:
    """Normalised inlet distribution f(r, theta) with area-weighted mean 1, shape (Nr, Nth)."""
    prof = cfg.feed.profile
    A = cfg.feed.profile_amplitude
    rr, th = np.meshgrid(grid.r_c, grid.th_c, indexing="ij")
    if prof == "uniform":
        f = np.ones_like(rr)
    elif prof == "cosine_theta":
        f = 1.0 + A * np.cos(th - np.radians(cfg.feed.profile_theta0_deg))
    elif prof == "linear_r":
        f = 1.0 + A * (2.0 * rr / grid.R - 1.0)
    elif prof == "central_disc":
        # gas enters only through r <= r_in (fraction of each ring inside r_in)
        r_in = cfg.feed.profile_radius_fraction * grid.R
        lo, hi = grid.r_f[:-1], grid.r_f[1:]
        frac = np.clip((np.minimum(hi, r_in) ** 2 - lo ** 2) / (hi ** 2 - lo ** 2), 0.0, 1.0)
        f = np.repeat(frac[:, None], grid.Nth, axis=1)
    elif prof == "sector":
        # gas enters only through an angular sector centred on profile_theta0_deg
        half = 0.5 * np.radians(cfg.feed.profile_sector_deg)
        t0 = np.radians(cfg.feed.profile_theta0_deg)
        lo = grid.th_f[:-1] - t0
        hi = grid.th_f[1:] - t0
        # wrap-aware overlap of [lo,hi] with [-half, half] (try shifts of 2 pi)
        ov = np.zeros(grid.Nth)
        for sh in (-2 * np.pi, 0.0, 2 * np.pi):
            ov += np.clip(np.minimum(hi + sh, half) - np.maximum(lo + sh, -half), 0.0, None)
        f = np.repeat((ov / grid.dth)[None, :], grid.Nr, axis=0)
    else:
        raise ValueError(prof)
    w = grid.Az * np.ones_like(f)
    return f * w.sum() / (f * w).sum()


def describe(step: Dict, conds: Dict[str, EndCondition]) -> str:
    parts = [f"{v}={step['valves'][v]}" for v in ("V1", "V2", "V3")]
    bc = [f"{e}:{c.kind}" + (f"@{c.P_target:.0f}Pa" if c.kind == "pressure" else "") for e, c in conds.items()]
    return f"[{step['name']}] " + " ".join(parts) + " | " + ", ".join(bc)
