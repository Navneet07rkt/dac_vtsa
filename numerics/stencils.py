"""Finite-difference / finite-volume stencil utilities.

* :func:`fd_weights` -- Fornberg (1988) algorithm: weights of any derivative
  order on any (non-uniform) set of nodes.  Used by the tests and available for
  higher-order post-processing (e.g. gradients of exported fields).
* TVD reconstruction (:func:`limited_slope`, :func:`face_values_z`) used by the
  axial convective fluxes (MUSCL with a choice of limiter).
* :func:`upwind` -- first-order upwind face value (radial/angular convection).
"""
from __future__ import annotations

import numpy as np


def fd_weights(x0: float, x: np.ndarray, m: int) -> np.ndarray:
    """Fornberg weights c[j] so that f^(m)(x0) ~= sum_j c[j] f(x[j])."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    c = np.zeros((n, m + 1))
    c1, c4 = 1.0, x[0] - x0
    c[0, 0] = 1.0
    for i in range(1, n):
        mn = min(i, m)
        c2, c5, c4 = 1.0, c4, x[i] - x0
        for j in range(i):
            c3 = x[i] - x[j]
            c2 *= c3
            if j == i - 1:
                for k in range(mn, 0, -1):
                    c[i, k] = c1 * (k * c[i - 1, k - 1] - c5 * c[i - 1, k]) / c2
                c[i, 0] = -c1 * c5 * c[i - 1, 0] / c2
            for k in range(mn, 0, -1):
                c[j, k] = (c4 * c[j, k] - k * c[j, k - 1]) / c3
            c[j, 0] = c4 * c[j, 0] / c3
        c1 = c2
    return c[:, m]


# Width of the smooth upwind switch [m/s].  Below |u| ~ U_SMOOTH the face value blends to the
# central average.  That is accurate AND stable there: upwinding is only needed when the cell
# Peclet number v dx / D exceeds ~2, i.e. above ~1e-2 m/s for this model's dispersion
# (D ~ 1e-5 m2/s) and cell sizes (0.3-3 mm); at 1e-4 m/s the cell Peclet number is < 0.1.
# A much sharper switch (1e-6) made near-stagnant flows (end of pressurisation, closed ends)
# chatter between upwind directions from one time step to the next.
U_SMOOTH = 1e-4


def upwind_weight(u, us=U_SMOOTH):
    """Smooth upwind weight w(u) = 0.5 (1 + tanh(u/us)): 1 for u >> us, 0 for u << -us.

    A hard switch np.where(u >= 0, ...) has a jump in d(flux)/d(state) at u = 0; Newton and
    the finite-difference Jacobian-vector products used by Newton-Krylov then see
    inconsistent derivatives wherever a velocity changes sign (pressure waves, stagnant
    zones).  The blend changes the flux by < |u| |dphi| for |u| < us -- central differencing
    there, which is stable because dispersion dominates at such small velocities."""
    return 0.5 * (1.0 + np.tanh(u / us))


def limited_slope(dm: np.ndarray, dp: np.ndarray, limiter: str) -> np.ndarray:
    """Limited slope (per cell) from backward (dm) and forward (dp) differences."""
    if limiter == "upwind":
        return np.zeros_like(dm)
    prod = dm * dp
    if limiter == "van_albada":
        # original (unclipped) van Albada: s = dm dp (dm + dp)/(dm^2 + dp^2): smooth everywhere
        # except the origin, TVD-like away from extrema -- Newton friendly
        den = dm * dm + dp * dp + 1e-300
        return prod * (dm + dp) / den
    if limiter == "van_albada_tvd":
        den = dm * dm + dp * dp + 1e-300
        return np.where(prod > 0.0, prod * (dm + dp) / den, 0.0)
    if limiter == "van_leer":
        return np.where(prod > 0.0, 2.0 * prod / (dm + dp + 1e-300), 0.0)
    if limiter == "minmod":
        return np.where(prod > 0.0, np.sign(dm) * np.minimum(np.abs(dm), np.abs(dp)), 0.0)
    raise ValueError(limiter)


def face_values_z(phi: np.ndarray, u: np.ndarray, limiter: str) -> np.ndarray:
    """MUSCL face values on interior z-faces (last axis).

    phi : (..., Nz) cell values;   u : (..., Nz-1) face velocity (+ = +z)
    returns (..., Nz-1) upwind-biased, limited face values.  Cells adjacent to a
    boundary use zero slope (first order), so no ghost cells are needed.
    """
    nz = phi.shape[-1]
    w = upwind_weight(u)
    if nz < 3 or limiter == "upwind":
        return w * phi[..., :-1] + (1.0 - w) * phi[..., 1:]
    d = np.diff(phi, axis=-1)                      # (..., Nz-1)
    s = np.zeros_like(phi)
    s[..., 1:-1] = limited_slope(d[..., :-1], d[..., 1:], limiter)
    left = phi[..., :-1] + 0.5 * s[..., :-1]       # value at face from the cell on its left
    right = phi[..., 1:] - 0.5 * s[..., 1:]        # value at face from the cell on its right
    return w * left + (1.0 - w) * right


def upwind(phi_minus: np.ndarray, phi_plus: np.ndarray, u: np.ndarray) -> np.ndarray:
    """First-order (smoothly switched) upwind face value between phi_minus and phi_plus."""
    w = upwind_weight(u)
    return w * phi_minus + (1.0 - w) * phi_plus
