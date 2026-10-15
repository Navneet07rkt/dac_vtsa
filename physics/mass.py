"""Gas-phase species mass balance (finite-volume flux assembly).

Conservative form for every gas species k in every cell (volume V):

    d(eps_t c_k V)/dt = - sum_faces (N_k . n A) - rho_b V dq_k/dt

    N_k = u c_k  - eps_b D C grad(y_k)      (convection + dispersion)

* convection: MUSCL/TVD along z (limiter configurable), first-order upwind in
  r and theta (transverse velocities are small);
* dispersion acts on mole-fraction gradients with ONE coefficient for all
  species, so sum_k J_k = 0 exactly: dispersion never creates net molar flow;
* the divergence telescopes, so the sum over cells equals the boundary flux
  -> species inventories are conserved to round-off (checked in validation/).
"""
from __future__ import annotations

import numpy as np

from numerics.stencils import face_values_z, upwind


def theta_next(a):
    """Periodic neighbour j+1 along the theta axis (axis -2); faster than np.roll."""
    n = a.shape[-2]
    return np.take(a, np.r_[1:n, 0], axis=-2)


def theta_prev(a):
    n = a.shape[-2]
    return np.take(a, np.r_[n - 1, 0:n - 1], axis=-2)


def convective_z(phi, uz, Az, limiter):
    """phi (n, Nr, Nth, Nz), uz (Nr, Nth, Nz-1) -> flux (n, Nr, Nth, Nz-1) in +z [phi*m3/s]."""
    return uz * Az[:, :, None] * face_values_z(phi, uz, limiter)


def convective_r(phi, ur, Ar):
    """phi (n, Nr, ...), ur (Nr-1, Nth, Nz) -> flux through interior r faces (+ = outward)."""
    return ur * Ar * upwind(phi[:, :-1], phi[:, 1:], ur)


def convective_theta(phi, uth, Ath):
    """Flux through the face between theta-cells j and j+1 (periodic), + = increasing theta."""
    return uth * Ath * upwind(phi, theta_next(phi), uth)


def dispersive_z(y, Cf, D, eps_b, Az, dz):
    return -eps_b * D * Cf * np.diff(y, axis=-1) / dz * Az[:, :, None]


def dispersive_r(y, Cf, D, eps_b, Ar, drc):
    return -eps_b * D * Cf * np.diff(y, axis=1) / drc * Ar


def dispersive_theta(y, Cf, D, eps_b, Ath, arc):
    return -eps_b * D * Cf * (theta_next(y) - y) / arc * Ath


def divergence(Fz_full, Fr_int, Fth, nth):
    """Net outflow per cell.

    Fz_full : (..., Nr, Nth, Nz+1) fluxes in +z incl. boundary faces
    Fr_int  : (..., Nr-1, Nth, Nz) interior radial fluxes (r=0 and r=R faces carry no mass)
    Fth     : (..., Nr, Nth, Nz) fluxes through face j+1/2, or None if Nth == 1
    """
    out = Fz_full[..., 1:] - Fz_full[..., :-1]
    if Fr_int is not None and Fr_int.shape[-3] > 0:
        out[..., :-1, :, :] += Fr_int
        out[..., 1:, :, :] -= Fr_int
    if Fth is not None and nth > 1:
        out += Fth - theta_prev(Fth)
    return out
