"""Cylindrical (r, theta, z) finite-volume grid.

Cell (i, j, k) occupies r in [r_f[i], r_f[i+1]], theta in [th_f[j], th_f[j+1]],
z in [z_f[k], z_f[k+1]].  z = 0 is the TOP of the column (feed end, V1),
z = L the BOTTOM (V2).  All geometric factors needed by a conservative
finite-volume discretisation are precomputed here:

* V        cell volumes                      (Nr, 1, 1)  -> broadcast
* Az       axial face areas                  (Nr, 1)     (same for every z face)
* Ar       radial face areas at r_f[1..Nr-1] (Nr-1, 1, 1)
* Ath      angular face areas                (Nr, 1, 1)
* dzc, drc, arc   centre-to-centre distances used in gradients

The axis r = 0 is handled by the finite-volume construction itself: the inner
face of the innermost ring has zero area, so no flux crosses it and no 1/r
singularity ever appears (the 'negative axis temperature' problem of a
finite-difference r-discretisation cannot occur).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Grid:
    R: float
    L: float
    Nr: int
    Nth: int
    Nz: int
    wall_thickness: float = 0.002
    radial_spacing: str = "uniform"
    # filled in __post_init__
    r_f: np.ndarray = field(init=False)
    r_c: np.ndarray = field(init=False)
    th_f: np.ndarray = field(init=False)
    th_c: np.ndarray = field(init=False)
    z_f: np.ndarray = field(init=False)
    z_c: np.ndarray = field(init=False)

    def __post_init__(self):
        Nr, Nth, Nz = self.Nr, self.Nth, self.Nz
        if self.radial_spacing == "equal_area":
            self.r_f = self.R * np.sqrt(np.linspace(0.0, 1.0, Nr + 1))
        else:
            self.r_f = np.linspace(0.0, self.R, Nr + 1)
        self.r_c = 0.5 * (self.r_f[1:] + self.r_f[:-1])
        self.dr = np.diff(self.r_f)
        self.dth = 2.0 * np.pi / Nth
        self.th_f = np.linspace(0.0, 2.0 * np.pi, Nth + 1)
        self.th_c = 0.5 * (self.th_f[1:] + self.th_f[:-1])
        self.dz = self.L / Nz
        self.z_f = np.linspace(0.0, self.L, Nz + 1)
        self.z_c = 0.5 * (self.z_f[1:] + self.z_f[:-1])

        ring = 0.5 * (self.r_f[1:] ** 2 - self.r_f[:-1] ** 2)            # (Nr,)
        self.Az = (ring * self.dth)[:, None]                              # (Nr,1)
        self.V = (ring * self.dth * self.dz)[:, None, None]               # (Nr,1,1)
        self.Ar = (self.r_f[1:-1] * self.dth * self.dz)[:, None, None]    # (Nr-1,1,1)
        self.Ath = (self.dr * self.dz)[:, None, None]                     # (Nr,1,1)
        self.drc = np.diff(self.r_c)[:, None, None]                       # (Nr-1,1,1)
        self.arc = (self.r_c * self.dth)[:, None, None]                   # (Nr,1,1)
        self.dzc = self.dz
        # outer bed face (r = R) and the wall ring
        self.A_wall_in = self.R * self.dth * self.dz                      # per (theta, z) wall cell
        Ro = self.R + self.wall_thickness
        self.A_wall_out = Ro * self.dth * self.dz
        self.V_wall = 0.5 * (Ro ** 2 - self.R ** 2) * self.dth * self.dz
        self.A_wall_axial = 0.5 * (Ro ** 2 - self.R ** 2) * self.dth     # wall cross-section for axial conduction
        self.R_wall_mid = self.R + 0.5 * self.wall_thickness
        self.half_dr_outer = self.R - self.r_c[-1]
        self.cross_section = np.pi * self.R ** 2
        self.bed_volume = self.cross_section * self.L
        self.ncell = Nr * Nth * Nz
        self.shape = (Nr, Nth, Nz)
        self.wall_shape = (Nth, Nz)

    # ------------------------------------------------------------------ helpers
    @classmethod
    def from_config(cls, cfg) -> "Grid":
        return cls(R=cfg.geometry.bed_radius, L=cfg.geometry.bed_length, Nr=cfg.grid.Nr,
                   Nth=cfg.grid.Ntheta, Nz=cfg.grid.Nz, wall_thickness=cfg.geometry.wall_thickness,
                   radial_spacing=cfg.grid.radial_spacing)

    def cell_centres(self):
        """Broadcast arrays (r, theta, z) of cell centres, each (Nr, Nth, Nz)."""
        return np.meshgrid(self.r_c, self.th_c, self.z_c, indexing="ij")

    def volume_full(self) -> np.ndarray:
        return np.broadcast_to(self.V, self.shape)

    def integrate(self, field3d: np.ndarray) -> float:
        """Volume integral of a cell field (Nr, Nth, Nz) -> scalar."""
        return float(np.sum(field3d * self.V))

    def coverage_fraction(self, theta_deg, z_frac) -> np.ndarray:
        """Fraction of each wall cell (Nth, Nz) covered by an angular x axial window."""
        t0, t1 = np.radians(theta_deg[0]), np.radians(theta_deg[1])
        z0, z1 = z_frac[0] * self.L, z_frac[1] * self.L
        ov_t = np.clip(np.minimum(self.th_f[1:], t1) - np.maximum(self.th_f[:-1], t0), 0, None) / self.dth
        ov_z = np.clip(np.minimum(self.z_f[1:], z1) - np.maximum(self.z_f[:-1], z0), 0, None) / self.dz
        return ov_t[:, None] * ov_z[None, :]

    def describe(self) -> str:
        return (f"Grid: Nr={self.Nr} Ntheta={self.Nth} Nz={self.Nz} ({self.ncell} cells), "
                f"R={self.R} m, L={self.L} m, dr={self.dr.mean():.3e} m, dz={self.dz:.3e} m")
