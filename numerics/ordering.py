"""Geometric nested-dissection ordering for the structured (r, theta, z) grid.

SuperLU's built-in COLAMD ordering produces heavy fill for 3D stencils.  For a
structured grid a geometric nested dissection (George 1973) is easy to build
and much better: recursively cut the box with a separator plane, number the
two halves first and the separator last.

Separator thickness follows the stencil: 2 layers along z (MUSCL couples
k +- 2), 1 layer along r and theta.  theta is periodic, so the first cut along
theta removes two planes (j = 0 and j = Nth/2) and turns the ring into two
open segments.

The ordering is at the cell level; all unknowns of a cell (gas species,
loadings, energy) -- and, for outer-ring cells, the adjacent wall unknown --
are numbered consecutively.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np


def _nd_cells(Nr, Nth, Nz, leaf=16):
    order = []

    def rec(r0, r1, t_idx, z0, z1):
        # t_idx: list of theta indices (open segment, no wrap)
        nr, nt, nz = r1 - r0, len(t_idx), z1 - z0
        if nr * nt * nz == 0:
            return
        if nr * nt * nz <= leaf:
            for k in range(z0, z1):
                for j in t_idx:
                    for i in range(r0, r1):
                        order.append((i, j, k))
            return
        # extents in "separator cost" units: z separators are 2 layers thick
        cand = []
        if nz >= 5:
            cand.append((nz / 2.0, "z"))
        if nr >= 3:
            cand.append((nr, "r"))
        if nt >= 3:
            cand.append((nt, "t"))
        if not cand:
            for k in range(z0, z1):
                for j in t_idx:
                    for i in range(r0, r1):
                        order.append((i, j, k))
            return
        _, d = max(cand)
        if d == "z":
            ks = z0 + (nz - 2) // 2
            rec(r0, r1, t_idx, z0, ks)
            rec(r0, r1, t_idx, ks + 2, z1)
            for k in (ks, ks + 1):
                for j in t_idx:
                    for i in range(r0, r1):
                        order.append((i, j, k))
        elif d == "r":
            rs = r0 + nr // 2
            rec(r0, rs, t_idx, z0, z1)
            rec(rs + 1, r1, t_idx, z0, z1)
            for k in range(z0, z1):
                for j in t_idx:
                    order.append((rs, j, k))
        else:
            m = nt // 2
            rec(r0, r1, t_idx[:m], z0, z1)
            rec(r0, r1, t_idx[m + 1:], z0, z1)
            for k in range(z0, z1):
                for i in range(r0, r1):
                    order.append((i, t_idx[m], k))

    if Nth >= 4:
        # periodic ring: cut at j = 0 and j = Nth//2, giving two open segments
        half = Nth // 2
        seg1 = list(range(1, half))
        seg2 = list(range(half + 1, Nth))
        rec(0, Nr, seg1, 0, Nz)
        rec(0, Nr, seg2, 0, Nz)
        for jsep in (0, half):
            for k in range(Nz):
                for i in range(Nr):
                    order.append((i, jsep, k))
    else:
        rec(0, Nr, list(range(Nth)), 0, Nz)
    assert len(order) == Nr * Nth * Nz and len(set(order)) == len(order)
    return order


@lru_cache(maxsize=16)
def nested_dissection(Nr, Nth, Nz, nv, wall):
    """Permutation p (new -> old index) of the core unknown vector (variable-major layout)."""
    n = Nr * Nth * Nz
    cells = _nd_cells(Nr, Nth, Nz)
    perm = []
    for (i, j, k) in cells:
        cid = (i * Nth + j) * Nz + k
        perm.extend(v * n + cid for v in range(nv))
        if wall and i == Nr - 1:
            perm.append(nv * n + j * Nz + k)
    p = np.array(perm, dtype=np.int64)
    assert len(np.unique(p)) == len(p)
    return p
