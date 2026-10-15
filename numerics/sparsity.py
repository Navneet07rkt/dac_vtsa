"""Jacobian sparsity pattern and column colouring.

Which unknowns couple to which follows from the stencils:

* gas-species and energy rows of cell (i,j,k) depend on ALL cell unknowns of
  (i,j,k), (i,j,k+-1), (i,j,k+-2) [MUSCL in z], (i+-1,j,k), (i,j+-1 mod Nth,k);
  outer-ring cells also on the wall cell (j,k);
* loading rows (LDF + isotherm) depend only on their own cell;
* wall rows depend on wall neighbours (j+-1, k+-1) and the outer bed cell.

Columns that never share a row can be perturbed together in a finite-
difference Jacobian (Curtis-Powell-Reid).  :func:`color_columns` finds such
groups greedily on the column-intersection graph P^T P.  The number of groups
(~ n_vars_per_cell x ~20) -- not the number of unknowns -- sets the cost of
one Jacobian, which makes the implicit solver affordable in 3D.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
import scipy.sparse as sp


def cell_index(Nr, Nth, Nz):
    return np.arange(Nr * Nth * Nz).reshape(Nr, Nth, Nz)


def build_pattern(Nr: int, Nth: int, Nz: int, ng: int, na: int, wall: bool) -> sp.csc_matrix:
    n = Nr * Nth * Nz
    nv = ng + na + 1
    idx = cell_index(Nr, Nth, Nz)
    # ---- cell-to-cell stencil pairs (row cell, col cell)
    pairs = [(idx.ravel(), idx.ravel())]
    for dk in (1, 2):
        if Nz > dk:
            pairs.append((idx[:, :, dk:].ravel(), idx[:, :, :-dk].ravel()))
            pairs.append((idx[:, :, :-dk].ravel(), idx[:, :, dk:].ravel()))
    if Nr > 1:
        pairs.append((idx[1:].ravel(), idx[:-1].ravel()))
        pairs.append((idx[:-1].ravel(), idx[1:].ravel()))
    if Nth > 1:
        nb = np.roll(idx, -1, axis=1)
        pairs.append((idx.ravel(), nb.ravel()))
        pairs.append((nb.ravel(), idx.ravel()))
    rc = np.concatenate([p[0] for p in pairs])
    cc = np.concatenate([p[1] for p in pairs])
    rows, cols = [], []
    flux_vars = list(range(ng)) + [ng + na]            # gas species + energy rows
    for vr in flux_vars:
        for vc in range(nv):
            rows.append(vr * n + rc)
            cols.append(vc * n + cc)
    own = idx.ravel()
    for vr in range(ng, ng + na):                       # loading rows: own cell only
        for vc in range(nv):
            rows.append(vr * n + own)
            cols.append(vc * n + own)
    if wall:
        base = nv * n
        widx = np.arange(Nth * Nz).reshape(Nth, Nz)
        outer = idx[-1]                                 # (Nth, Nz) outer-ring cells
        wr, wc = [widx.ravel()], [widx.ravel()]
        if Nz > 1:
            wr += [widx[:, 1:].ravel(), widx[:, :-1].ravel()]
            wc += [widx[:, :-1].ravel(), widx[:, 1:].ravel()]
        if Nth > 1:
            nbw = np.roll(widx, -1, axis=0)
            wr += [widx.ravel(), nbw.ravel()]
            wc += [nbw.ravel(), widx.ravel()]
        rows.append(base + np.concatenate(wr))
        cols.append(base + np.concatenate(wc))
        for vc in range(nv):                            # wall rows <- outer cell
            rows.append(base + widx.ravel())
            cols.append(vc * n + outer.ravel())
        for vr in flux_vars:                            # outer-cell rows <- wall
            rows.append(vr * n + outer.ravel())
            cols.append(base + widx.ravel())
    r = np.concatenate(rows)
    c = np.concatenate(cols)
    N = nv * n + (Nth * Nz if wall else 0)
    P = sp.coo_matrix((np.ones(len(r), dtype=bool), (r, c)), shape=(N, N)).tocsc()
    P.sum_duplicates()
    P.data[:] = True
    return P


def color_columns(P: sp.csc_matrix) -> np.ndarray:
    """Greedy distance-2 colouring: group id per column (no two columns of a group share a row)."""
    P = P.astype(np.int8).tocsc()
    G = (P.T @ P).tocsr()                   # columns i, j conflict iff G[i, j] != 0
    n = P.shape[1]
    color = -np.ones(n, dtype=np.int64)
    indptr, indices = G.indptr, G.indices
    # order: largest degree first tends to use fewer colours
    order = np.argsort(-np.diff(indptr), kind="stable")
    for j in order:
        nb = indices[indptr[j]:indptr[j + 1]]
        used = color[nb]
        used = used[used >= 0]
        if used.size == 0:
            color[j] = 0
            continue
        mark = np.zeros(used.max() + 2, dtype=bool)
        mark[used] = True
        color[j] = int(np.argmin(mark))
    return color


def verify_coloring(P: sp.csc_matrix, color: np.ndarray) -> bool:
    P = P.tocsr()
    for r in range(P.shape[0]):
        cols = P.indices[P.indptr[r]:P.indptr[r + 1]]
        cc = color[cols]
        if len(np.unique(cc)) != len(cc):
            return False
    return True


@lru_cache(maxsize=16)
def pattern_and_colors(Nr, Nth, Nz, ng, na, wall):
    P = build_pattern(Nr, Nth, Nz, ng, na, wall)
    return P, color_columns(P)
