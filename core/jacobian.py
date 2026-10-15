"""dR/dy for the implicit stages.

R(y) = y - gamma h f(t, y) - rhs  =>  dR/dy = I - gamma h J,  J = df/dy.

J is obtained by coloured forward differences on the exact sparsity pattern of
numerics/sparsity.py (one RHS evaluation per colour group, typically ~100-150
for the 3D model regardless of grid size).  Every entry of the pattern is
computed -- no structural approximation -- so modified Newton converges at its
theoretical rate.  An analytical Jacobian is not used: the RHS contains flux
limiters, Ergun inversion and co-adsorption isotherms whose hand-coded
derivatives would be a large, error-prone code path for little gain.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from numerics.sparsity import pattern_and_colors


class FDJacobian:
    def __init__(self, model):
        g = model.grid
        self.model = model
        self.P, self.color = pattern_and_colors(g.Nr, g.Nth, g.Nz, model.ng, model.na, model.n_w > 0)
        assert self.P.shape[0] == model.n_core, (self.P.shape, model.n_core)
        P = self.P.tocsc()
        self.indptr, self.indices = P.indptr, P.indices
        self.n = P.shape[0]
        self.ngroups = int(self.color.max()) + 1
        col_of_entry = np.repeat(np.arange(self.n), np.diff(self.indptr))
        self.entry_col = col_of_entry
        self.entry_row = self.indices
        self.groups = [np.where(self.color == k)[0] for k in range(self.ngroups)]
        ent_group = self.color[col_of_entry]
        order = np.argsort(ent_group, kind="stable")
        bounds = np.searchsorted(ent_group[order], np.arange(self.ngroups + 1))
        self.group_entries = [order[bounds[k]:bounds[k + 1]] for k in range(self.ngroups)]
        # typical magnitudes for the perturbation size
        n = g.ncell
        typ = []
        for s in model.species:
            typ.append(np.full(n, {"CO2": 1e-3, "H2O": 1e-2, "N2": 1.0}[s]))
        typ += [np.full(n, 1e-3)] * model.na
        typ.append(np.full(n, model.rho_b * model.cfg.sorbent.cp_solid))
        if model.n_w:
            typ.append(np.ones(model.n_w))
        self.typ = np.concatenate(typ)
        self.nevals = 0
        # ---- theta-block reduction (cheap preconditioner): drop couplings between
        # different theta sectors.  Exact for axisymmetric states (all sectors equal).
        nv = model.ng + model.na + 1
        idx = np.arange(self.n)
        jsec = np.empty(self.n, dtype=np.int64)
        cellpart = idx < nv * n
        jsec[cellpart] = ((idx[cellpart] % n) // g.Nz) % g.Nth
        jsec[~cellpart] = (idx[~cellpart] - nv * n) // g.Nz
        ring = np.empty(self.n, dtype=np.int64)
        ring[cellpart] = (idx[cellpart] % n) // (g.Nth * g.Nz)
        ring[~cellpart] = g.Nr - 1
        cut = model.cfg.solver.precond_theta_rings - 1
        # keep couplings inside a sector, plus ALL couplings among the innermost rings: the wedge
        # cells around the axis have small arcs r*dtheta and are stiffly coupled in theta
        keep = (jsec[self.entry_row] == jsec[self.entry_col]) | \
               ((ring[self.entry_row] <= cut) & (ring[self.entry_col] <= cut))
        self.keep_theta = np.where(keep)[0]
        kept_cols = self.entry_col[self.keep_theta]
        self.red_indices = self.entry_row[self.keep_theta]
        self.red_indptr = np.concatenate([[0], np.cumsum(np.bincount(kept_cols, minlength=self.n))])

    def __call__(self, t, y, f0=None) -> sp.csc_matrix:
        m = self.model
        if f0 is None:
            f0, _ = m.rhs(t, y, want_acc=False)
        h = 1.5e-8 * np.maximum(np.abs(y), self.typ)
        h = (y + h) - y                          # exactly representable steps
        data = np.empty(len(self.entry_row))
        for k in range(self.ngroups):
            cols = self.groups[k]
            yp = y.copy()
            yp[cols] += h[cols]
            fp, _ = m.rhs(t, yp, want_acc=False)
            ent = self.group_entries[k]
            data[ent] = (fp[self.entry_row[ent]] - f0[self.entry_row[ent]]) / h[self.entry_col[ent]]
        self.nevals += 1
        return sp.csc_matrix((data, self.indices, self.indptr), shape=(self.n, self.n))

    def theta_block(self, J: sp.csc_matrix) -> sp.csc_matrix:
        """J without couplings between different theta sectors (same structure order)."""
        return sp.csc_matrix((J.data[self.keep_theta], self.red_indices, self.red_indptr),
                             shape=(self.n, self.n))
