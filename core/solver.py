"""Implicit Newton-Raphson time integration (production path).

Scheme: TR-BDF2 (Bank et al. 1985; Hosea & Shampine 1996) -- one-step,
L-stable, second order; two implicit stages share the iteration matrix
M = I - d h J (d = gamma/2, gamma = 2 - sqrt 2):

    stage 1 (trapezoid to t + gamma h):  z1 - d h f(z1) = y_n + d h f_n
    stage 2 (BDF2 to t + h):             z2 - d h f(z2) = c1 z1 - c2 y_n

Each stage is a nonlinear system R(z) = 0 (core/residual.py) solved by Newton:

1. *modified Newton* with the stored sparse LU of M (cheap: one RHS + one
   triangular solve per iteration);
2. if that stalls (moving fronts make the stored Jacobian stale, or h has
   drifted from the factorised h), *Newton-Krylov*: the exact Newton
   correction is obtained with GMRES on the true operator (I - d h J(z)) v,
   whose action is a finite-difference directional derivative, preconditioned
   by the stored LU.  This converges where modified Newton does not, without
   a new Jacobian or factorisation;
3. only if both fail is J re-evaluated (coloured FD, core/jacobian.py), M
   re-factorised (SuperLU), and finally h reduced.

Local error: difference to a 3rd-order quadrature of the three stage
derivatives, filtered by M^-1 (Shampine) so that stiff components do not
trigger spurious rejections; h_new = h * clip(0.9 err^-1/3, 0.2, 4).

L-stability matters: the pressure field relaxes in ~1e-7 s while the CO2
front moves over ~1e3 s.  Explicit RK4 would need ~1e10 steps for one
adsorption step (numerics/cfl.py); TR-BDF2 needs a few hundred.

Accumulators (boundary-flow and heat integrals) are advanced with the same
stage formulas after each converged stage -> mass and energy balances close
to the Newton tolerance.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Callable, List, Optional

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .jacobian import FDJacobian

GAMMA = 2.0 - np.sqrt(2.0)
D = GAMMA / 2.0
C1 = 1.0 / (GAMMA * (2.0 - GAMMA))
C2 = (1.0 - GAMMA) ** 2 / (GAMMA * (2.0 - GAMMA))
W1 = (3.0 * GAMMA - 1.0) / (6.0 * GAMMA)
W2 = 1.0 / (6.0 * GAMMA * (1.0 - GAMMA))
W3 = (2.0 - 3.0 * GAMMA) / (6.0 * (1.0 - GAMMA))
SQRT_EPS = np.sqrt(np.finfo(float).eps)
# Alexander's L-stable, stiffly accurate SDIRK2 (both stages damp stiff modes)
G_SD = 1.0 - 1.0 / np.sqrt(2.0)
SW1 = (3.0 * G_SD - 1.0) / (6.0 * G_SD)
SW2 = 1.0 / (6.0 * G_SD * (1.0 - G_SD))
SW3 = (2.0 - 3.0 * G_SD) / (6.0 * (1.0 - G_SD))


class SolverError(RuntimeError):
    pass


@dataclass
class IntegrationResult:
    y: np.ndarray
    acc: np.ndarray
    t: float
    n_steps: int = 0
    n_reject: int = 0
    n_newton_fail: int = 0
    n_krylov: int = 0
    n_jac: int = 0
    n_lu: int = 0
    nfev: int = 0
    wall_time: float = 0.0
    snapshots: List = field(default_factory=list)       # (t, y, acc)
    stopped_early: bool = False
    stop_value: Optional[float] = None


def _rms(x):
    return float(np.sqrt(np.mean(x * x)))


def _cons_ok(cons, R, s, tol):
    """Global balance test of a stage solution: for every conserved quantity m (moles of each gas
    species incl. its adsorbed phase, bed + wall energy)  |W_m R| <= tol * W_m s.  Internal fluxes
    telescope, so W_m R is exactly the stage's mismatch between inventory change and boundary flows.
    A step-size convergence test alone is not enough for this: in stiff directions (pressure) a
    state error far below tolerance still leaves a large residual, i.e. a balance error."""
    if cons is None or tol <= 0.0:
        return True
    return bool(np.all(np.abs(cons @ R) <= tol * (cons @ s)))


def _truncate(z_new, z_old, pos, floor=0.1):
    """Keep concentrations / loadings positive along the Newton path: a component may drop by at
    most (1 - floor) of its current value per iteration.  Iterates then never enter the
    non-physical, non-smooth region c < 0 (where isotherms, upwinding and FD Jacobian-vector
    products misbehave); the converged solution is unchanged because it is positive."""
    if pos is None:
        return z_new
    lo = floor * np.maximum(z_old[pos], 0.0)
    z_new[pos] = np.maximum(z_new[pos], lo)
    return z_new


def newton_modified(model, t, z0, rhs_known, dh, lu, scale_fn, max_iter, kappa, fail_rate=0.9, pos=None,
                    cons=None, cons_tol=0.0):
    """Modified Newton with a fixed LU.  Returns (z, ok, iters, rate)."""
    z = z0.copy()
    nrm_prev = None
    rate = 0.0
    done = False
    for it in range(max_iter + 2):
        f, _ = model.rhs(t, z, want_acc=False)
        R = z - dh * f - rhs_known
        if not np.all(np.isfinite(R)):
            return z0, False, it, np.inf
        if done:                                   # converged in the step norm: check the balances
            if _cons_ok(cons, R, scale_fn(z), cons_tol) or it >= max_iter + 1:
                return z, True, it, rate
        elif it >= max_iter:
            break
        dz = lu.solve(-R)
        z_old = z
        z = _truncate(z + dz, z_old, pos)
        nrm = _rms((z - z_old) / scale_fn(z))
        if not np.isfinite(nrm):
            return z0, False, it + 1, np.inf
        if done:
            continue
        # no acceptance on the first iteration: a stale LU can produce a small update while
        # the residual is still large (false convergence).  The rate estimate needs >= 2 iterates.
        if nrm_prev is not None:
            rate = nrm / max(nrm_prev, 1e-300)
            if rate >= fail_rate:
                return z, False, it + 1, rate
            if rate / (1.0 - rate) * nrm <= kappa:
                if cons is None or cons_tol <= 0.0:
                    return z, True, it + 1, rate
                done = True
        elif nrm < 1e-8:
            if cons is None or cons_tol <= 0.0:
                return z, True, it + 1, 0.0
            done = True
        nrm_prev = nrm
    return z, False, max_iter, rate


def newton_krylov(model, t, z0, rhs_known, dh, lu, scale_fn, max_iter, kappa, krylov_dim=30,
                  eta=1e-2, counter=None, pos=None, cons=None, cons_tol=0.0, stats=None, trace=None,
                  refresh=None):
    """Inexact Newton with backtracking line search.

    Newton direction: GMRES on (I - dh J(z)) dz = -R, with J v by finite differences
    (true Jacobian action) preconditioned by the stored LU, in tolerance-scaled variables.
    Globalisation: Armijo backtracking on ||R/s|| (lambda = 1, 1/2, 1/4, 1/8) and positivity
    truncation of concentrations/loadings.  `counter` accumulates GMRES iterations.

    refresh(z) -> new LU (optional): called at most once per solve, when GMRES cannot solve the
    linear system with the current preconditioner (fast transients make a stored Jacobian a poor
    preconditioner); the caller keeps the new factorisation for the following stages/steps.

    Returns (z, ok, iterations, rate, new_lu_or_None)."""
    z = z0.copy()
    f, _ = model.rhs(t, z, want_acc=False)
    s = scale_fn(z)
    R = z - dh * f - rhs_known
    rn = _rms(R / s)
    if not np.isfinite(rn):
        return z0, False, 0, np.inf, None
    # Progress is measured AFFINE-INVARIANTLY (Deuflhard): by the simplified Newton correction
    # ||M^-1 R / s|| rather than by ||R / s||.  In this model a trace species (CO2 at 400 ppm,
    # scale ~1e-6 mol/m3) is carried by the stiff, pressure-driven bulk flow: scaled Jacobian
    # entries reach ~1e7, so a 1e-3 relative error in a finite-difference derivative shows up as
    # a residual of ~1e4 that corresponds to a state error far below tolerance.  The raw residual
    # norm then rejects good Newton steps (line search stalls); the correction norm does not.
    pnorm = lambda RR, ss: _rms(lu.solve(RR) / ss)          # noqa: E731
    pn = pnorm(R, s)
    n = len(z)
    step_prev = None
    polish = 0                 # extra accurate iterations spent closing the stage balances
    best = None                # (rn, z) of the best state that already passed the convergence test
    new_lu = None
    it = -1
    while it + 1 < max_iter + polish:
        it += 1
        s = scale_fn(z)
        b = -R / s
        zs_norm = max(1.0, _rms(z / s))

        def matvec(u, z=z, f=f, s=s, zs_norm=zs_norm):
            un = _rms(u)
            if un == 0.0:
                return np.zeros_like(u)
            sig = SQRT_EPS * zs_norm / un
            v = u * s
            fp, _ = model.rhs(t, z + sig * v, want_acc=False)
            return (v - dh * (fp - f) / sig) / s

        # RIGHT preconditioning: GMRES on (A P^-1) w = b, x = P^-1 w.  Unlike scipy's built-in
        # left preconditioning, the minimised residual is then the TRUE linear residual, so an
        # inaccurate (stale) preconditioner cannot make GMRES stop early with a useless step.
        def prec(r, s=s, lu=lu):
            return lu.solve(r * s) / s

        AP = spla.LinearOperator((n, n), matvec=lambda w, prec=prec: matvec(prec(w)), dtype=float)
        cb = None
        if polish and stats is not None:            # polish iterations do not signal a stale Jacobian
            def cb(_x):
                stats["gmres_polish"] = stats.get("gmres_polish", 0) + 1
        elif counter is not None:
            def cb(_x):
                counter[0] += 1
        # absolute floor: Newton declares convergence at RMS(R/s) <= kappa, so there is no point
        # in solving the linear system below ~0.1 kappa -- near the round-off floor of the stiff
        # fluxes GMRES would otherwise burn its whole Krylov space chasing noise.
        floor = (1e-4 if polish else 0.1 * kappa) * np.sqrt(n)
        w, info = spla.gmres(AP, b, rtol=(1e-3 if polish else eta), atol=floor, restart=krylov_dim,
                             maxiter=1, callback=cb, callback_type="pr_norm")
        if info != 0 and refresh is not None and new_lu is None and not polish:
            # the stored Jacobian is a poor preconditioner here: re-linearise at the current iterate
            lu = new_lu = refresh(z)
            pnorm = lambda RR, ss, lu=lu: _rms(lu.solve(RR) / ss)   # noqa: E731
            pn = pnorm(R, s)
            if stats is not None:
                stats["jac_refresh_in_newton"] = stats.get("jac_refresh_in_newton", 0) + 1

            def prec(r, s=s, lu=lu):
                return lu.solve(r * s) / s

            AP = spla.LinearOperator((n, n), matvec=lambda w, prec=prec: matvec(prec(w)), dtype=float)
            w, info = spla.gmres(AP, b, rtol=eta, atol=floor, restart=krylov_dim,
                                 maxiter=1, callback=cb, callback_type="pr_norm")
        u = prec(w)
        dz = u * s
        if not np.all(np.isfinite(dz)):
            return z0, False, it + 1, np.inf, new_lu
        lam = 1.0
        for ls in range(4):
            z_try = _truncate(z + lam * dz, z, pos)
            f_try, _ = model.rhs(t, z_try, want_acc=False)
            R_try = z_try - dh * f_try - rhs_known
            rn_try = _rms(R_try / s)
            pn_try = pnorm(R_try, s) if np.isfinite(rn_try) else np.inf
            if np.isfinite(pn_try) and (pn_try <= (1.0 - 1e-4 * lam) * pn or ls == 3):
                break
            lam *= 0.5
        if not np.isfinite(pn_try):
            return z0, False, it + 1, np.inf, new_lu
        step = _rms((z_try - z) / s)
        if trace is not None:
            trace.append({"it": it, "polish": polish, "gmres_info": info, "lam": lam, "step": step,
                          "res_before": rn, "res_after": rn_try, "pres_before": pn, "pres_after": pn_try})
        z, f, R, rn, pn = z_try, f_try, R_try, rn_try, pn_try
        rate = 0.0 if step_prev is None else step / max(step_prev, 1e-300)
        # converged if (a) the last full Newton step is below tolerance, (b) the contraction
        # estimate (only meaningful for accurate linear solves) says the remaining error is, or
        # (c) the simplified Newton correction ||M^-1 R / s|| is -- the affine-invariant measure
        # of the remaining state error, which stays meaningful where the raw residual is
        # dominated by finite-difference noise in stiff couplings.
        converged = ((lam == 1.0 and step <= kappa and (info == 0 or step_prev is not None))
                     or (step_prev is not None and lam == 1.0 and info == 0 and rate < 1.0
                         and rate / (1.0 - rate) * step <= kappa)
                     or pn <= kappa)
        if converged:
            if best is None or rn < best[0]:
                best = (rn, z.copy())
            if _cons_ok(cons, R, scale_fn(z), cons_tol):
                return z, True, it + 1, rate, new_lu
            if polish < 2:
                polish += 1        # state converged, balances not yet: one more accurate iteration
                step_prev = None
                continue
            if stats is not None:
                stats["balance_unclosed"] = stats.get("balance_unclosed", 0) + 1
            return best[1], True, it + 1, rate, new_lu
        if best is not None:       # a polish iteration made things worse: keep the converged state
            if stats is not None:
                stats["balance_unclosed"] = stats.get("balance_unclosed", 0) + 1
            return best[1], True, it + 1, rate, new_lu
        if step_prev is not None and lam == 1.0 and info == 0:
            # early exit (Hairer & Wanner, RADAU5): if the observed contraction cannot reach the
            # tolerance within the remaining iterations, stop now -- the caller then refreshes the
            # Jacobian or cuts h, instead of spending ~30 GMRES iterations per useless Newton step.
            left = max_iter - it - 1
            if it >= 2 and (rate >= 1.0 or rate ** left / (1.0 - rate) * step > kappa * 1e3):
                return z, False, it + 1, rate, new_lu
        step_prev = step
    return z, False, it + 1, 1.0, new_lu


class TRBDF2:
    def __init__(self, model, cfg, h_max: float, h_init: Optional[float] = None, logger=None):
        self.m = model
        s = cfg.solver
        self.rtol = s.rtol
        self.atol = model.atol
        self.h_max = h_max
        self.h_init = h_init if h_init is not None else s.h_initial
        self.max_iter = s.max_newton_iter
        self.kappa = s.newton_kappa
        self.jac_max_age = s.jacobian_max_age
        self.max_steps = s.max_steps_per_step
        self.jac = FDJacobian(model)
        self.log = logger or (lambda *a, **k: None)
        self.I = sp.identity(model.n_core, format="csc")
        self.stats = {"modified_fail": 0, "krylov_ok": 0, "krylov_fail": 0}
        self._gm = [0]
        self.pos = np.arange(model.off_e)          # gas concentrations + loadings stay >= 0
        self.scheme = s.scheme
        self.mode = s.newton_mode
        self.gmres_refresh = s.gmres_refresh     # GMRES iterations per step that trigger a Jacobian refresh
        self.log_refactor = np.log(s.refactor_ratio)
        self.precond = s.preconditioner if model.grid.Nth > 1 else "full"
        self.cons = model.conservation_matrix() if hasattr(model, "conservation_matrix") else None
        self.cons_tol = float(getattr(s, "conservation_tol", 0.0) or 0.0)
        self.newton_refresh = True              # re-linearise inside a stage when GMRES stalls
        self._new_lu = None
        self._new_J = None
        self._fresh = False
        self.cap_after_fail = 0.5               # h limit after a Newton failure with a fresh Jacobian
        self.cap_relax = 1.15                   # ... relaxed by this factor per accepted step
        if self.precond == "theta_block" and self.mode == "modified":
            self.mode = "hybrid"     # an approximate LU needs the Krylov fallback

    def scale(self, y):
        return self.atol + self.rtol * np.abs(y)

    def _factor(self, J, h):
        Jp = self.jac.theta_block(J) if self.precond == "theta_block" else J
        M = (self.I - (D * h) * Jp).tocsc()
        return spla.splu(M, permc_spec="COLAMD")

    def _solve_stage(self, t, z0, rhs_known, dh, lu, res):
        if self.mode in ("hybrid", "modified"):
            z, ok, it, rate = newton_modified(self.m, t, z0, rhs_known, dh, lu, self.scale,
                                              self.max_iter if self.mode == "modified" else 4,
                                              self.kappa, fail_rate=0.9 if self.mode == "modified" else 0.5,
                                              pos=self.pos, cons=self.cons, cons_tol=self.cons_tol)
            if ok or self.mode == "modified":
                return z, ok, rate
            self.stats["modified_fail"] += 1
        h = dh / D

        def refresh(zc):
            fc, _ = self.m.rhs(t, zc, want_acc=False)
            Jn = self.jac(t, zc, fc)
            res.n_jac += 1
            res.n_lu += 1
            self._new_J = (Jn, h)
            return self._factor(Jn, h)

        self._new_J = None
        z, ok, it, rate, new_lu = newton_krylov(self.m, t, z0, rhs_known, dh, lu, self.scale,
                                                self.max_iter, self.kappa, counter=self._gm, pos=self.pos,
                                                cons=self.cons, cons_tol=self.cons_tol, stats=self.stats,
                                                refresh=refresh if (self.newton_refresh and not self._fresh) else None)
        if new_lu is not None:
            self._new_lu = new_lu
        res.n_krylov += 1
        self._krylov_used = True
        self.stats["krylov_ok" if ok else "krylov_fail"] += 1
        return z, ok, (0.0 if ok else rate)

    def integrate(self, y0: np.ndarray, t_end: float, t_out: np.ndarray,
                  on_step: Optional[Callable] = None, stop_fn: Optional[Callable] = None,
                  stop_value: Optional[float] = None, stop_direction: int = 1) -> IntegrationResult:
        m = self.m
        wall0 = time.time()
        t = 0.0
        y = y0.copy()
        acc = np.zeros(m.n_acc)
        f_n, g_n = m.rhs(t, y)
        res = IntegrationResult(y=y, acc=acc, t=t)
        t_out = np.asarray(sorted(set(np.clip(t_out, 0.0, t_end)) | {t_end}))
        k_out = 0
        if t_out[0] <= 0.0:
            res.snapshots.append((0.0, y.copy(), acc.copy()))
            k_out = 1
        h = min(self.h_init, self.h_max, t_end)
        J = None
        lu = None
        h_lu = None
        jac_age = 0
        fresh = False
        slope_prev = None
        h_cap = np.inf          # Newton-convergence limit on h (set by failures, relaxed on success)
        self._gm = [0]
        while t < t_end * (1.0 - 1e-12):
            self._krylov_used = False
            gm0 = self._gm[0]
            if res.n_steps + res.n_reject + res.n_newton_fail > self.max_steps:
                raise SolverError(f"too many steps in '{m.step['name']}' at t={t:.4g}s")
            t_target = t_out[k_out] if k_out < len(t_out) else t_end
            h_desired = h
            h = min(h, self.h_max, t_target - t)
            truncated = h < h_desired * (1.0 - 1e-9)
            if h < 1e-10 * max(1.0, t_end):
                raise SolverError(f"step size underflow (h={h:.3e}) in '{m.step['name']}' at t={t:.6g}s")
            if J is None or jac_age >= self.jac_max_age:
                J = self.jac(t, y, f_n)
                res.n_jac += 1
                jac_age = 0
                fresh = True
                lu = None
            if lu is None or abs(np.log(h / h_lu)) > (self.log_refactor if self.mode != "modified" else 0.35):
                lu = self._factor(J, h)
                h_lu = h
                res.n_lu += 1
            dh = D * h                                   # = G_SD * h as well (same diagonal)
            sd = self.scheme == "sdirk2"
            c_mid = G_SD if sd else GAMMA
            # ---- stage 1 (TR for TR-BDF2, backward Euler to t + gamma h for SDIRK2)
            rhs1 = y if sd else y + dh * f_n
            z_pred = y if slope_prev is None else _truncate(y + c_mid * h * slope_prev, y, self.pos)
            self._new_lu = None
            self._fresh = fresh                     # J evaluated at y_n of this step: no in-stage refresh
            z1, ok, r1 = self._solve_stage(t + c_mid * h, z_pred, rhs1, dh, lu, res)
            if self._new_lu is not None:            # stage 1 re-linearised: keep it for stage 2 on
                lu, (J, h_lu) = self._new_lu, self._new_J
                jac_age, fresh, self._new_lu = 0, True, None
            r2 = 0.0
            if ok:
                if sd:
                    f1_id = (z1 - y) / dh
                    rhs2 = y + (1.0 - G_SD) * h * f1_id
                    z2_pred = _truncate(y + (z1 - y) / G_SD, y, self.pos)
                else:
                    rhs2 = C1 * z1 - C2 * y
                    z2_pred = _truncate(y + (z1 - y) / GAMMA, y, self.pos)
                z2, ok, r2 = self._solve_stage(t + h, z2_pred, rhs2, dh, lu, res)
                if self._new_lu is not None:
                    lu, (J, h_lu) = self._new_lu, self._new_J
                    jac_age, fresh, self._new_lu = 0, True, None
            if not ok:
                res.n_newton_fail += 1
                if not fresh:
                    J = None             # refresh Jacobian (+ LU) and retry the same h
                else:
                    # Newton cannot converge at this h even with a fresh Jacobian: the step is
                    # limited by the non-linearity (flow reversals, fronts), not by accuracy.
                    # Remember the limit, so the error controller does not immediately grow h
                    # back into the same failure (it allows growth up to x4 per step).
                    h_cap = min(h_cap, self.cap_after_fail * h)
                    h *= 0.3
                    lu = None
                continue
            # ---- error estimate (filtered by M^-1)
            f1, g1 = m.rhs(t + c_mid * h, z1)
            f2, g2 = m.rhs(t + h, z2)
            if sd:
                y3 = y + h * (SW1 * f_n + SW2 * f1 + SW3 * f2)
            else:
                y3 = y + h * (W1 * f_n + W2 * f1 + W3 * f2)
            err = lu.solve(z2 - y3)
            sc = self.atol + self.rtol * np.maximum(np.abs(y), np.abs(z2))
            self.last_err = err / sc              # (diagnostics) scaled local error of the last attempt
            en = _rms(self.last_err)
            if not np.isfinite(en):
                en = 1e6
            if en > 1.0:
                res.n_reject += 1
                h *= max(0.2, 0.9 * en ** (-1.0 / 3.0))
                continue
            # ---- accept
            if sd:
                acc = acc + h * ((1.0 - G_SD) * g1 + G_SD * g2)
            else:
                acc_g = acc + D * h * (g_n + g1)
                acc = C1 * acc_g - C2 * acc + D * h * g2
            slope_prev = (z2 - y) / h
            t_new = t + h
            if stop_fn is not None and stop_value is not None:
                sv = stop_fn(t_new, z2)
                if sv is not None and (sv - stop_value) * stop_direction >= 0.0:
                    t, y, f_n, g_n = t_new, z2, f2, g2
                    res.stopped_early = True
                    res.stop_value = sv
                    res.n_steps += 1
                    if on_step:
                        on_step(t, y, acc)
                    res.snapshots.append((t, y.copy(), acc.copy()))
                    break
            t, y, f_n, g_n = t_new, z2, f2, g2
            res.n_steps += 1
            jac_age += 1
            fresh = False
            if on_step:
                on_step(t, y, acc)
            if k_out < len(t_out) and t >= t_out[k_out] * (1.0 - 1e-12):
                res.snapshots.append((t, y.copy(), acc.copy()))
                k_out += 1
            gm_step = self._gm[0] - gm0
            if (self.mode == "modified" and max(r1, r2) > 0.5) or gm_step > self.gmres_refresh:
                jac_age = self.jac_max_age          # stale Jacobian / preconditioner -> refresh next step
            fac = 4.0 if en < 1e-8 else min(4.0, max(0.2, 0.9 * en ** (-1.0 / 3.0)))
            if 1.0 <= fac < 1.25:
                fac = 1.0                           # hysteresis: avoid needless refactorisations
            h = max(h * fac, h_desired) if truncated else h * fac
            h_cap *= self.cap_relax                 # the Newton limit is forgotten gradually
            h = min(h, h_cap)
        res.y, res.acc, res.t = y, acc, t
        self.stats["gmres_iters"] = self._gm[0]
        res.nfev = m.nfev
        res.wall_time = time.time() - wall0
        return res
