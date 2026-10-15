"""Semi-discrete 3D model  dy/dt = f(t, y)  and the implicit-stage residual.

`TVSAModel` assembles, for ONE cycle step (fixed valve states / jacket),
the finite-volume right-hand side from the physics modules:

    thermo (T from e) -> properties -> isotherms + LDF -> Ergun face velocities
    -> species fluxes (mass.py) -> energy fluxes (energy.py) -> wall -> accumulators

State vector (variable-major):
    y = [ c_1..c_ng (ncell each) | q_1..q_na | e | Tw (if wall) ]      "core"
    accumulators (time integrals of boundary/heat flows) are NOT part of the
    Newton system: they do not feed back, so the integrator advances them
    exactly with the same stage formula after each converged stage
    (see core/solver.py).  This keeps the Jacobian sparse.

`stage_residual` implements R(y) = y - gamma h f(t, y) - rhs_known, the
equation solved by Newton in every implicit stage (backward Euler is the
special case gamma = 1, rhs_known = y_old).
"""
from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np

from boundaries.conditions import end_face_fluxes
from physics import mass
from physics.mass import theta_next
from physics.dispersion import Dispersion
from physics.energy import Thermo, Wall, enthalpy_flux
from physics.heat_transfer import HeatTransfer
from physics.isotherms import EquilibriumModel
from physics.kinetics import LDFKinetics
from physics.pressure import Ergun
from numerics.stencils import face_values_z, upwind
from properties import R_GAS
from valves.valve_logic import end_conditions, feed_molar_flux, feed_profile


class TVSAModel:
    def __init__(self, cfg, grid, props, step: Dict, state0, T_jacket_prev: Optional[float] = None):
        self.cfg, self.grid, self.props, self.step = cfg, grid, props, step
        g = grid
        self.species = list(props.species)
        self.ng = len(self.species)
        self.water = bool(cfg.components.water)
        self.adsorbed = ["CO2", "H2O"] if self.water else ["CO2"]
        self.na = len(self.adsorbed)
        self.i_co2 = self.species.index("CO2")
        self.i_h2o = self.species.index("H2O") if self.water else None
        self.i_inert = self.species.index("N2")
        self.ads_to_gas = [self.species.index(s) for s in self.adsorbed]

        # ---- physics modules
        self.thermo = Thermo(cfg, g, props)
        self.eq = EquilibriumModel(cfg, g, props)
        self.kin = LDFKinetics(cfg, g, props)
        self.disp = Dispersion(cfg, g, props)
        self.ht = HeatTransfer(cfg, g, props)
        self.erg = Ergun(cfg, g, props)
        self.wall = Wall(cfg, g, props)
        self.eps_b = cfg.sorbent.bed_voidage
        self.eps_t = cfg.sorbent.total_voidage
        self.rho_b = cfg.derived.bulk_density
        self.cpg = np.asarray(props.cp)
        self.T_ref = cfg.thermo.T_ref
        self.limiter = cfg.solver.limiter
        self.lumped = cfg.transport.heat_transfer_mode == "lumped"

        # ---- boundary conditions of this step
        self.conds = end_conditions(cfg, step, self.species, props)
        self.tau_valve = cfg.operation.valve_time_constant
        P0 = state0.P
        self.P_start = {"top": float(np.mean(P0[:, :, 0])), "bottom": float(np.mean(P0[:, :, -1]))}
        prof = feed_profile(cfg, g)                                   # (Nr, Nth), mean 1
        self.feed_flux = feed_molar_flux(cfg) * prof * g.Az           # mol/s per inlet face
        has_feed = any(c.kind == "feed" for c in self.conds.values())
        self.u_nom = cfg.feed.velocity if has_feed else 0.0           # transverse dispersion scale
        # Nominal molar flux of each (r, theta) column [mol/m2/s], used for the VELOCITY-DEPENDENT
        # parts of axial dispersion and axial thermal dispersion (u = N_nom / C at local
        # conditions).  Evaluating them with the instantaneous, pressure-driven face velocity
        # couples D to grad P, whose response is extremely stiff (du/dP ~ 1/(a dz)): the
        # semi-discrete system then has spurious UNSTABLE eigenvalues (+1e3 1/s were found at a
        # near-stagnant open end) and Newton/step control collapse.  See MODEL.md section 4.
        self.N_nom = (feed_molar_flux(cfg) * prof)[..., None] if has_feed else np.zeros(g.Az.shape + (1,))

        # ---- jacket / wall / internal exchanger
        jk = step["jacket"]
        self.jacket_on = jk["mode"] == "temperature"
        self.T_jac_set = jk.get("T", cfg.ambient.T)
        self.T_jac_prev = self.T_jac_set if T_jacket_prev is None else T_jacket_prev
        self.jac_tau = cfg.jacket.ramp_time
        self.h_jac = cfg.jacket.h_jacket
        self.h_amb = cfg.jacket.h_ambient
        self.T_amb = cfg.ambient.T
        self.P_amb = cfg.ambient.P
        self.cov = g.coverage_fraction(cfg.jacket.theta_coverage_deg, cfg.jacket.z_coverage)
        self.hx_on = bool(cfg.internal_hx.enabled) and self.jacket_on
        self.hA_v = cfg.internal_hx.hA_per_volume
        self.U_lumped = cfg.transport.lumped_U

        # ---- energy accounting
        en = cfg.energy
        self.gam = en.gamma_gas
        self.eta_vac = en.eta_vacuum
        self.vac_model = en.vacuum_efficiency_model
        self.eta_bl = en.eta_blower
        self.pump_water = bool(en.pump_handles_water)

        # ---- layout
        n = g.ncell
        self.n_c = self.ng * n
        self.n_q = self.na * n
        self.off_q = self.n_c
        self.off_e = self.n_c + self.n_q
        self.off_w = self.off_e + n
        self.n_w = g.Nth * g.Nz if self.wall.enabled else 0
        self.n_core = self.off_w + self.n_w
        self.acc_names = self._acc_names()
        self.n_acc = len(self.acc_names)
        self.acc_index = {nm: i for i, nm in enumerate(self.acc_names)}
        self.atol = self._atol()
        self.nfev = 0

    # ------------------------------------------------------------------ layout helpers
    def _acc_names(self) -> List[str]:
        names = []
        for v in ("V1", "V2", "V3"):
            names += [f"n_in_{v}_{s}" for s in self.species]
            names += [f"n_out_{v}_{s}" for s in self.species]
            names += [f"H_in_{v}", f"H_out_{v}", f"W_pump_{v}", f"W_fan_{v}"]
        names += ["Q_jac_pos", "Q_jac_neg", "Q_amb", "Q_wb", "Q_hx_pos", "Q_hx_neg"]
        return names

    def _atol(self) -> np.ndarray:
        s = self.cfg.solver
        n = self.grid.ncell
        tol_c = {"CO2": s.atol_c_co2, "H2O": s.atol_c_h2o, "N2": s.atol_c_inert}
        parts = [np.full(n, tol_c[sp]) for sp in self.species]
        parts += [np.full(n, s.atol_q) for _ in range(self.na)]
        parts.append(np.full(n, self.rho_b * self.cfg.sorbent.cp_solid * s.atol_T))
        if self.n_w:
            parts.append(np.full(self.n_w, s.atol_T))
        return np.concatenate(parts)

    def conservation_matrix(self):
        """Sparse W (ng + 1, n_core): W @ y = inventories [mol of each gas species (gas + adsorbed),
        J of bed + wall energy].  Applied to a stage residual R it gives the stage's global balance
        error (internal fluxes telescope), which the Newton solver drives below tolerance -- see
        core/solver.py, 'conservation polish'."""
        import scipy.sparse as sp
        g = self.grid
        n = g.ncell
        Vc = np.broadcast_to(g.V, g.shape).ravel()
        rows, cols, vals = [], [], []
        for k, s in enumerate(self.species):
            rows.append(np.full(n, k)); cols.append(np.arange(k * n, (k + 1) * n)); vals.append(self.eps_t * Vc)
            if s in self.adsorbed:
                j = self.adsorbed.index(s)
                rows.append(np.full(n, k)); cols.append(self.off_q + np.arange(j * n, (j + 1) * n))
                vals.append(self.rho_b * Vc)
        ke = self.ng
        rows.append(np.full(n, ke)); cols.append(self.off_e + np.arange(n)); vals.append(Vc)
        if self.n_w:
            rows.append(np.full(self.n_w, ke)); cols.append(self.off_w + np.arange(self.n_w))
            vals.append(np.ravel(self.wall.C) * np.ones(self.n_w))
        return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(self.ng + 1, self.n_core))

    def pack(self, st) -> np.ndarray:
        e = self.thermo.energy(st.T, st.c, st.q)
        parts = [st.c.ravel(), st.q.ravel(), e.ravel()]
        if self.n_w:
            parts.append(st.Tw.ravel())
        return np.concatenate(parts)

    def unpack(self, y):
        g = self.grid
        c = y[:self.n_c].reshape((self.ng,) + g.shape)
        q = y[self.off_q:self.off_e].reshape((self.na,) + g.shape)
        e = y[self.off_e:self.off_w].reshape(g.shape)
        Tw = y[self.off_w:self.n_core].reshape(g.wall_shape) if self.n_w else None
        return c, q, e, Tw

    def to_state(self, y, template, t_abs: float, step_name: str):
        c, q, e, Tw = self.unpack(y)
        T = self.thermo.temperature(e, c, q)
        st = template.copy()
        st.c, st.q, st.T = c.copy(), q.copy(), T.copy()
        if Tw is not None:
            st.Tw = Tw.copy()
        else:   # wall not resolved: report the jacket-side boundary temperature of the outer ring
            st.Tw = T[-1].copy()
        st.t = t_abs
        st.step = step_name
        return st

    # ------------------------------------------------------------------ jacket temperature
    def T_jacket(self, t):
        if self.jac_tau <= 0:
            return self.T_jac_set
        return self.T_jac_set + (self.T_jac_prev - self.T_jac_set) * np.exp(-t / self.jac_tau)

    def eta_vacuum(self, P):
        if self.vac_model == "krishnamurthy":
            pb = P / 1e5
            return np.maximum(0.8 * 19.55 * pb / (1.0 + 19.55 * pb), 1e-3)
        return self.eta_vac

    # ------------------------------------------------------------------ RHS
    def rhs(self, t: float, y: np.ndarray, want_acc: bool = True, diag: bool = False):
        """Return (dy_core/dt, dacc/dt) -- or a diagnostics dict if diag=True."""
        self.nfev += 1
        g = self.grid
        c, q, e, Tw = self.unpack(y)
        th = self.thermo
        ng = self.ng
        # ---------- thermodynamic state
        Ctot = c.sum(axis=0)
        T = th.temperature(e, c, q)
        P = Ctot * R_GAS * T
        Tp = np.clip(T, 150.0, 1500.0)      # property/isotherm evaluation guard for wild Newton iterates
        Cpos = np.maximum(Ctot, 1e-9)
        yk = c / Cpos
        ycl = np.clip(yk, 0.0, 1.0)
        mu, kg = self.props.mu_k_mix(Tp, ycl)
        rho = Cpos * self.props.molar_mass_mix(ycl)
        cpm = self.props.cp_mix(ycl)
        Dm = self.props.D_co2_mix(Tp, np.maximum(P, 1.0), ycl)
        lam0 = self.ht.stagnant(kg)
        # ---------- equilibrium + LDF
        pk = c * R_GAS * Tp      # isotherms are smooth/safe for slightly negative iterates
        psat = self.props.p_sat(Tp)
        if self.water:
            qstar = self.eq.equilibrium(pk[self.i_co2], pk[self.i_h2o], Tp, q[1], psat)
        else:
            qstar = self.eq.equilibrium(pk[self.i_co2], None, Tp, None, psat)
        rate = self.kin.rates(qstar, q, Tp)                     # (na, ...)
        # ---------- face velocities (Ergun)
        erg = self.erg
        # z interior
        mu_z = 0.5 * (mu[..., 1:] + mu[..., :-1])
        rho_z = 0.5 * (rho[..., 1:] + rho[..., :-1])
        az, bz = erg.coefficients(mu_z, rho_z)
        uz = erg.velocity(np.diff(P, axis=-1) / g.dz, az, bz)          # (Nr,Nth,Nz-1)
        # r interior
        if g.Nr > 1:
            mu_r = 0.5 * (mu[1:] + mu[:-1])
            rho_r = 0.5 * (rho[1:] + rho[:-1])
            ar, br = erg.coefficients(mu_r, rho_r)
            ur = erg.velocity(np.diff(P, axis=0) / g.drc, ar, br)       # (Nr-1,Nth,Nz)
        else:
            ur = None
        # theta faces (periodic)
        if g.Nth > 1:
            mu_t = 0.5 * (mu + theta_next(mu))
            rho_t = 0.5 * (rho + theta_next(rho))
            at, bt = erg.coefficients(mu_t, rho_t)
            uth = erg.velocity((theta_next(P) - P) / g.arc, at, bt)
        else:
            uth = None
        # ---------- boundary (valve) fluxes
        bnd = {}
        for end, k in (("top", 0), ("bottom", -1)):
            layer = {"c": c[..., k], "T": T[..., k], "P": P[..., k], "mu": mu[..., k], "rho": rho[..., k]}
            bnd[end] = end_face_fluxes(self.conds[end], t, self.P_start[end], layer, g.Az, 0.5 * g.dz,
                                       erg, self.cpg, self.T_ref, self.feed_flux, self.tau_valve)
            if end == "bottom":   # outward normal is +z there: flux in +z = -F_in
                pass
        # ---------- species fluxes
        Cz = 0.5 * (Cpos[..., 1:] + Cpos[..., :-1])
        u_disp = self.N_nom / Cz                                       # nominal local velocity (see __init__)
        Dz = self.disp.axial(0.5 * (Dm[..., 1:] + Dm[..., :-1]), u_disp, rho_z, mu_z)
        Fz_conv = mass.convective_z(c, uz, g.Az, self.limiter)
        Jz = mass.dispersive_z(yk, Cz, Dz, self.eps_b, g.Az, g.dz)
        Fz_full = np.concatenate([bnd["top"]["F_in"][..., None], Fz_conv + Jz,
                                  -bnd["bottom"]["F_in"][..., None]], axis=-1)
        if ur is not None:
            Cr = 0.5 * (Cpos[1:] + Cpos[:-1])
            Dr = self.disp.transverse(0.5 * (Dm[1:] + Dm[:-1]), self.u_nom)
            Fr_conv = mass.convective_r(c, ur, g.Ar)
            Jr = mass.dispersive_r(yk, Cr, Dr, self.eps_b, g.Ar, g.drc)
            Fr = Fr_conv + Jr
        else:
            Fr_conv = Jr = Fr = None
        if uth is not None:
            Ct = 0.5 * (Cpos + theta_next(Cpos))
            Dt = self.disp.transverse(0.5 * (Dm + theta_next(Dm)), self.u_nom)
            Ft_conv = mass.convective_theta(c, uth, g.Ath)
            Jt = mass.dispersive_theta(yk, Ct, Dt, self.eps_b, g.Ath, g.arc)
            Ft = Ft_conv + Jt
        else:
            Ft_conv = Jt = Ft = None
        div_c = mass.divergence(Fz_full, Fr, Ft, g.Nth)
        sink = np.zeros_like(c)
        for j, kgas in enumerate(self.ads_to_gas):
            sink[kgas] = rate[j]
        dcdt = -div_c / (self.eps_t * g.V) - (self.rho_b / self.eps_t) * sink
        # ---------- energy fluxes
        cpg, Tr = self.cpg, self.T_ref
        Tz_f = face_values_z(T, uz, self.limiter)
        Tz_avg = 0.5 * (T[..., 1:] + T[..., :-1])
        cpz = 0.5 * (cpm[..., 1:] + cpm[..., :-1])
        lam_z = self.ht.axial(0.5 * (lam0[..., 1:] + lam0[..., :-1]), Cz, cpz, u_disp)
        Hz_int = (enthalpy_flux(Fz_conv, Tz_f, cpg, Tr) + enthalpy_flux(Jz, Tz_avg, cpg, Tr)
                  - lam_z * np.diff(T, axis=-1) / g.dz * g.Az[:, :, None])
        Hz_full = np.concatenate([bnd["top"]["H_in"][..., None], Hz_int,
                                  -bnd["bottom"]["H_in"][..., None]], axis=-1)
        if ur is not None:
            Tr_up = upwind(T[:-1], T[1:], ur)
            Tr_avg = 0.5 * (T[1:] + T[:-1])
            lam_r = self.ht.radial(0.5 * (lam0[1:] + lam0[:-1]), Cr, 0.5 * (cpm[1:] + cpm[:-1]), self.u_nom)
            Hr = (enthalpy_flux(Fr_conv, Tr_up, cpg, Tr) + enthalpy_flux(Jr, Tr_avg, cpg, Tr)
                  - lam_r * np.diff(T, axis=0) / g.drc * g.Ar)
        else:
            Hr = None
        if uth is not None:
            Tn = theta_next(T)
            Tt_up = upwind(T, Tn, uth)
            lam_t = self.ht.radial(0.5 * (lam0 + theta_next(lam0)), Ct,
                                   0.5 * (cpm + theta_next(cpm)), self.u_nom)
            Ht = (enthalpy_flux(Ft_conv, Tt_up, cpg, Tr) + enthalpy_flux(Jt, 0.5 * (T + Tn), cpg, Tr)
                  - lam_t * (Tn - T) / g.arc * g.Ath)
        else:
            Ht = None
        div_e = mass.divergence(Hz_full, Hr, Ht, g.Nth)
        # ---------- wall / jacket / exchanger
        To = T[-1]                                                         # outer ring (Nth, Nz)
        lam_out = self.ht.radial(lam0[-1], Cpos[-1], cpm[-1], self.u_nom)
        if self.lumped:
            R_bed = 1.0 / self.U_lumped
        else:
            alpha_w = self.ht.wall_coefficient(lam0[-1], kg[-1], rho[-1], mu[-1],
                                               cpm[-1] / np.maximum(rho[-1] / Cpos[-1], 1e-6), self.u_nom)
            R_bed = 1.0 / alpha_w + g.half_dr_outer / lam_out
        Tj = self.T_jacket(t)
        if self.wall.enabled:
            Q_jac = self.cov * self.h_jac * g.A_wall_out * (Tj - Tw) if self.jacket_on else np.zeros_like(Tw)
            Q_amb = (1.0 - self.cov) * self.h_amb * g.A_wall_out * (self.T_amb - Tw)
            Q_wb = g.A_wall_in * (Tw - To) / R_bed
            dTw = (Q_jac + Q_amb - Q_wb + self.wall.conduction(Tw)) / self.wall.C
        else:
            if self.jacket_on and self.h_jac > 0:
                R_j = R_bed if self.lumped else R_bed + 1.0 / self.h_jac
                Q_jac = self.cov * g.A_wall_in * (Tj - To) / R_j
            else:
                Q_jac = np.zeros_like(To)
            Q_amb = ((1.0 - self.cov) * g.A_wall_in * (self.T_amb - To) / (R_bed + 1.0 / self.h_amb)
                     if self.h_amb > 0 else np.zeros_like(To))
            Q_wb = Q_jac + Q_amb
            dTw = None
        src = np.zeros_like(T)
        src[-1] += Q_wb
        if self.hx_on:
            Q_hx = self.hA_v * g.V * (Tj - T)
            src += Q_hx
        else:
            Q_hx = None
        dedt = (-div_e + src) / g.V
        # ---------- outputs
        dq = rate
        parts = [dcdt.ravel(), dq.ravel(), dedt.ravel()]
        if dTw is not None:
            parts.append(dTw.ravel())
        f = np.concatenate(parts)
        if diag:
            return {"T": T, "P": P, "y": yk, "qstar": qstar, "uz": uz, "ur": ur, "uth": uth, "bnd": bnd,
                    "Q_jac": Q_jac, "Q_wb": Q_wb, "Q_amb": Q_amb, "Q_hx": Q_hx, "R_bed": R_bed}
        if not want_acc:
            return f, None
        # ---------- accumulator rates
        acc = np.zeros(self.n_acc)
        ai = self.acc_index
        for end in ("top", "bottom"):
            cond = self.conds[end]
            if cond.valve is None:
                continue
            v = cond.valve
            b = bnd[end]
            inflow = b["u_out"] < 0.0
            Fin = b["F_in"]
            if "F_inpart" in b:          # flux-split valve: in- and out-going streams counted separately
                Fin_p, Fout_p = b["F_inpart"], b["F_outpart"]
                for k, s in enumerate(self.species):
                    acc[ai[f"n_in_{v}_{s}"]] = np.sum(Fin_p[k])
                    acc[ai[f"n_out_{v}_{s}"]] = np.sum(Fout_p[k])
                acc[ai[f"H_in_{v}"]] = np.sum(b["H_inpart"])
                acc[ai[f"H_out_{v}"]] = np.sum(b["H_outpart"])
            else:
                for k, s in enumerate(self.species):
                    acc[ai[f"n_in_{v}_{s}"]] = np.sum(np.where(inflow, Fin[k], 0.0))
                    acc[ai[f"n_out_{v}_{s}"]] = np.sum(np.where(inflow, 0.0, -Fin[k]))
                acc[ai[f"H_in_{v}"]] = np.sum(np.where(inflow, b["H_in"], 0.0))
                acc[ai[f"H_out_{v}"]] = np.sum(np.where(inflow, 0.0, -b["H_in"]))
            ex = (self.gam - 1.0) / self.gam
            pre = self.gam / (self.gam - 1.0) * R_GAS
            if cond.is_vacuum:
                Fout = np.where(inflow, 0.0, -Fin)
                if not self.pump_water and self.i_h2o is not None:
                    Fout = np.delete(Fout, self.i_h2o, axis=0)
                Ftot = Fout.sum(axis=0)
                k_lay = 0 if end == "top" else -1
                ratio = np.maximum(self.P_amb / np.maximum(b["P_face"], 1.0), 1.0)
                acc[ai[f"W_pump_{v}"]] = np.sum(pre * T[..., k_lay] * Ftot * (ratio ** ex - 1.0)
                                                / self.eta_vacuum(b["P_face"]))
            if cond.is_feed:
                Ftot = b["F_inpart"].sum(axis=0) if "F_inpart" in b else np.where(inflow, Fin, 0.0).sum(axis=0)
                ratio = np.maximum(b["P_face"] / self.P_amb, 1.0)
                acc[ai[f"W_fan_{v}"]] = np.sum(pre * cond.T_ext * Ftot * (ratio ** ex - 1.0) / self.eta_bl)
        acc[ai["Q_jac_pos"]] = np.sum(np.maximum(Q_jac, 0.0))
        acc[ai["Q_jac_neg"]] = np.sum(np.minimum(Q_jac, 0.0))
        acc[ai["Q_amb"]] = np.sum(Q_amb)
        acc[ai["Q_wb"]] = np.sum(Q_wb)
        if Q_hx is not None:
            acc[ai["Q_hx_pos"]] = np.sum(np.maximum(Q_hx, 0.0))
            acc[ai["Q_hx_neg"]] = np.sum(np.minimum(Q_hx, 0.0))
        return f, acc

    # ------------------------------------------------------------------ residual for implicit stages
    def stage_residual(self, t, y, gamma_h, rhs_known):
        """R(y) = y - gamma_h * f(t, y) - rhs_known   (core variables only)."""
        f, _ = self.rhs(t, y, want_acc=False)
        return y - gamma_h * f - rhs_known

    # ------------------------------------------------------------------ inventories / diagnostics
    def inventories(self, y) -> Dict[str, float]:
        """Moles of each species (gas + adsorbed) and energy of bed + wall."""
        g = self.grid
        c, q, e, Tw = self.unpack(y)
        out = {}
        V = g.V
        for k, s in enumerate(self.species):
            gas = float(np.sum(self.eps_t * c[k] * V))
            ads = 0.0
            if s in self.adsorbed:
                ads = float(np.sum(self.rho_b * q[self.adsorbed.index(s)] * V))
            out[f"gas_{s}"] = gas
            out[f"ads_{s}"] = ads
            out[f"tot_{s}"] = gas + ads
        out["E_bed"] = float(np.sum(e * V))
        out["E_wall"] = self.wall.energy(Tw) if Tw is not None else 0.0
        out["E_total"] = out["E_bed"] + out["E_wall"]
        return out

    def stop_metric(self, t, y) -> Optional[float]:
        stop = self.step.get("stop_on")
        if not stop:
            return None
        crit = stop["criterion"]
        g = self.grid
        c, q, e, Tw = self.unpack(y)
        if crit in ("min_temperature", "max_temperature"):
            T = self.thermo.temperature(e, c, q)
            return float(np.sum(T * g.V) / np.sum(g.V * np.ones_like(T)))
        if crit == "end_saturation":
            # mean CO2 loading in the last axial layer relative to equilibrium with the feed
            return float(np.mean(q[0][..., -1]) / self.q_feed_eq())
        if crit == "outlet_ratio":
            d = self.rhs(t, y, diag=True)
            b = d["bnd"]["bottom"] if self.conds["top"].kind == "feed" else d["bnd"]["top"]
            Fout = -b["F_in"]
            tot = Fout.sum(axis=0)
            m = tot > 0
            if not np.any(m):
                return 0.0
            y_out = Fout[self.i_co2][m].sum() / tot[m].sum()
            yf = self.conds["top"].y_ext[self.i_co2] if self.conds["top"].y_ext is not None else 4e-4
            return float(y_out / yf)
        return None

    def q_feed_eq(self) -> float:
        cfg = self.cfg
        from valves.valve_logic import humid_air_composition
        T = cfg.feed.T
        psat = float(self.props.p_sat(T))
        yv = humid_air_composition(self.species, cfg.feed.P, T, cfg.feed.RH if self.water else 0.0,
                                   cfg.feed.y_co2_dry, psat)
        p = yv * cfg.feed.P
        if self.water:
            qw = self.eq.q_h2o(p[self.i_h2o], T, psat)
            return float(self.eq.q_co2(p[self.i_co2], T, qw))
        return float(self.eq.q_co2_dry(p[self.i_co2], T))
