"""Unit tests: configuration, isotherms (vs hand calculations from Young et al. 2021),
conservation of the semi-discrete operator, Jacobian colouring, CSV round trip,
axisymmetry, explicit vs implicit agreement, Fornberg stencils.

Run:  python -m pytest -q tests/
"""
import math
import os
import tempfile

import numpy as np
import pytest

from config.parameters import load_config, with_overrides
from config.schema import ConfigError
from core.grid import Grid
from core.residual import TVSAModel
from core.state import initial_state
from numerics.sparsity import build_pattern, color_columns, verify_coloring
from numerics.stencils import fd_weights
from physics.isotherms import EquilibriumModel
from properties import arden_buck, get_backend
from simio.checkpoint import load_state_csv, save_state_csv


def small_cfg(**over):
    o = {"grid.Nr": 3, "grid.Ntheta": 4, "grid.Nz": 6}
    o.update(over)
    return load_config(override_dict=o)


# ---------------------------------------------------------------- configuration
def test_defaults_load_and_derived():
    cfg = load_config()
    assert math.isclose(cfg.derived.bulk_density, 880 * 0.6)
    assert math.isclose(cfg.derived.particle_porosity, (0.54 - 0.40) / 0.60)
    assert cfg.cycle[2]["jacket"]["T"] == cfg.operation.T_heat


@pytest.mark.parametrize("key,val", [("geometry.bed_radius", -1.0), ("sorbent.bed_voidage", 1.2),
                                     ("grid.Nr", 0), ("isotherm.coadsorption.model", "magic")])
def test_invalid_values_rejected(key, val):
    with pytest.raises(ConfigError):
        load_config(override_dict={key: val})


def test_unknown_key_rejected():
    with pytest.raises(ConfigError):
        load_config(override_dict={"sorbent.bed_vodage": 0.4})


# ---------------------------------------------------------------- isotherms (Young 2021 parameters)
def test_isotherm_values_match_hand_calculation():
    cfg = load_config()
    eq = EquilibriumModel(cfg)
    T = 288.15
    ps = float(arden_buck(T))
    pw = 0.55 * ps
    pc = 4e-4 * (101325 - pw)
    qw = float(eq.q_h2o(pw, T, ps))
    assert abs(qw - 4.091) < 5e-3                                   # GAB, 55 % RH, 15 C
    assert abs(float(eq.q_co2_dry(pc, T)) - 1.0226) < 2e-3          # dry Toth
    assert abs(float(eq.q_co2(pc, T, qw)) - 1.511) < 3e-3           # WADST
    mech = EquilibriumModel(with_overrides(cfg, {"isotherm.coadsorption.model": "mechanistic"}))
    assert abs(float(mech.q_co2(pc, T, qw)) - 1.462) < 3e-3


def test_isotherm_safe_for_negative_iterates():
    eq = EquilibriumModel(load_config())
    q = eq.q_co2(np.array([-1.0, 0.0, 1e-9]), 300.0, np.array([0.0, 1.0, 2.0]))
    assert np.all(np.isfinite(q)) and q[0] == 0.0 and q[1] == 0.0


def test_saturation_pressure():
    assert abs(float(arden_buck(373.15)) - 101325) / 101325 < 3e-3


# ---------------------------------------------------------------- semi-discrete conservation
@pytest.mark.parametrize("step_index", [0, 1, 2, 3])
def test_rhs_conserves_species_and_energy(step_index):
    cfg = small_cfg(**{"jacket.theta_coverage_deg": [0.0, 180.0], "feed.profile": "cosine_theta",
                       "feed.profile_amplitude": 0.3})
    g = Grid.from_config(cfg)
    props = get_backend(cfg)
    st = initial_state(cfg, g, props)
    st.T = st.T + 5.0 * np.random.default_rng(1).random(st.T.shape)    # make it non-trivial
    m = TVSAModel(cfg, g, props, cfg.cycle[step_index], st)
    y = m.pack(st)
    f, acc = m.rhs(1.0, y)
    ai = m.acc_index
    dc = f[:m.n_c].reshape((m.ng,) + g.shape)
    dq = f[m.off_q:m.off_e].reshape((m.na,) + g.shape)
    for k, s in enumerate(m.species):
        dinv = np.sum(m.eps_t * dc[k] * g.V)
        if s in m.adsorbed:
            dinv += np.sum(m.rho_b * dq[m.adsorbed.index(s)] * g.V)
        net = sum(acc[ai[f"n_in_{v}_{s}"]] - acc[ai[f"n_out_{v}_{s}"]] for v in ("V1", "V2", "V3"))
        assert abs(dinv - net) <= 1e-12 * max(1.0, abs(net)) + 1e-18
    de = f[m.off_e:m.off_w].reshape(g.shape)
    dTw = f[m.off_w:m.n_core].reshape(g.wall_shape)
    dE = np.sum(de * g.V) + np.sum(m.wall.C * dTw)
    netH = (sum(acc[ai[f"H_in_{v}"]] - acc[ai[f"H_out_{v}"]] for v in ("V1", "V2", "V3"))
            + acc[ai["Q_jac_pos"]] + acc[ai["Q_jac_neg"]] + acc[ai["Q_amb"]] + acc[ai["Q_hx_pos"]] + acc[ai["Q_hx_neg"]])
    assert abs(dE - netH) <= 1e-9 * max(1.0, abs(netH))


def test_temperature_energy_roundtrip():
    cfg = small_cfg()
    g = Grid.from_config(cfg)
    props = get_backend(cfg)
    st = initial_state(cfg, g, props)
    m = TVSAModel(cfg, g, props, cfg.cycle[0], st)
    c, q, e, _ = m.unpack(m.pack(st))
    assert np.max(np.abs(m.thermo.temperature(e, c, q) - st.T)) < 1e-9


# ---------------------------------------------------------------- Jacobian structure
def test_coloring_is_valid_and_jacobian_matches_dense_fd():
    cfg = small_cfg()
    g = Grid.from_config(cfg)
    props = get_backend(cfg)
    P = build_pattern(g.Nr, g.Nth, g.Nz, 3, 2, True)
    col = color_columns(P)
    assert verify_coloring(P, col)
    from core.jacobian import FDJacobian
    st = initial_state(cfg, g, props)
    st.T = st.T + np.random.default_rng(0).random(st.T.shape)
    m = TVSAModel(cfg, g, props, cfg.cycle[0], st)
    y = m.pack(st)
    J = FDJacobian(m)(5.0, y).toarray()
    f0, _ = m.rhs(5.0, y, want_acc=False)
    for j in np.random.default_rng(2).choice(len(y), 12, replace=False):
        h = 1e-7 * max(abs(y[j]), 1e-3)
        yp = y.copy()
        yp[j] += h
        col_fd = (m.rhs(5.0, yp, want_acc=False)[0] - f0) / h
        scale = max(np.abs(col_fd).max(), 1e-12)
        assert np.abs(col_fd - J[:, j]).max() / scale < 1e-3


# ---------------------------------------------------------------- I/O
def test_state_csv_roundtrip_is_exact():
    cfg = small_cfg()
    g = Grid.from_config(cfg)
    props = get_backend(cfg)
    st = initial_state(cfg, g, props)
    st.T = st.T + np.random.default_rng(3).random(st.T.shape)
    st.q = st.q + np.random.default_rng(4).random(st.q.shape)
    with tempfile.TemporaryDirectory() as d:
        p = save_state_csv(st, g, os.path.join(d, "s.csv"))
        st2 = load_state_csv(p, g, props.species, st.adsorbed)
    assert np.array_equal(st.T, st2.T) and np.array_equal(st.c, st2.c) and np.array_equal(st.q, st2.q)
    assert np.array_equal(st.Tw, st2.Tw)


# ---------------------------------------------------------------- stencils
def test_fornberg_weights():
    w = fd_weights(0.0, np.array([-1.0, 0.0, 1.0]), 2)
    assert np.allclose(w, [1.0, -2.0, 1.0])
    w = fd_weights(0.0, np.array([-0.5, 0.0, 1.5]), 1)
    x = np.array([-0.5, 0.0, 1.5])
    assert abs(np.dot(w, x ** 2)) < 1e-12 and abs(np.dot(w, x) - 1.0) < 1e-12


# ---------------------------------------------------------------- symmetry and integrator cross-check
def test_axisymmetric_state_stays_axisymmetric():
    """Uniform feed + full jacket: a theta-uniform state must have a theta-uniform time derivative."""
    cfg = small_cfg()
    g = Grid.from_config(cfg)
    props = get_backend(cfg)
    st = initial_state(cfg, g, props)
    rng = np.random.default_rng(5)
    st.T = st.T + 5.0 * rng.random((g.Nr, 1, g.Nz))          # varies in r and z only
    m = TVSAModel(cfg, g, props, cfg.cycle[0], st)
    f, _ = m.rhs(10.0, m.pack(st), want_acc=False)
    c, q, e, Tw = m.unpack(f)
    for arr in (c, q, e[None]):
        spread = np.abs(arr - arr.mean(axis=-2, keepdims=True)).max()
        assert spread <= 1e-9 * max(np.abs(arr).max(), 1e-30)


def test_implicit_solver_matches_rk4_on_short_transient():
    """SDIRK2/Newton-Krylov vs classical RK4 (tiny grid, 2 ms of blowdown, explicit-stable dt)."""
    from core.ode_solver import integrate_rk4
    from core.solver import TRBDF2
    from numerics.cfl import explicit_dt
    cfg = load_config(override_dict={"grid.Nr": 2, "grid.Ntheta": 1, "grid.Nz": 4,
                                     "solver.rtol": 1e-6, "solver.newton_kappa": 1e-3})
    g = Grid.from_config(cfg)
    props = get_backend(cfg)
    st = initial_state(cfg, g, props)
    m = TVSAModel(cfg, g, props, cfg.cycle[1], st)          # blowdown: pump on, fast pressure transient
    y0 = m.pack(st)
    t_end = 2e-3
    dt = explicit_dt(m, y0, safety=0.25)
    y_rk = integrate_rk4(m, y0, t_end, dt)[0]
    res = TRBDF2(m, cfg, h_max=t_end, h_init=1e-6).integrate(y0, t_end, np.array([t_end]))
    sc = m.atol + 1e-4 * np.abs(y_rk)
    err = np.sqrt(np.mean(((res.y - y_rk) / sc) ** 2))
    assert err < 1.0, err
