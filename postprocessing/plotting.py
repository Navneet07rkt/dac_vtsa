"""All plotting helpers (matplotlib, PNG).  Pure presentation -- reads results only.

* outlet / breakthrough curves                (plot_outlet)
* z-r slice maps of q, T, P at chosen times     (plot_zr_slice)
* polar r-theta maps at a chosen height z       (plot_polar)
* time histories at probe points over a cycle   (plot_probe_histories)
* CSS convergence and energy breakdown          (plot_css, plot_energy_breakdown)
* tornado chart of sensitivities                (plot_tornado)

Colours: categorical slots in fixed order, single-hue sequential ramp for
fields (magnitude), one y-axis per panel.
"""
from __future__ import annotations

import os
from typing import Dict, List, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
SEQ = LinearSegmentedColormap.from_list(
    "seq_blue", ["#f4f8fd", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])
SEQ_WARM = LinearSegmentedColormap.from_list(
    "seq_orange", ["#fdf6f1", "#fbd9c6", "#f6b08e", "#ef8a5c", "#eb6834", "#c4501f", "#933a14", "#5e240c"])

plt.rcParams.update({
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "axes.edgecolor": INK2,
    "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 9, "axes.titlesize": 10, "legend.frameon": False, "lines.linewidth": 1.6,
    "axes.prop_cycle": matplotlib.cycler(color=SERIES),
})


def _save(fig, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    try:
        fig.tight_layout()
    except Exception:            # colorbars spanning several axes can defeat tight_layout
        pass
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_outlet(rows: List[Dict], path: str, title: str = ""):
    """Four stacked panels over the cycle: C_out/C_feed, flows, pressure, bed T.

    C_out/C_feed is drawn only while the valve at that end is open (a closed end has no
    outlet stream); during back-flow the end shows the gas sitting in the end cells."""
    if not rows:
        return None
    t = np.array([r["t_cycle_s"] for r in rows]) / 3600.0
    steps = [r["step"] for r in rows]
    fig, ax = plt.subplots(4, 1, figsize=(8, 9), sharex=True)
    for end, col, lab in (("bottom", SERIES[0], "bottom end (V2/V3)"), ("top", SERIES[1], "top end (V1)")):
        # only ends that are open and not a feed inlet carry an outlet stream
        ratio = np.array([r[f"{end}_C_over_Cfeed"]
                          if r[f"{end}_valve"] != "closed" and r.get(f"{end}_valve_state", "") != "OPEN"
                          else np.nan for r in rows], dtype=float)
        ax[0].plot(t, np.clip(ratio, 1e-5, None), color=col, label=lab)
    ax[0].set_yscale("log")
    ax[0].set_ylabel("C_CO2,out / C_CO2,feed")
    ax[0].legend(loc="lower right", frameon=False)
    for end, col in (("bottom", SERIES[0]), ("top", SERIES[1])):
        ax[1].plot(t, [r[f"{end}_flow_out_mol_s"] * 1e3 for r in rows], color=col, label=f"{end} (+ = leaving)")
    ax[1].set_ylabel("molar flow [mmol/s]")
    ax[1].legend(loc="lower right", frameon=False)
    ax[2].plot(t, [r["bed_P_mean_Pa"] / 1e3 for r in rows], color=SERIES[0], label="bed mean")
    ax[2].set_ylabel("bed mean P [kPa]")
    ax[3].plot(t, [r["bed_T_mean_K"] - 273.15 for r in rows], color=SERIES[1], label="bed mean T")
    ax[3].set_ylabel("bed mean T [°C]")
    ax[3].set_xlabel("time in cycle [h]")
    # step boundaries; labels staggered so that short steps (60 s blowdown) do not overprint
    last = None
    k = 0
    for ti, s_ in zip(t, steps):
        if s_ != last:
            for a in ax:
                a.axvline(ti, color=GRID, lw=1.0, zorder=0)
            ax[3].text(ti, 0.98 - 0.12 * (k % 3), f" {s_}", va="top", ha="left", color=INK2, fontsize=8,
                       transform=ax[3].get_xaxis_transform())
            last = s_
            k += 1
    if title:
        ax[0].set_title(title, loc="left", color=INK)
    return _save(fig, path)


def _field(state, name):
    if name == "T":
        return state.T - 273.15, "T [°C]"
    if name == "P":
        return state.P / 1e3, "P [kPa]"
    if name.startswith("q_"):
        return state.q[state.adsorbed.index(name[2:])], f"{name} [mol/kg]"
    if name.startswith("y_"):
        return state.y[state.species.index(name[2:])], name
    raise ValueError(name)


def plot_zr_slice(grid, snaps: Sequence[Dict], field: str, path: str, theta_index: int = 0):
    """Contour maps in the (z, r) half-plane at theta_j for several snapshots (+ the
    opposite half-plane theta_j + 180 deg mirrored, to show a full diameter)."""
    n = len(snaps)
    if n == 0:
        return None
    fig, axes = plt.subplots(1, n, figsize=(2.6 * n + 1, 3.2), squeeze=False)
    vals = [_field(s["state"], field)[0] for s in snaps]
    vmin = min(float(v.min()) for v in vals)
    vmax = max(float(v.max()) for v in vals)
    if vmax - vmin < 1e-12:
        vmax = vmin + 1e-12
    j_opp = (theta_index + grid.Nth // 2) % grid.Nth
    cmap = SEQ_WARM if field == "T" else SEQ
    for a, s, v in zip(axes[0], snaps, vals):
        half = v[:, theta_index, :]
        other = v[:, j_opp, :] if grid.Nth > 1 else half
        full = np.concatenate([other[::-1], half], axis=0)            # r from -R to R
        rr = np.concatenate([-grid.r_c[::-1], grid.r_c]) * 100
        zz = grid.z_c * 100
        m = a.pcolormesh(zz, rr, full, cmap=cmap, vmin=vmin, vmax=vmax, shading="nearest")
        a.set_title(f"{s['step']} t={s['t_step']:.0f}s", fontsize=8)
        a.set_xlabel("z [cm]")
        a.grid(False)
    axes[0][0].set_ylabel("r [cm]  (theta and theta+180)")
    fig.colorbar(m, ax=axes[0].tolist(), label=_field(snaps[0]["state"], field)[1], shrink=0.9)
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def plot_polar(grid, snap: Dict, field: str, path: str, z_index: int = None):
    k = grid.Nz // 2 if z_index is None else z_index
    v, label = _field(snap["state"], field)
    sl = v[:, :, k]
    th_edges = grid.th_f
    r_edges = grid.r_f * 100
    fig = plt.figure(figsize=(4.4, 4))
    a = fig.add_subplot(111, projection="polar")
    TH, RR = np.meshgrid(th_edges, r_edges)
    m = a.pcolormesh(TH, RR, sl, cmap=SEQ_WARM if field == "T" else SEQ, shading="flat")
    a.set_title(f"{field} at z = {grid.z_c[k] * 100:.2f} cm ({snap['step']}, t={snap['t_step']:.0f}s)", fontsize=9)
    fig.colorbar(m, ax=a, label=label, shrink=0.8)
    return _save(fig, path)


def plot_probe_histories(grid, snaps: Sequence[Dict], path: str, field: str = "T"):
    """Field at axis / mid-radius / wall-side cells (mid-height, theta=0) across the cycle."""
    if not snaps:
        return None
    k = grid.Nz // 2
    picks = [(0, "axis"), (grid.Nr // 2, "mid-radius"), (grid.Nr - 1, "wall side")]
    t = np.array([s["t_cycle"] for s in snaps]) / 3600.0
    fig, ax = plt.subplots(figsize=(7.5, 3.4))
    for (i, lab), col in zip(picks, SERIES):
        ax.plot(t, [_field(s["state"], field)[0][i, 0, k] for s in snaps], color=col,
                label=f"{lab} (r={grid.r_c[i] * 100:.2f} cm)", marker="o", ms=3)
    ax.set_xlabel("time in cycle [h]")
    ax.set_ylabel(_field(snaps[0]["state"], field)[1])
    ax.legend()
    return _save(fig, path)


def plot_css(cycles, path: str):
    n = [c.cycle for c in cycles]
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.2))
    ax[0].semilogy(n, [max(c.change["max"], 1e-16) for c in cycles], marker="o", color=SERIES[0], label="state change")
    ax[0].semilogy(n, [max(c.co2_balance, 1e-16) for c in cycles], marker="s", color=SERIES[1], label="CO2 cycle balance")
    ax[0].set_xlabel("cycle")
    ax[0].set_ylabel("relative measure")
    ax[0].legend()
    ax[1].plot(n, [c.kpis["recovery_pct"] for c in cycles], marker="o", color=SERIES[0], label="recovery [%]")
    ax[1].plot(n, [c.kpis["purity_dry_pct"] for c in cycles], marker="s", color=SERIES[1], label="purity [%]")
    ax[1].set_xlabel("cycle")
    ax[1].legend()
    return _save(fig, path)


def plot_energy_breakdown(row: Dict, path: str, kg_co2: float):
    keys = ["desorption_CO2_J", "desorption_H2O_J", "bed_sensible_J", "wall_sensible_J",
            "gas_enthalpy_out_J", "ambient_loss_J"]
    labels = ["CO2 desorption", "H2O desorption", "bed sensible", "wall sensible", "gas enthalpy out", "ambient loss"]
    vals = [row[k] / max(kg_co2, 1e-30) / 1e3 for k in keys]
    fig, ax = plt.subplots(figsize=(7, 3.2))
    ax.barh(labels[::-1], vals[::-1], color=SERIES[0], height=0.6)
    for y, v in enumerate(vals[::-1]):
        ax.text(v, y, f" {v:.0f}", va="center", color=INK2, fontsize=8)
    ax.set_xlabel(f"kJ per kg CO2 product  (step '{row['step']}', total {row['heat_input_J'] / max(kg_co2, 1e-30) / 1e3:.0f})")
    ax.grid(axis="y", visible=False)
    return _save(fig, path)


def plot_tornado(rows: List[Dict], kpi: str, path: str):
    rows = [r for r in rows if r["kpi"] == kpi]
    rows.sort(key=lambda r: abs(r["normalized_sensitivity"]))
    if not rows:
        return None
    fig, ax = plt.subplots(figsize=(7, 0.35 * len(rows) + 1.2))
    vals = [r["normalized_sensitivity"] for r in rows]
    ax.barh([r["parameter"] for r in rows], vals, color=[SERIES[0] if v >= 0 else SERIES[1] for v in vals],
            height=0.6)
    ax.axvline(0, color=INK2, lw=0.8)
    ax.set_xlabel(f"normalised sensitivity  d ln({kpi}) / d ln(p)")
    ax.grid(axis="y", visible=False)
    return _save(fig, path)
