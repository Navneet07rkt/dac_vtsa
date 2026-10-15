"""Adsorption: V1 = OPEN (feed blower), V2 = OPEN_ATMOSPHERIC (passive vent), V3 = CLOSED.

Humid air enters at z = 0 with a fixed molar flow; the bottom end is a pressure port at P_atm
(exponential ramp from the bed pressure left by the previous step, so the bed also re-pressurises
through V2 when it starts under vacuum).  Jacket at T_adsorption_jacket.
Optional stop criteria: outlet_ratio (C_out/C_feed) or end_saturation (q_CO2 in the last layer / q*_feed).

run_adsorption(state_in, cfg, grid, props, ...) -> (state_out, StepResult)
Standalone:  python -m steps.adsorption --state results/<run>/cycle_001/state_after_<previous>.csv
"""
from __future__ import annotations

from ._runner import run_step, step_cli


def _spec(cfg, name=None):
    for s in cfg.cycle:
        if s["type"] == "adsorption" and (name is None or s["name"] == name):
            return s
    raise KeyError("no 'adsorption' step in cfg.cycle")


def run_adsorption(state_in, cfg, grid, props, step=None, **kw):
    """Run the adsorption step (spec taken from cfg.cycle unless given) and return (state_out, result)."""
    return run_step(state_in, cfg, grid, props, step or _spec(cfg), **kw)


if __name__ == "__main__":
    step_cli("adsorption")
