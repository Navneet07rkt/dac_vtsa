"""Blowdown: V1 = CLOSED, V2 = OPEN_VACUUM (ramps P_atm -> P_vacuum), V3 = CLOSED.

Removes the air (N2/O2) held in the voids before heating; everything leaving through V2 is vented
(counted as a CO2 loss).  Jacket normally still at the adsorption temperature.

run_blowdown(state_in, cfg, grid, props, ...) -> (state_out, StepResult)
Standalone:  python -m steps.blowdown --state results/<run>/cycle_001/state_after_<previous>.csv
"""
from __future__ import annotations

from ._runner import run_step, step_cli


def _spec(cfg, name=None):
    for s in cfg.cycle:
        if s["type"] == "blowdown" and (name is None or s["name"] == name):
            return s
    raise KeyError("no 'blowdown' step in cfg.cycle")


def run_blowdown(state_in, cfg, grid, props, step=None, **kw):
    """Run the blowdown step (spec taken from cfg.cycle unless given) and return (state_out, result)."""
    return run_step(state_in, cfg, grid, props, step or _spec(cfg), **kw)


if __name__ == "__main__":
    step_cli("blowdown")
