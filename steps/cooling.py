"""Cooling: V1 = CLOSED, V2 = CLOSED, V3 = CLOSED (default) or OPEN_VACUUM (design choice).

Jacket at T_cool.  With V3 closed the bed pressure falls as the gas cools and CO2/H2O re-adsorb;
the next adsorption step re-pressurises it.

run_cooling(state_in, cfg, grid, props, ...) -> (state_out, StepResult)
Standalone:  python -m steps.cooling --state results/<run>/cycle_001/state_after_<previous>.csv
"""
from __future__ import annotations

from ._runner import run_step, step_cli


def _spec(cfg, name=None):
    for s in cfg.cycle:
        if s["type"] == "cooling" and (name is None or s["name"] == name):
            return s
    raise KeyError("no 'cooling' step in cfg.cycle")


def run_cooling(state_in, cfg, grid, props, step=None, **kw):
    """Run the cooling step (spec taken from cfg.cycle unless given) and return (state_out, result)."""
    return run_step(state_in, cfg, grid, props, step or _spec(cfg), **kw)


if __name__ == "__main__":
    step_cli("cooling")
