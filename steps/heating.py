"""Heating / desorption: V1 = CLOSED, V2 = CLOSED, V3 = OPEN_VACUUM (product).

Jacket (and optional internal exchanger) at T_heat; desorbed CO2 + H2O leave through V3 at P_vacuum
and are counted as product.  a.k.a. desorption.py

run_heating(state_in, cfg, grid, props, ...) -> (state_out, StepResult)
Standalone:  python -m steps.heating --state results/<run>/cycle_001/state_after_<previous>.csv
"""
from __future__ import annotations

from ._runner import run_step, step_cli


def _spec(cfg, name=None):
    for s in cfg.cycle:
        if s["type"] == "heating" and (name is None or s["name"] == name):
            return s
    raise KeyError("no 'heating' step in cfg.cycle")


def run_heating(state_in, cfg, grid, props, step=None, **kw):
    """Run the heating step (spec taken from cfg.cycle unless given) and return (state_out, result)."""
    return run_step(state_in, cfg, grid, props, step or _spec(cfg), **kw)


if __name__ == "__main__":
    step_cli("heating")
