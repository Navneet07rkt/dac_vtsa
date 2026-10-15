"""Pressurisation (optional, not in the default 4-step cycle): V1 = OPEN_ATMOSPHERIC with feed gas,
V2 = CLOSED, V3 = CLOSED.

Used to reproduce Young et al. (2021)'s 6-step cycle.

run_pressurization(state_in, cfg, grid, props, ...) -> (state_out, StepResult)
Standalone:  python -m steps.pressurization --state results/<run>/cycle_001/state_after_<previous>.csv
"""
from __future__ import annotations

from ._runner import run_step, step_cli


def _spec(cfg, name=None):
    for s in cfg.cycle:
        if s["type"] == "pressurization" and (name is None or s["name"] == name):
            return s
    raise KeyError("no 'pressurization' step in cfg.cycle")


def run_pressurization(state_in, cfg, grid, props, step=None, **kw):
    """Run the pressurization step (spec taken from cfg.cycle unless given) and return (state_out, result)."""
    return run_step(state_in, cfg, grid, props, step or _spec(cfg), **kw)


if __name__ == "__main__":
    step_cli("pressurization")
