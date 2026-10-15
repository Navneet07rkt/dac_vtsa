"""Cycle steps.  Each module exposes run_<type>(state_in, ...) -> (state_out, StepResult)."""
from .adsorption import run_adsorption
from .blowdown import run_blowdown
from .heating import run_heating
from .cooling import run_cooling
from .pressurization import run_pressurization

RUNNERS = {"adsorption": run_adsorption, "blowdown": run_blowdown, "heating": run_heating,
           "cooling": run_cooling, "pressurization": run_pressurization}
