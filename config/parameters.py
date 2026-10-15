"""Expose the validated Config object to the rest of the code.

Usage
-----
>>> from config.parameters import load_config
>>> cfg = load_config()                                  # defaults.yaml
>>> cfg = load_config("config/cases/quick_test.yaml", overrides=["grid.Nr=4"])
>>> cfg.operation.T_heat
373.15

``CONFIG`` is a lazily-built default instance for interactive work.  Simulation
code never reads a global: every function receives the Config it should use,
so sensitivity runs with perturbed parameters cannot leak into each other.
"""
from __future__ import annotations

from typing import Any, Dict, Optional, Sequence

from . import schema
from .schema import Config, ConfigError  # noqa: F401  (re-export)

_DEFAULT: Optional[Config] = None


def load_config(user_yaml: Optional[str] = None, overrides: Optional[Sequence[str]] = None,
                override_dict: Optional[Dict[str, Any]] = None) -> Config:
    """Load defaults.yaml (+ optional user YAML and key=value overrides) and validate."""
    return schema.load(user_yaml, overrides, override_dict)


def default_config() -> Config:
    global _DEFAULT
    if _DEFAULT is None:
        _DEFAULT = schema.load()
    return _DEFAULT


def with_overrides(cfg: Config, override_dict: Dict[str, Any]) -> Config:
    """Return a new validated Config = cfg with dotted-key overrides applied."""
    raw = cfg.to_dict()
    raw.pop("derived", None)
    for k, v in override_dict.items():
        schema.set_dotted(raw, k, v)
    return schema.validate(raw)


class _LazyConfig:
    def __getattr__(self, item):
        return getattr(default_config(), item)


CONFIG = _LazyConfig()
