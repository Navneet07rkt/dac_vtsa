"""Gas property backends selected by ``properties.backend``."""
from .ideal_gas_backend import IdealGasBackend, PropertyBackend, R_GAS, arden_buck  # noqa: F401


def species_for(cfg):
    return ["CO2", "H2O", "N2"] if cfg.components.water else ["CO2", "N2"]


def get_backend(cfg, species=None):
    species = species or species_for(cfg)
    if cfg.properties.backend == "coolprop":
        from .coolprop_backend import CoolPropBackend
        return CoolPropBackend(species, cfg.properties.cp_reference_T)
    return IdealGasBackend(species, cfg.properties.cp_reference_T)
