from dataclasses import dataclass
from types import ModuleType

from . import helium_scattering, calcium_scattering, rubidium_scattering_Fabrikant

AMU = 1.660539068e-27

@dataclass
class ElementData:
    mass: float
    ionisation_energy: float
    energy_data: str
    scattering_data: ModuleType

ELEMENTS = {
    "He": ElementData(
        mass = 4.00260325413 * AMU,
        ionisation_energy = 0.9035698802,
        energy_data = "helium_energy_levels.csv",
        scattering_data = helium_scattering
    ),
    "Ca": ElementData(
        mass = 39.962590863 * AMU,
        ionisation_energy = 0.224654299701,
        energy_data = "calcium_energy_levels.csv",
        scattering_data = calcium_scattering
    ),
    "Rb": ElementData(
        mass = 84.9117897379 * AMU,
        ionisation_energy = 0.15350663,
        energy_data = "rubidium_energy_levels.csv",
        scattering_data = rubidium_scattering_Fabrikant
    )
}