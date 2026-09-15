from dataclasses import dataclass
from types import ModuleType

from .scattering import helium_scattering, calcium_scattering, rubidium_scattering
from .quantum_defects import rubidium_defects

AMU = 1.660539068e-27

@dataclass
class ElementData:
    mass: float
    ionisation_energy: float
    energy_data: str
    defect_data: ModuleType | None  # Rydberg-Ritz parameters; None => use energy levels
    scattering_data: ModuleType

ELEMENTS = {
    "He": ElementData(
        mass = 4.00260325413 * AMU,
        ionisation_energy = 0.9035698802,
        energy_data = "helium_energy_levels.csv",
        defect_data = None,
        scattering_data = helium_scattering
    ),
    "Ca": ElementData(
        mass = 39.962590863 * AMU,
        ionisation_energy = 0.224654299701,
        energy_data = "calcium_energy_levels.csv",
        defect_data = None,
        scattering_data = calcium_scattering
    ),
    "Rb": ElementData(
        mass = 84.9117897379 * AMU,
        ionisation_energy = 0.15350663,
        energy_data = "rubidium_energy_levels.csv",
        defect_data = rubidium_defects,
        scattering_data = rubidium_scattering
    )
}
