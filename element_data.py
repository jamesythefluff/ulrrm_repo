from dataclasses import dataclass
from types import ModuleType

from quantum_defects.data import rubidium_defects, helium_defects
from scattering.data import (
    calcium_scattering, 
    helium_scattering, 
    rubidium_scattering, 
    krypton_scattering
)

AMU = 1.660539068e-27

@dataclass
class ElementData:
    mass: float
    ionisation_energy: float
    ground_state_polarisability: float | None
    defect_data: ModuleType | None  # Rydberg-Ritz parameters; None => use energy levels
    scattering_data: ModuleType

ELEMENTS = {
    "He": ElementData(
        mass = 4.00260325413 * AMU,
        ionisation_energy = 0.9035698802,
        ground_state_polarisability = None,
        defect_data = helium_defects,
        scattering_data = helium_scattering
    ),
    "Ca": ElementData(
        mass = 39.962590863 * AMU,
        ionisation_energy = 0.224654299701,
        ground_state_polarisability = None,
        defect_data = None,
        scattering_data = calcium_scattering
    ),
    "Rb": ElementData(
        mass = 84.9117897379 * AMU,
        ionisation_energy = 0.15350663,
        ground_state_polarisability = 319.2,
        defect_data = rubidium_defects,
        scattering_data = rubidium_scattering
    ),
    "Kr": ElementData(
        mass = 83.798 * AMU,
        ionisation_energy = None,
        ground_state_polarisability = None,
        defect_data = None,
        scattering_data = krypton_scattering
    )
}
