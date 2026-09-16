import numpy as np
from dataclasses import dataclass, field

from data.energy_levels.energy_levels import load_energy_levels
from physics.quantum_defects import (
    calculate_quantum_defect_from_levels,
    calculate_quantum_defect_rydberg_ritz,
)
from data.element_data import ELEMENTS

ELECTRON_MASS = 9.1093837139e-31

@dataclass
class RydbergConfig:
    """
    Pure user-facing input.
    """
    rydberg_element : str = "Rb" # chemical symbol (str)
    perturber_element : str = "Rb"

    R_arr : np.ndarray = field(default_factory=lambda : np.linspace(200, 1600, 500))

    theta : float = 0.0 # Polar angle
    phi : float = 0.0 # Azimuthal angle

    n : int = 30
    n_range : tuple[int, int] = (-1, 4)

    l_min : int = 0
    l_max : int = 1000

    m_arr : list = field(default_factory=lambda : [-1, 0, 1])
    m : int | None = None

    rydberg_spin : float = 0.0 # spin quantum number of valence Rydberg electrons

    # there is only one singlet channel 1S0 1P1 ("singlet")
    # there are three singlet channels 3S1 3PJ where J=0,1,2 ("tripletJ=0" etc.)
    scattering_channels : tuple = ("singlet",)

    # "linear" plots R, "root" plots sqrt(R) as Eiles does. Note that
    # xaxis limits are in whatever units the chosen scale produces.
    x_scale : str = "linear"
    energy_units : str = "GHz" # GHz or Hartree

    x_lims : list[float] = field(default_factory=list)
    y_lims : list[float] = field(default_factory=list)

    # Digitised literature curves in data/potential_energy_curves/ to
    # overlay on the PECs for comparison.
    literature_pecs : list[str] = field(default_factory=list)

    
    @classmethod
    def eiles_rb(cls, x_scale="linear"):

        return cls(
            rydberg_element = "Rb",
            perturber_element = "Rb",

            R_arr = np.linspace(200, 1600, 500),

            n = 30,
            n_range = (-1, 4),

            l_min = 0,
            l_max = 34,

            scattering_channels = ("singlet", "tripletJ=2"),

            x_scale = x_scale,

            x_lims = [200, 1600],
            y_lims = [-260, 10],

            literature_pecs = [
                "eiles_rb_trilobite_lower.csv",
                "eiles_rb_trilobite_upper.csv"
            ]
        )


class RydbergSystem():

    def __init__(self, config):

        self.config = config

        # n_range holds offsets either side of config.n, not absolute n values.
        self.n_arr = np.arange(config.n + config.n_range[0], config.n + config.n_range[1])
        self.basis = [
            (n, l)
            for n in self.n_arr
            for l in range(config.l_min, min(n, config.l_max+1))
        ]

        self.reference_level = config.n # This hydrogenic level is set as the zero ref. for PECs

        self.classically_forbidden = False

        self.rydberg_element_data = ELEMENTS[config.rydberg_element]
        self.perturber_element_data = ELEMENTS[config.perturber_element]

        self.reduced_mass = (ELECTRON_MASS * self.rydberg_element_data.mass / 
                             (ELECTRON_MASS + self.rydberg_element_data.mass))
        self.rydberg_factor = self.reduced_mass / ELECTRON_MASS

        # Prefer fitted Rydberg-Ritz parameters where we have them, and fall
        # back to inverting measured energy levels where we do not.
        if self.rydberg_element_data.defect_data is not None:
            self.quantum_defects = self.build_quantum_defects_rydberg_ritz()

        else:
            self.levels = load_energy_levels(self.rydberg_element_data.energy_data)
            self.quantum_defects = self.build_quantum_defects_from_levels()

        self.reference_energy = - self.rydberg_factor / (2 * self.reference_level**2)

    def build_quantum_defects_from_levels(self):

        return {
            (n, l): calculate_quantum_defect_from_levels(
                n,
                l,
                self.config.rydberg_spin,
                self.levels,
                self.rydberg_element_data.ionisation_energy,
                self.rydberg_factor,
            )
            for n, l in self.basis
        }

    def build_quantum_defects_rydberg_ritz(self):

        return {
            (n, l): calculate_quantum_defect_rydberg_ritz(
                n,
                l,
                self.rydberg_element_data.defect_data,
            )
            for n, l in self.basis
        }
