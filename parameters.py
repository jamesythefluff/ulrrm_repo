import numpy as np
from dataclasses import dataclass, field
from pathlib import Path

from quantum_defects.quantum_defects import calculate_quantum_defect
from element_data import ELEMENTS

ELECTRON_MASS = 9.1093837139e-31
CALCULATED_DIR = Path("potential_energy_curves") / "data" / "calculated_pecs"

@dataclass
class RydbergConfig:
    """
    Pure user-facing input.
    """
    rydberg_element : str = "He" # chemical symbol (str)
    perturber_element : str = "Kr"

    R_arr : np.ndarray = field(default_factory=lambda : np.linspace(50, 1600, 500))

    theta : float = 0.0 # Polar angle
    phi : float = 0.0 # Azimuthal angle

    n : int = 30
    n_range : tuple[int, int] = (-3, 6)

    l_min : int = 0
    l_max : int = 1000

    m_arr : list = field(default_factory=lambda : [0, 1])
    m : int | None = None

    reference_n : int = n
    reference_l : int | None = None

    # spin quantum number of valence Rydberg electrons
    # this modifies the quantum defect for atoms with more than
    # one valence electron such as helium
    rydberg_spin : int = 1 

    # there is only one singlet channel 1S0 1P1 ("singlet")
    # there are three singlet channels 3S1 3PJ where J=0,1,2 ("tripletJ=0" etc.)
    scattering_channels : tuple = ("singlet",)

    energy_units : str = "GHz" # GHz or Hartree

    fix_zero_energy_scattering_length : bool = True

    # Set false to turn off p-wave scattering
    p_wave_scattering : bool = True

    filename : str | None = None


    def __post_init__(self):
        if self.filename is None:
            self.filename = f"{self.rydberg_element}_{self.perturber_element}_n{n}.npz"
        self.pecs_file = CALCULATED_DIR / self.filename
    
    
    @classmethod
    def eiles_rb(cls):

        return cls(
            rydberg_element = "Rb",
            perturber_element = "Rb",

            R_arr = np.linspace(200, 1600, 500),

            n = 30,
            n_range = (-3, 6),

            l_min = 0,
            l_max = 100,

            scattering_channels = ("singlet", "tripletJ=2"),
        )

    @classmethod
    def he_kr(cls):

        return cls(
            rydberg_element = "He",
            perturber_element = "Kr",

            R_arr = np.linspace(100, 750, 500),

            n = 20,
            n_range = (-3, 6),

            l_min = 0,
            l_max = 1000,

            m_arr = [0],

            reference_n = 20,
            reference_l = 2,

            scattering_channels = ("singlet",),

            energy_units = "MHz",

            p_wave_scattering = True,

            filename = "He20d_Kr.npz"
        )


class RydbergSystem():

    def __init__(self, config):

        self.config = config

        # n_range holds offsets either side of config.n, not absolute n values.
        self.n_arr = np.arange(config.n + config.n_range[0], config.n + config.n_range[1] + 1)
        self.basis = [
            (n, l)
            for n in self.n_arr
            for l in range(config.l_min, min(n, config.l_max+1))
        ]

        self.classically_forbidden = False

        self.rydberg_element_data = ELEMENTS[config.rydberg_element]
        self.perturber_element_data = ELEMENTS[config.perturber_element]

        self.reduced_mass = (ELECTRON_MASS * self.rydberg_element_data.mass / 
                             (ELECTRON_MASS + self.rydberg_element_data.mass))
        self.rydberg_factor = self.reduced_mass / ELECTRON_MASS

        self.vibrational_reduced_mass = (
            self.rydberg_element_data.mass * self.perturber_element_data.mass
            / (self.rydberg_element_data.mass + self.perturber_element_data.mass)
        )

        self.quantum_defects = self.build_quantum_defects()

        self.reference_energy = self.reference_energy()

    def build_quantum_defects(self):

        return {
            (n, l): calculate_quantum_defect(
                n,
                l,
                self.config.rydberg_spin,
                self.rydberg_element_data.defect_data,
                self.config.rydberg_element
            )
            for n, l in self.basis
        }

    def reference_energy(self):

        config = self.config

        if config.reference_l:
            return (- self.rydberg_factor / 
             (2 * (config.reference_n - 
                   self.quantum_defects[(config.reference_n, config.reference_l)])**2)
            )
        else:
            return - self.rydberg_factor / (2 * config.reference_n**2)

