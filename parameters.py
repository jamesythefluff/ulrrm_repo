import numpy as np
from dataclasses import dataclass, field

from quantum_defects.quantum_defects import calculate_quantum_defect
from data.element_data import ELEMENTS

ELECTRON_MASS = 9.1093837139e-31

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

    # spin quantum number of valence Rydberg electrons
    # this modifies the quantum defect for atoms with more than
    # one valence electron such as helium
    rydberg_spin : int = 1 


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

    fix_zero_energy_scattering_length : bool = True

    # Set false to turn off p-wave scattering
    p_wave_scattering : bool = True

    filename : str = f"{rydberg_element}_{perturber_element}_n{n}"

    
    @classmethod
    def eiles_rb(cls, x_scale="linear"):

        return cls(
            rydberg_element = "Rb",
            perturber_element = "Rb",

            R_arr = np.linspace(200, 1600, 500),

            n = 30,
            n_range = (-3, 6),

            l_min = 0,
            l_max = 100,

            scattering_channels = ("singlet", "tripletJ=2"),

            x_scale = x_scale,

            x_lims = [200, 1600],
            y_lims = [-260, 10],

            literature_pecs = [
                "eiles_rb_singlet_trilobite.csv",
                "eiles_rb_butterfly.csv",
                "eiles_rb_triplet_trilobite_zoom.csv"
            ]
        )

    @classmethod
    def he_kr(cls, x_scale="linear"):

        return cls(
            rydberg_element = "He",
            perturber_element = "Kr",

            R_arr = np.linspace(100, 1750, 500),

            n = 30,
            n_range = (0, 0),

            l_min = 1,
            l_max = 1,

            m_arr = [0],

            scattering_channels = ("singlet",),

            x_scale = x_scale,

            energy_units = "MHz",

            x_lims = [0, 2000],
            y_lims = [-40, 20],

            p_wave_scattering = False
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

        self.reference_level = config.n # This hydrogenic level is set as the zero ref. for PECs

        self.classically_forbidden = False

        self.rydberg_element_data = ELEMENTS[config.rydberg_element]
        self.perturber_element_data = ELEMENTS[config.perturber_element]

        self.reduced_mass = (ELECTRON_MASS * self.rydberg_element_data.mass / 
                             (ELECTRON_MASS + self.rydberg_element_data.mass))
        self.rydberg_factor = self.reduced_mass / ELECTRON_MASS

        self.quantum_defects = self.build_quantum_defects()

        self.reference_energy = - self.rydberg_factor / (2 * self.reference_level**2)


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
