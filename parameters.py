import numpy as np

from data.energy_levels import load_energy_levels
from physics.quantum_defects import calculate_quantum_defect, get_quantum_defects_Rb
from data.element_data import ELEMENTS

ELECTRON_MASS = 9.1093837139e-31

class Params:

    def __init__(self):

        # SETUP PARAMETERS

        self.rydberg_element = "Rb" # chemical symbol (str)
        self.perturber_element = "Rb"
        self.energy_units = "GHz"

        self.R_arr = np.linspace(3, 1000, 500) # Ground perturber position relative to Rydberg ionic core

        self.phi = 0
        self.theta = 0

        self.reference_level = 30 # This hydrogenic level is set as the zero ref. for PECs

        # QUANTUM NUMBERS
        
        self.n = 30
        self.n_min = self.n - 1
        self.n_max = self.n + 4
        self.n_arr = np.arange(self.n_min, self.n_max)

        self.l = 8
        self.l_min = 0

        self.basis = [
            (n, l)
            for n in self.n_arr
            for l in range(self.l_min, n)
        ]

        self.m = 0
        self.rydberg_spin = 0 # spin quantum number of valence Rydberg electrons

        self.scattering_channels = ["singlet", "triplet"]

        # AUTOMATED SETUP STUFF

        self.classically_forbidden = False

        self.rydberg_element_data = ELEMENTS[self.rydberg_element]
        self.perturber_element_data = ELEMENTS[self.perturber_element]

        self.reduced_mass = ELECTRON_MASS * self.rydberg_element_data.mass / (ELECTRON_MASS + self.rydberg_element_data.mass)
        self.rydberg_factor = self.reduced_mass / ELECTRON_MASS

        if self.rydberg_element == "Rb":
            self.quantum_defects = self.build_quantum_defects_Rb()

        else:
            self.levels = load_energy_levels(self.rydberg_element_data.energy_data)
            self.quantum_defects = self.build_quantum_defects()

        self.reference_energy = - self.rydberg_factor / (2 * self.reference_level**2)

    def build_quantum_defects(self):

        quantum_defects = {}

        for n, l in self.basis:

            quantum_defects[(n, l)] = calculate_quantum_defect(
                 n,
                 l,
                 self.rydberg_spin,
                 self.levels,
                 self.rydberg_element_data.ionisation_energy,
                 self.rydberg_factor,
                 self.rydberg_element
            )

        return quantum_defects

    def build_quantum_defects_Rb(self):

        quantum_defects = {}

        for n, l in self.basis:

            quantum_defects[(n, l)] = get_quantum_defects_Rb(
                n,
                l,
                self.rydberg_spin
            )

        return quantum_defects


        