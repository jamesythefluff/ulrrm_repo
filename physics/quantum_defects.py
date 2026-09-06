"""
Creates a table of quantum defects using energy level data.
"""
import numpy as np

from data.energy_levels import get_energy

def calculate_quantum_defect(n, l, spin, levels, ionisation_energy, rydberg_factor):

    energy = get_energy(n, l, spin, levels)

    if energy is None:
        return 0.0

    diff = ionisation_energy - energy

    quantum_defect = n - np.sqrt(rydberg_factor / (2*diff))

    return quantum_defect