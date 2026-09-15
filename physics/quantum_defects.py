"""
Creates a table of quantum defects using energy level data.
"""
import numpy as np
import pandas as pd

from data.energy_levels.energy_levels import get_energy

def calculate_quantum_defect(n, l, rydberg_spin, levels, ionisation_energy, rydberg_factor, rydberg_element):

    if l >= 4:
        return 0.0

    energy = get_energy(n, l, rydberg_spin, levels)

    if energy is None:
        print(f"Energy level for n = {n}, l = {l}, spin = {rydberg_spin} not found in data")
        return 0.0

    diff = ionisation_energy - energy

    quantum_defect = n - np.sqrt(rydberg_factor / (2*diff))

    return quantum_defect

def get_quantum_defects_Rb(n, l, s):

    if l >= 4:
        return 0.0

    quantum_defects = pd.read_csv("data/energy_levels/rubidium_quantum_defects.csv")

    row = quantum_defects[
        (quantum_defects["n"] == n) &
        (quantum_defects["l"] == l) &
        (quantum_defects["s"] == s)
    ]

    if row.empty:
        raise ValueError(f"Quantum defect not found for n={n}, l={l}, s={s}")

    return row.iloc[0]["quantum_defect"]