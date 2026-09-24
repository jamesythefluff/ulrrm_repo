"""
Uses fitted parameters to compute the vibrational states of Rydberg molecule
bound states in the parabolic well approx.
"""
import numpy as np

def harmonic_frequency(curvature, reduced_mass):
    omega = np.sqrt(curvature / reduced_mass)
    return omega

def vibrational_energies(curvature, reduced_mass, v_max):
    """
    Calculates the vibrational energies of states up to v_max 
    in Hartree relative to the minimum of the potential well.
    """
    omega = harmonic_frequency(curvature, reduced_mass)
    vibrational_energies = [omega*(1/2 + v) for v in np.arange(v_max)]
    return vibrational_energies