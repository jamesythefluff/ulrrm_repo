"""
Uses fitted parameters to compute the vibrational states of Rydberg molecule
bound states in the parabolic well approx.
"""
import numpy as np
import matplotlib.pyplot as plt

from analysis.pec_fitting import ParabolaFit, parabola

HA_TO_MHZ = 6.5796839205e9
HA_TO_GHZ = 6.5796839205e6
ELECTRON_MASS = 9.1093837139e-31


def harmonic_frequency(curvature, reduced_mass):
    omega = np.sqrt(2 * curvature / reduced_mass)
    return omega


def vibrational_energies(curvature, vibrational_reduced_mass, v_max, energy_units):
    """
    Calculates the vibrational energies of states up to v_max 
    relative to the minimum of the potential well.
    """
    if energy_units == "Ha":
        curvature /= 1
    elif energy_units == "MHz":
        curvature /= HA_TO_MHZ
    elif energy_units == "GHz":
        curvature /= HA_TO_GHZ
    else:
        raise ValueError("energy_units must be 'Ha', 'MHz' or 'GHz'")

    vib_red_mass_au = vibrational_reduced_mass / ELECTRON_MASS

    omega = harmonic_frequency(curvature, vib_red_mass_au)
    vibrational_energies = omega * (0.5 + np.arange(v_max + 1))

    if energy_units == "Ha":
        vibrational_energies *= 1
    elif energy_units == "MHz":
        vibrational_energies *= HA_TO_MHZ
    else:
        vibrational_energies *= HA_TO_GHZ

    return vibrational_energies


def plot_vibrational_states(para_fit : ParabolaFit, vib_states : np.ndarray) -> None:

    R_para = np.linspace(para_fit.R_min, para_fit.R_max, 500)
    E_para = parabola(R_para, para_fit.E0, para_fit.R0, para_fit.a)

    label = (
        f"E0 = {para_fit.E0} \n"
        f"R0 = {para_fit.R0} \n"
        f"a = {para_fit.a}"
    )

    plt.plot(para_fit.R, para_fit.pec, "--")
    plt.plot(para_fit.fitted_R, para_fit.fitted_pec, "o")
    plt.plot(R_para, E_para, label = label)

   # Plot vibrational energy levels
    for v, energy in enumerate(vib_states):
        plt.hlines(
            energy + para_fit.E0,
            para_fit.R_min,
            para_fit.R_max,
            linewidth=1
        )

    plt.grid()
    plt.legend()
    plt.show()