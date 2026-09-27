"""
Plotting functions for visualising scattering phase shifts and lengths.

Contains separate functions that work on either Scattering or
MERTScattering objects.
"""
import numpy as np
import matplotlib.pyplot as plt

from scattering import Scattering
from scattering_noble_gases import MERTScattering
from physics.hamiltonian import calculate_k


def plot_phase_shifts(scattering_data, channels, system):
    """
    s- and p-wave phase shifts for each named scattering channel.
    """
    config = system.config
    alpha = system.perturber_element_data.ground_state_polarisability

    fig, (ax_s, ax_p) = plt.subplots(1, 2, figsize=(11, 4))

    for channel in channels:

        scattering = Scattering(
            scattering_data, 
            channel,
            alpha,
            config.fix_zero_energy_scattering_length
            )

        ax_s.plot(scattering.k_s_data, scattering.delta_0_data, label=channel)
        ax_p.plot(scattering.k_p_data, scattering.delta_1_data, label=channel)

    for ax, name in ((ax_s, "s-wave"), (ax_p, "p-wave")):
        ax.grid()
        ax.set_xlabel(r"Electron momentum $k$ (a.u.)")
        ax.set_ylabel("Phase shift (rad)")
        ax.set_title(f"{name} scattering phase shift")
        ax.legend()

    fig.tight_layout()
    plt.show()


def plot_scattering_lengths(scattering_data, channels, system):
    """
    Energy dependent s-wave scattering length and p-wave scattering volume,
    over the momentum range for which phase shift data exists.
    """
    config = system.config

    alpha = system.perturber_element_data.ground_state_polarisability

    fig, (ax_s, ax_p) = plt.subplots(1, 2, figsize=(11, 4))

    for channel in channels:

        scattering = Scattering(
            scattering_data, 
            channel, 
            alpha, 
            config.fix_zero_energy_scattering_length
        )

        a_s = [scattering.a_s(k) for k in scattering.k_s_data]
        a_p = [scattering.a_p(k) for k in scattering.k_p_data]

        ax_s.plot(scattering.k_s_data, a_s, label=channel)
        ax_p.plot(scattering.k_p_data, a_p, label=channel)

    ax_s.set_ylabel(r"$a_s(k)$ ($a_0$)")
    ax_s.set_title("s-wave scattering length")
    ax_p.set_ylabel(r"$a_p^3(k)$ ($a_0^3$)")
    ax_p.set_title("p-wave scattering volume")

    for ax in (ax_s, ax_p):
        ax.grid()
        ax.set_xlabel(r"Electron momentum $k$ (a.u.)")
        ax.legend()

    fig.tight_layout()
    plt.show()


def plot_mert_phase_shifts(system):
    """
    s- and p-wave phase shifts for each named scattering channel.
    """
    config = system.config
    channels = config.scattering_channels
    scattering_data = system.perturber_element_data.scattering_data

    k_arr = [calculate_k(R, system) for R in config.R_arr]

    fig, (ax_s, ax_p) = plt.subplots(1, 2, figsize=(11, 4))

    for channel in channels:

        scattering = MERTScattering(scattering_data, channel)

        tan_delta_0 = scattering.tan_delta_0(k_arr)
        tan_delta_1 = scattering.tan_delta_1(k_arr)

        delta_0 = np.arctan(tan_delta_0)
        delta_1 = np.arctan(tan_delta_1)

        ax_s.plot(k_arr, delta_0, label=channel)
        ax_p.plot(k_arr, delta_1, label=channel)

    for ax, name in ((ax_s, "s-wave"), (ax_p, "p-wave")):
        ax.grid()
        ax.set_xlabel(r"Electron momentum $k$ (a.u.)")
        ax.set_ylabel("Phase shift (rad)")
        ax.set_title(f"{name} scattering phase shift")
        ax.legend()

    fig.tight_layout()
    plt.show()


def plot_mert_scattering_lengths(system):
    """
    Energy dependent s-wave scattering length and p-wave scattering volume,
    over the momentum range for which phase shift data exists.
    """
    config = system.config
    channels = config.scattering_channels
    scattering_data = system.perturber_element_data.scattering_data

    k_arr = [calculate_k(R, system) for R in config.R_arr]

    fig, (ax_s, ax_p) = plt.subplots(1, 2, figsize=(11, 4))

    for channel in channels:

        scattering = MERTScattering(scattering_data, channel)

        a_s = scattering.a_s(k_arr)
        a_p = scattering.a_p(k_arr)

        ax_s.plot(k_arr, a_s, label=channel)
        ax_p.plot(k_arr, a_p, label=channel)

    ax_s.set_ylabel(r"$a_s(k)$ ($a_0$)")
    ax_s.set_title("s-wave scattering length")
    ax_p.set_ylabel(r"$a_p^3(k)$ ($a_0^3$)")
    ax_p.set_title("p-wave scattering volume")

    for ax in (ax_s, ax_p):
        ax.grid()
        ax.set_xlabel(r"Electron momentum $k$ (a.u.)")
        ax.legend()

    fig.tight_layout()
    plt.show()