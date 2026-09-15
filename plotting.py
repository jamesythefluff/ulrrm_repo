import matplotlib.pyplot as plt
import numpy as np

from physics.scattering import Scattering

l_letters = {
    0: "s",
    1: "p",
    2: "d",
    3: "f",
    4: "g",
    5: "h",
    6: "i",
    7: "k",
    8: "l",
    9: "m",
    10: "n",
    11: "o",
    12: "q",
    13: "r",
    14: "t",
    15: "u",
    16: "v",
    17: "w",
    18: "x",
    19: "y",
    20: "z"
}

def plot_potential(R_arr, potential, labels, params):

    fig, ax = plt.subplots()

    for m in potential.keys():

        for channel in potential[m].keys():

            channel_potential = potential[m][channel]

            rydberg_element = params.rydberg_element
            perturber_element = params.perturber_element
            energy_units = params.energy_units

            for i, (n, l) in enumerate(labels[m][channel]):

                ax.plot(R_arr, channel_potential[:, i])

                # Put label at the start of each curve
                if l in l_letters:
                    ax.text(
                    R_arr[-1],
                    channel_potential[-1, i],
                    f"{n}{l_letters[l]}",
                    fontsize=8
                    )

    if params.xaxis is not None:
        ax.set_xlim(params.xaxis[0], params.xaxis[1])
    if params.yaxis is not None:
        ax.set_ylim(params.yaxis[0], params.yaxis[1])

    ax.grid()
    ax.set_xlabel(r"Radial distance ($a_0$)")
    ax.set_ylabel(f"Potential energy ({energy_units})")
    ax.set_title(f"{rydberg_element}*-{perturber_element} ULRRM bonding potential")

    plt.show()

def plot_phase_shifts(scattering_data, channels):
    """
    s- and p-wave phase shifts for each named scattering channel.
    """
    fig, (ax_s, ax_p) = plt.subplots(1, 2, figsize=(11, 4))

    for channel in channels:

        scattering = Scattering(scattering_data, channel)

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

def plot_scattering_lengths(scattering_data, channels):
    """
    Energy dependent s-wave scattering length and p-wave scattering volume,
    over the momentum range for which phase shift data exists.
    """
    fig, (ax_s, ax_p) = plt.subplots(1, 2, figsize=(11, 4))

    for channel in channels:

        scattering = Scattering(scattering_data, channel)

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
