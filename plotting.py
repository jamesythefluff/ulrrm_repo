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

states_to_plot = [
    (29, 8),  
    (30, 8),  
    (31, 2),  
    (32, 1),  
    (33, 0)
]

def plot_potential(R_arr, potential, labels, params):

    fig, ax = plt.subplots()

    for channel in potential.keys():

        channel_potential = potential[channel]

        states_to_plot = labels[channel]

        rydberg_element = params.rydberg_element
        perturber_element = params.perturber_element
        energy_units = params.energy_units

        for i, (n, l) in enumerate(labels[channel]):

            if (n,l) not in states_to_plot:
                continue

            ax.plot(R_arr, channel_potential[:, i])

            # Put label at the start of each curve
            if l in l_letters:
                ax.text(
                R_arr[-1],
                channel_potential[-1, i],
                f"{n}{l_letters[l]}",
                fontsize=8
                )

    ax.grid()
    ax.set_xlabel(r"Radial distance ($a_0$)")
    ax.set_ylabel(f"Potential energy ({energy_units})")
    ax.set_title(f"{rydberg_element}*-{perturber_element} ULRRM bonding potential")

    plt.show()

def plot_phase_shifts(scattering_data):

    k_s = scattering_data.k_s_data 
    delta_0 = scattering_data.delta_0_data

    plt.plot(k_s, delta_0)
    plt.grid()
    plt.xlabel("Momentum k")
    plt.ylabel("Phase Shift (Rad)")
    plt.title("s-wave scattering phase shift")
    plt.show()

    k_p = scattering_data.k_p_data
    delta_1 = scattering_data.delta_1_data

    plt.plot(k_p, delta_1)
    plt.grid()
    plt.xlabel("Momentum k")
    plt.ylabel("Phase Shift (Rad)")
    plt.title("p-wave scattering phase shift")
    plt.show()

def plot_scattering_lengths(scattering_data):

    scattering = Scattering(scattering_data)

    k_s_arr = scattering_data.k_s_data
    k_p_arr = scattering_data.k_p_data

    a_s_arr = [scattering.a_s(k_s) for k_s in k_s_arr]
    a_p_arr = [scattering.a_p(k_p) for k_p in k_p_arr]

    plt.plot(k_s_arr, a_s_arr)
    plt.grid() 
    plt.show()

    plt.plot(k_p_arr, a_p_arr)
    plt.grid()
    plt.show()