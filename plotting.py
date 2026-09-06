import matplotlib.pyplot as plt
import numpy as np

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

    states_to_plot = labels

    rydberg_element = params.rydberg_element
    perturber_element = params.perturber_element
    energy_units = params.energy_units

    fig, ax = plt.subplots()

    for i, (n, l) in enumerate(labels):
        print(
            f"{n}{l_letters[l]}:",
            "finite =", np.all(np.isfinite(potential[:, i])),
            "last =", potential[-1, i]
        )

    for i, (n, l) in enumerate(labels):

        if (n,l) not in states_to_plot:
            continue

        ax.plot(R_arr, potential[:, i])

        # Put label at the start of each curve
        if l in l_letters:
            ax.text(
            R_arr[-1],
            potential[-1, i],
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