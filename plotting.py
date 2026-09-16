import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

from physics.scattering import Scattering

LITERATURE_DIR = Path(__file__).parent / "data" / "potential_energy_curves"

# Cycled so overlaid literature curves stay tellable apart in the legend.
LITERATURE_STYLES = ["--", "-.", ":"]

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

def scale_x(R, x_scale):
    """
    Eiles plots these curves against sqrt(R), which spreads out the small R
    region where the wells are. "linear" plots R itself.
    """
    if x_scale == "root":
        return np.sqrt(R)

    return R


def plot_potential(R_arr, potential, labels, system):

    config = system.config

    fig, ax = plt.subplots()

    energy_units = config.energy_units
    x_scale = getattr(config, "x_scale", "linear")

    # Digitised literature curves for comparison, plotted underneath.
    for i, filename in enumerate(getattr(config, "literature_pecs", [])):
        plot_literature_pec(
            ax, filename, x_scale, LITERATURE_STYLES[i % len(LITERATURE_STYLES)]
        )

    x = scale_x(np.asarray(R_arr, dtype=float), x_scale)

    for m in potential.keys():

        for channel in potential[m].keys():

            channel_potential = potential[m][channel]

            for i, (n, l) in enumerate(labels[m][channel]):

                ax.plot(x, channel_potential[:, i], linewidth=0.8)

                # Put label at the end of each curve
                if l in l_letters:
                    ax.text(
                    x[-1],
                    channel_potential[-1, i],
                    f"{n}{l_letters[l]}",
                    fontsize=8
                    )

    if config.x_lims:
        ax.set_xlim(config.x_lims[0], config.x_lims[1])
    if config.y_lims:
        ax.set_ylim(config.y_lims[0], config.y_lims[1])

    ax.grid()
    if x_scale == "root":
        ax.set_xlabel(r"$\sqrt{R}$ ($\sqrt{a_0}$)")
    else:
        ax.set_xlabel(r"Radial distance $R$ ($a_0$)")
    ax.set_ylabel(f"Potential energy ({energy_units})")
    ax.set_title(
        f"{config.rydberg_element}*-{config.perturber_element} ULRRM bonding potential"
    )

    if getattr(config, "literature_pecs", []):
        ax.legend(loc="lower right", fontsize=8)

    plt.show()


def plot_literature_pec(ax, filename, x_scale="linear", linestyle="--"):
    """
    Overlay a digitised literature curve. The CSVs are headerless, column 0 the
    internuclear distance in a0 and column 1 the energy in GHz, so they must be
    plotted against the same axis scaling as the calculated curves.
    """
    pecs = pd.read_csv(LITERATURE_DIR / filename, header=None)

    # WebPlotDigitizer returns points in extraction order, not in order of
    # increasing R, so joining them up as read produces a zigzag rather than a
    # curve. Point order carries no meaning here, so sort on R.
    pecs = pecs.sort_values(by=0)

    distance = scale_x(pecs.iloc[:, 0].to_numpy(dtype=float), x_scale)
    energy = pecs.iloc[:, 1].to_numpy(dtype=float)

    ax.plot(
        distance,
        energy,
        color="k",
        linestyle=linestyle,
        linewidth=1.2,
        label=Path(filename).stem,
    )


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
