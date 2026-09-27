import matplotlib.pyplot as plt
from matplotlib.pyplot import Axes
import numpy as np
import pandas as pd
from pathlib import Path
from dataclasses import dataclass, field

from potential_energy_curves.store_pecs import PECData, load_pecs

@dataclass
class PlotConfig:
    xlims : tuple | None = None
    ylims : tuple | None = None

    xscale : str = "linear"

    literature_pecs : list[str] = field(default_factory=list)

    @classmethod
    def eiles_rb(cls):
        return cls(
            xlims = (200, 1600),
            ylims = (-260, 10),

            literature_pecs = [
                "eiles_rb_singlet_trilobite.csv",
                "eiles_rb_butterfly.csv",
                "eiles_rb_triplet_trilobite_zoom.csv"
            ]
        )

    @classmethod
    def he_kr(cls):
        return cls(
            xlims = (0, 1000),
            ylims = (-60000, 5000),

            literature_pecs = ["deiglmayr_He_Kr_n20.csv"]
        )


LITERATURE_DIR = Path(__file__).parent / "data" / "literature_pecs"

# Cycled so overlaid literature curves stay tellable apart in the legend.
LITERATURE_STYLES = ["--", "-.", ":"]

l_letters = {
    0: "s", 1: "p", 2: "d", 3: "f", 4: "g", 5: "h", 6: "i", 7: "k", 8: "l",
    9: "m", 10: "n", 11: "o", 12: "q", 13: "r", 14: "t", 15: "u", 16: "v", 
    17: "w", 18: "x", 19: "y", 20: "z"
}


def scale_x(R : np.ndarray, x_scale : str) -> np.ndarray:
    """
    Eiles plots these curves against sqrt(R), which spreads out the small R
    region where the wells are. "linear" plots R itself.
    """
    if x_scale == "root":
        return np.sqrt(R)

    return R


def plot_pecs(pecs_file : Path, plot_config : PlotConfig | None = None) -> None:

    if not plot_config:
        plot_config = PlotConfig()

    pec_data = load_pecs(pecs_file)

    R_arr = pec_data.R
    potential = pec_data.energies
    labels = pec_data.labels

    fig, ax = plt.subplots()

    energy_units = pec_data.metadata["energy_units"]
    x_scale = plot_config.xscale

    # Digitised literature curves for comparison, plotted underneath.
    for i, filename in enumerate(plot_config.literature_pecs):
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

    if plot_config.xlims:
        ax.set_xlim(plot_config.xlims[0], plot_config.xlims[1])
    if plot_config.ylims:
        ax.set_ylim(plot_config.ylims[0], plot_config.ylims[1])

    ax.grid()
    if x_scale == "root":
        ax.set_xlabel(r"$\sqrt{R}$ ($\sqrt{a_0}$)")
    else:
        ax.set_xlabel(r"Radial distance $R$ ($a_0$)")

    rydberg_element = pec_data.metadata["rydberg_element"]
    perturber_element = pec_data.metadata["perturber_element"]

    ax.set_ylabel(f"Potential energy ({energy_units})")
    ax.set_title(
        f"{rydberg_element}*-{perturber_element} ULRRM bonding potential"
    )

    if getattr(plot_config, "literature_pecs", []):
        ax.legend(loc="lower right", fontsize=8)

    plt.show()


def plot_literature_pec(
        ax : Axes, 
        filename : str, 
        x_scale : str = "linear", 
        linestyle : str = "--"
        ) -> None:
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