from pathlib import Path

import potential_energy_curves.calculate_pecs as peccal
import potential_energy_curves.plot_pecs as pecplt
from potential_energy_curves.store_pecs import load_pecs
from parameters import RydbergConfig, RydbergSystem
import analysis.pec_fitting as pecfit
from analysis.vibrational_states import vibrational_energies, plot_vibrational_states

pec_file = Path("potential_energy_curves") / "data" / "calculated_pecs" / "He20d_Kr.npz"


def calculate_pecs():
    config = RydbergConfig.he_kr()
    system = RydbergSystem(config)

    peccal.calculate_pecs(system)


def plot_pecs():
    plot_config = pecplt.PlotConfig.he_kr()
    pecplt.plot_pecs(
        pec_file, 
        plot_config
    )


def fit_parabola():

    pec_data = load_pecs(pec_file)
    para_fit = pecfit.fit_parabola(pec_data, (20, 7), R_min=650, R_max=670)

    vib_lvls = vibrational_energies(
        para_fit.a, 
        pec_data.metadata["vibrational_reduced_mass"], 
        10, 
        pec_data.metadata["energy_units"]
    )

    plot_vibrational_states(para_fit, vib_lvls)


if __name__ == "__main_":
    calculate_pecs()

if __name__ == "__main__":
    plot_pecs()

if __name__ == "__main_":
    fit_parabola()