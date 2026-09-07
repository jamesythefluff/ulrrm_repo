from physics.pecs import calculate_pecs
from plotting import plot_potential, plot_phase_shifts, plot_scattering_lengths
from parameters import Params

params = Params()

#plot_phase_shifts(params.perturber_element_data.scattering_data)

#plot_scattering_lengths(params.perturber_element_data.scattering_data)

potential, labels = calculate_pecs(params)
plot_potential(params.R_arr, potential, labels, params)