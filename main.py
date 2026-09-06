from physics.pecs import calculate_pecs
from plotting import plot_potential, plot_phase_shifts
from parameters import Params

params = Params()

#plot_phase_shifts(params.perturber_element_data.scattering_data)

potential, labels = calculate_pecs(params)
plot_potential(params.R_arr, potential, labels, params)