from physics.pecs import calculate_pecs
from plotting import plot_potential, plot_phase_shifts, plot_scattering_lengths
from parameters import RydbergConfig, RydbergSystem

config = RydbergConfig.eiles_rb()
system = RydbergSystem(config)

potential, labels = calculate_pecs(system)
plot_potential(config.R_arr, potential, labels, system)