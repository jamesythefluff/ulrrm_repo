from physics.potential_energy_curves import calculate_pecs
from plotting import (
    plot_potential, 
    plot_phase_shifts, 
    plot_scattering_lengths,
    plot_mert_scattering_lengths,
    plot_mert_phase_shifts
)
from parameters import RydbergConfig, RydbergSystem

config = RydbergConfig.he_kr()
system = RydbergSystem(config)

#plot_mert_phase_shifts(system)
#plot_mert_scattering_lengths(system)

potential, labels = calculate_pecs(system)
plot_potential(system)