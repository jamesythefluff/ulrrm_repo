from calculations import calculate_potential
from plotting import plot_potential

R_arr, potential_energy = calculate_potential()
plot_potential(R_arr, potential_energy)