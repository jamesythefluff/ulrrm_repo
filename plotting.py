import matplotlib.pyplot as plt

def plot_potential(R_arr, potential_energy):
    plt.plot(R_arr, potential_energy)
    plt.grid()
    plt.xlabel("Radial distance")
    plt.ylabel("Potential energy")
    plt.title("He*-He ULRRM bonding potential")
    plt.show()