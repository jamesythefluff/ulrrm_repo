"""
Singlet and triplet phase shifts from https://journals.aps.org/pr/pdf/10.1103/PhysRev.179.78
"""
import numpy as np

# SINGLET SCATTERING CHANNEL 1S0 1P1

E_Ha_data = np.array([0.000125, 0.000250, 0.000375, 0.000500, 0.000750, 0.001000,
 0.001250, 0.001500, 0.002000, 0.002500, 0.003000, 0.003750,
 0.005000, 0.010000, 0.015000, 0.025000, 0.032500, 0.037500,
 0.042500, 0.050000, 0.062500, 0.075000, 0.125000, 0.250000,
 0.375000, 0.500000])

k_s_singlet_data = np.sqrt(2 * E_Ha_data)
k_p_singlet_data = k_s_singlet_data

delta_0_singlet_data = np.array([-0.035, -0.044, -0.055, -0.070, -0.105, -0.143, -0.182, -0.220,
 -0.291, -0.355, -0.412, -0.487, -0.595, -0.913, -1.126, -1.412,
 -1.558, -1.636, -1.703, -1.788, -1.902, -1.994, -2.273, -2.761,
 -3.113, -3.391])

delta_1_singlet_data = np.array([0.001, 0.014, 0.036, 0.063, 0.128, 0.203, 0.281, 0.361,
 0.517, 0.664, 0.796, 0.963, 1.160, 1.422, 1.470, 1.474,
 1.466, 1.457, 1.448, 1.433, 1.404, 1.373, 1.229, 0.891,
 0.629, 0.418])


# TRIPLET SCATTERING CHANNELS 3S1 3PJ

# J=0

k_s_triplet0_data = k_s_singlet_data
k_p_triplet0_data = k_s_singlet_data

delta_0_triplet0_data = np.array([
    0.590, 0.671, 0.689, 0.685, 0.656, 0.619, 0.581, 0.546,
    0.483, 0.428, 0.379, 0.312, 0.212, -0.087, -0.298, -0.603,
    -0.777, -0.876, -0.965, -1.085, -1.256, -1.401, -1.839,
    -2.513, -2.947, -3.270
])

delta_1_triplet0_data = np.array([
    -0.004, -0.000, 0.008, 0.018, 0.042, 0.064, 0.084, 0.100,
    0.123, 0.134, 0.139, 0.139, 0.128, 0.028, -0.096, -0.321,
    -0.464, -0.550, -0.629, -0.736, -0.893, -1.026, -1.425,
    -2.004, -2.358, -2.615
])

# J=1



k_s_data = {
    "singlet": k_s_singlet_data,
    "tripletJ=0": k_s_triplet0_data,
    "tripletJ=1": k_s_triplet1_data,
    "tripletJ=2": k_s_triplet2_data
}

k_p_data = {
    "singlet": k_p_singlet_data,
    "tripletJ=0": k_p_triplet0_data,
    "tripletJ=1": k_p_triplet1_data,
    "tripletJ=2": k_p_triplet2_data
}

delta_0_data = {
    "singlet": delta_0_singlet_data,
    "tripletJ=0": delta_0_triplet0_data,
    "tripletJ=1": delta_0_triplet1_data,
    "tripletJ=2": delta_0_triplet2_data
}

delta_1_data = {
    "singlet": delta_1_singlet_data,
    "tripletJ=0": delta_1_triplet0_data,
    "tripletJ=1": delta_1_triplet1_data,
    "tripletJ=2": delta_1_triplet2_data
}