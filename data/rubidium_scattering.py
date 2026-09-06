"""
Singlet phase shifts from https://journals.aps.org/pr/pdf/10.1103/PhysRev.179.78
"""

import numpy as np
from scipy.interpolate import CubicSpline

k_s_data = np.array([0.000125, 0.000250, 0.000375, 0.000500, 0.000750, 0.001000,
 0.001250, 0.001500, 0.002000, 0.002500, 0.003000, 0.003750,
 0.005000, 0.010000, 0.015000, 0.025000, 0.032500, 0.037500,
 0.042500, 0.050000, 0.062500, 0.075000, 0.125000, 0.250000,
 0.375000, 0.500000])

delta_0_data = np.array([-0.035, -0.044, -0.055, -0.070, -0.105, -0.143, -0.182, -0.220,
 -0.291, -0.355, -0.412, -0.487, -0.595, -0.913, -1.126, -1.412,
 -1.558, -1.636, -1.703, -1.788, -1.902, -1.994, -2.273, -2.761,
 -3.113, -3.391])

k_p_data = np.array([0.000125, 0.000250, 0.000375, 0.000500, 0.000750, 0.001000,
 0.001250, 0.001500, 0.002000, 0.002500, 0.003000, 0.003750,
 0.005000, 0.010000, 0.015000, 0.025000, 0.032500, 0.037500,
 0.042500, 0.050000, 0.062500, 0.075000, 0.125000, 0.250000,
 0.375000, 0.500000])

delta_1_data = np.array([0.001, 0.014, 0.036, 0.063, 0.128, 0.203, 0.281, 0.361,
 0.517, 0.664, 0.796, 0.963, 1.160, 1.422, 1.470, 1.474,
 1.466, 1.457, 1.448, 1.433, 1.404, 1.373, 1.229, 0.891,
 0.629, 0.418])

delta_0_interp = CubicSpline(k_s_data, delta_0_data)
delta_1_interp = CubicSpline(k_p_data, delta_1_data)

def calculate_delta_0(k):

    if k < k_s_data.min():
        index_min = np.argmin(k_s_data)
        a_s = delta_0_data[index_min] / k_s_data[index_min]
        delta_0 = - a_s * k

    elif k > k_s_data.max():
        raise ValueError(f"k = {k} outside s-scattering data range")

    else:
        delta_0 = delta_0_interp(k)

    return delta_0

def calculate_delta_1(k):

    if k < k_p_data.min():
        index_min = np.argmin(k_p_data)
        a_p_cubed = delta_1_data[index_min] / k_p_data[index_min]**3
        delta_1 = - a_p_cubed * k**3

    elif k > k_p_data.max():
        raise ValueError(f"k = {k} outside p-scattering data range")

    else:
        delta_1 = delta_1_interp(k)

    return delta_1