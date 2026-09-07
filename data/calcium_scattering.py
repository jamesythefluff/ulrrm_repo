import numpy as np
from scipy.interpolate import CubicSpline

E_eV_s_data = np.array([
    0.01, 0.03, 0.05, 0.07, 0.09, 0.10, 0.11, 0.13, 0.15, 0.17,
    0.20, 0.40, 0.60, 0.80, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00,
    2.20, 2.40, 2.60, 2.80, 3.00, 3.40, 3.80, 4.00, 4.40, 4.80,
    5.00
])

k_s_data = np.sqrt(2 * E_eV_s_data / 27.211386)

delta_0_data = np.array([
    0.2015, 0.1880, 0.1436, 0.0962, 0.0502, 0.0279, 0.0062, -0.0355,
    -0.0751, -0.1128, -0.1659, -0.4542, -0.6720, -0.8502, -1.0023,
    -1.1354, -1.2541, -1.3614, -1.4592, -1.5492, 1.5091, 1.4316,
    1.3591, 1.2911, 1.2271, 1.1094, 1.0034, 0.9541, 0.8620, 0.7775,
    0.737
])

E_eV_p_data = np.array([
    0.01, 0.03, 0.05, 0.07, 0.09, 0.10, 0.11, 0.13, 0.15, 0.17,
    0.20, 0.40, 0.60, 0.80, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00,
    2.20, 2.40, 2.60, 2.80, 3.00, 3.40, 3.80, 4.00, 4.40, 4.80,
    5.00
])

k_p_data = np.sqrt(2 * E_eV_p_data / 27.211386)

delta_1_data = np.array([
    0.0081, -0.0063, -0.0349, -0.0684, -0.1035, -0.1211, -0.1386,
    -0.1730, -0.2064, -0.2386, -0.2848, -0.5351, -0.7183, -0.8623,
    -0.9807, -1.0812, -1.1687, -1.2460, -1.3155, -1.3787, -1.4366,
    -1.4902, -1.5402, 1.5545, 1.5104, 1.4289, 1.3549, 1.3203,
    1.2550, 1.1946, 1.1659
])

delta_0_interp = CubicSpline(k_s_data, delta_0_data)
delta_1_interp = CubicSpline(k_p_data, delta_1_data)

def calculate_delta_0(k):

    """ if k < k_s_data.min():
        index_min = np.argmin(k_s_data)
        a_s = delta_0_data[index_min] / k_s_data[index_min]
        delta_0 = - a_s * k"""

    if k < k_s_data.min():
        raise ValueError(f"k = {k} < {k_s_data.min()} below s-scattering data range")

    elif k > k_s_data.max():
        raise ValueError(f"k = {k} > {k_s_data.max()} above s-scattering data range")

    else:
        delta_0 = delta_0_interp(k)

    return delta_0

def calculate_delta_1(k):

    """if k < k_p_data.min():
        index_min = np.argmin(k_p_data)
        a_p_cubed = delta_1_data[index_min] / k_p_data[index_min]**3
        delta_1 = - a_p_cubed * k**3"""

    if k < k_p_data.min():
        raise ValueError(f"k = {k} < {k_p_data.min()} below p-scattering data range")

    elif k > k_p_data.max():
        raise ValueError(f"k = {k} > {k_p_data.max()} above p-scattering data range")

    else:
        delta_1 = delta_1_interp(k)

    return delta_1
