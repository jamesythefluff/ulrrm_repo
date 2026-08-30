import numpy as np
from scipy.interpolate import CubicSpline

"""
Scattering phase shift data from https://www.mdpi.com/2218-2004/9/4/82 used for interpolation.
"""

k_s_data = np.array([
    0.01, 0.02, 0.03, 0.04, 0.05,
    0.06, 0.07, 0.08, 0.10,
    0.20, 0.30, 0.40, 0.50,
    0.60, 0.70, 0.80, 0.90,
    1.00, 1.10
])

delta_0_data = np.array([
    -0.01181, -0.02386, -0.03612, -0.04858,
    -0.06119, -0.07396, -0.08686, -0.09987,
    -0.1262, -0.2614, -0.3978, -0.5319,
    -0.6609, -0.7835, -0.8990, -1.0069,
    -1.1071, -1.1998, -1.2854
])

k_p_data = np.array([
    0.10, 0.20, 0.30, 0.40, 0.50,
    0.60, 0.70, 0.80, 0.90, 1.00, 1.10
])

delta_1_data = np.array([
    0.00304, 0.01264, 0.02944, 0.05348,
    0.08380, 0.11856, 0.15535, 0.19175,
    0.22576, 0.25608, 0.28211
])

delta_0_interp = CubicSpline(k_s_data, delta_0_data)
delta_1_interp = CubicSpline(k_p_data, delta_1_data)

def calculate_delta_0(k):

    if not k_s_data.min() <= k <= k_s_data.max():
        raise ValueError(f"k = {k} outside s-scattering data range")
    
    return delta_0_interp(k)

def calculate_delta_1(k):

    if k < k_p_data.min():
        index_min = np.argmin(k_p_data)
        a_s = delta_1_data[index_min] / k_p_data[index_min]**3
        delta_1 = - a_s * k**3

    elif k > k_p_data.max():
        raise ValueError(f"k = {k} outside p-scattering data range")

    else:
        delta_1 = delta_1_interp(k)

    return delta_1