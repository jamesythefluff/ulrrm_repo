import numpy as np
from scipy.interpolate import CubicSpline

class Scattering:

    def __init__(self, k_s_data, delta_0_data, k_p_data, delta_1_data):

        self.k_s_data = k_s_data
        self.k_p_data = k_p_data

        self.delta_0_spline = CubicSpline(k_s_data, delta_0_data)
        self.delta_1_spline = CubicSpline(k_p_data, delta_1_data)

    def delta_0(self, k):

        if k < self.k_s_data.min():
            raise ValueError(f"k = {k} below s-scattering data range")

        elif k > self.k_s_data.max():
            raise ValueError(f"k = {k} above s-scattering data range")
        
        return self.delta_0_spline(k)

    def delta_1(self, k):
        
        if k < self.k_p_data.min():
            raise ValueError(f"k = {k} below p-scattering data range")

        elif k > self.k_p_data.max():
            raise ValueError(f"k = {k} above p-scattering data range")

        return self.delta_1_spline(k)

    def a_s(self, k):

        a_s = - np.tan(self.delta_0(k)) / k

        return a_s

    def a_p(self, k):

        a_p = - np.tan(self.delta_1(k)) / k**3

        return a_p