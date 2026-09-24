import numpy as np
from scipy.interpolate import CubicSpline

P_WAVE_CAP = 0.02

class Scattering:

    def __init__(self, data, channel, alpha, fix_zero_len):
        """
        Attributes
        --------------------
        data : ModuleType
            a module from data/scattering/, whose k_*_data and delta_*_data
                  are dicts keyed by scattering channel

        channel : str
            e.g. "singlet", "tripletJ=2"
        """
        self.alpha = alpha

        self.channel = channel

        self.k_s_data = data.k_s_data[channel]
        self.k_p_data = data.k_p_data[channel]
        self.delta_0_data = data.delta_0_data[channel]
        self.delta_1_data = data.delta_1_data[channel]

        self.delta_0_spline = CubicSpline(self.k_s_data, self.delta_0_data)
        self.delta_1_spline = CubicSpline(self.k_p_data, self.delta_1_data)

        if fix_zero_len:
            if channel == "singlet":
                self.fixed_a0 = data.zero_energy_scattering_lengths[channel]
            else:
                self.fixed_a0 = data.zero_energy_scattering_lengths["triplet"]
        else:
            self.fixed_a0 = None

        self.a0, self.a0_std, self.Re, self.Re_std = self.fit_mert()


    def delta_0(self, k):
        k = np.asarray(k)

        if k < self.k_s_data.min():
            raise ValueError(f"k = {k} below s-scattering data range")

        elif k > self.k_s_data.max():
            raise ValueError(f"k = {k} above s-scattering data range")
        
        return self.delta_0_spline(k)


    def delta_1(self, k):
        k = np.asarray(k)
        
        if k < self.k_p_data.min():
            raise ValueError(f"k = {k} below p-scattering data range")

        elif k > self.k_p_data.max():
            raise ValueError(f"k = {k} above p-scattering data range")

        return self.delta_1_spline(k)


    def a_s(self, k):

        a_s = - np.tan(self.delta_0(k)) / k
        
        return a_s


    def a_p(self, k):

        if k < P_WAVE_CAP:
            k = P_WAVE_CAP

        a_p = - np.tan(self.delta_1(k)) / k**3

        return a_p