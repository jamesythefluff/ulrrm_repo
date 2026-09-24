"""
Calculation of s- and p-wave scattering lengths using the
MERT expansion given in 
https://iopscience.iop.org/article/10.1088/0953-4075/22/9/010.
"""
import numpy as np

class MERTScattering:

    def __init__(self, data, channel):
        self.mert_params = data.mert_params
        self.channel = channel

        self.A = self.mert_params["A"]
        self.alpha_d = self.mert_params["alpha_d"]
        self.D = self.mert_params["D"]
        self.F = self.mert_params["F"]
        self.A1 = self.mert_params["A1"]
        self.alpha_q = self.mert_params["alpha_q"]
        self.H = self.mert_params["H"]

        self.a1 = np.pi / 15
        self.b1 = -83*np.pi / 23625
        self.c1 = -3*self.a1 / 7

    def tan_delta_0(self, k):
        k = np.asarray(k)

        tan_delta_0 = (-self.A*k*(1+4*self.alpha_d/3*k**2*np.log(k))
                  -np.pi*self.alpha_d/3*k**2
                  +self.D*k**3+self.F*k**4)

        return tan_delta_0

    def tan_delta_1(self, k):
        k = np.asarray(k)

        tan_delta_1 = (self.a1*self.alpha_d*k**2
                  - self.A1*k**3
                  + (self.b1*self.alpha_d**2+self.c1*self.alpha_q)*k**4
                  + self.H*k**5)

        return tan_delta_1

    def a_s(self, k):
        k = np.asarray(k)
        
        a_s = - self.tan_delta_0(k) / k
        
        return a_s


    def a_p(self, k):
        k = np.asarray(k)

        a_p = - self.tan_delta_1(k) / k**3

        return a_p