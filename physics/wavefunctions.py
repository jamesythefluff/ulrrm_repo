import numpy as np
import mpmath as mp
from scipy.special import sph_harm_y

class Wavefunction:

    def __init__(self,n,l,params):
        self.n = n
        self.l = l
        self.m = params.m

        self.quantum_defect = params.quantum_defects[(n,l)]
        self.n_eff = self.n - self.quantum_defect
        self.rf = params.rydberg_factor

        self.normalisation = self.calculate_normalisation()

        self.theta = params.theta
        self.phi = params.phi
        self.Y_lm = sph_harm_y(self.l, self.m, self.theta, self.phi)

    def calculate_normalisation(self):

        n_eff = self.n_eff
        l = self.l
        rf = self.rf

        N = np.sqrt(rf) / (
            n_eff
            * np.sqrt(
                mp.gamma(n_eff + l + 1)
                * mp.gamma(n_eff - l)
            )
        )

        return N

    def radial(self, r):

        W = mp.whitw(self.n_eff, self.l + 1/2, 2 * self.rf * r / self.n_eff)
        R_nl = self.normalisation * W / r

        return R_nl

    def angular(self):
        return self.Y_lm

    def psi(self, r):

        psi = self.radial(r) * self.angular()

        return psi

    def radial_derivative(self, r):

        n_eff = self.n_eff
        l = self.l
        z = 2 * self.rf * r / n_eff

        W = mp.whitw(n_eff, l + 0.5, z)
        W_plus = mp.whitw(n_eff + 1, l + 0.5, z)

        dW_dz = (0.5 - n_eff / z) * W - W_plus / z
        dW_dR = (2 * self.rf / n_eff) * dW_dz

        radial_derivative = self.normalisation * (dW_dR / r - W / r**2)

        return radial_derivative

    def gradient_m_0(self, r):

        Y_lm = self.angular()
        dR_dr = self.radial_derivative(r)

        grad = dR_dr * Y_lm

        return np.array([
            0.0,
            0.0,
            grad
        ])

    def gradient_m_mod_1(self, r):

        grad = []

        derivative = (self.radial(r) / (2*r) 
        * np.sqrt((2*self.l+1)*self.l*(self.l+1)/(4*np.pi))
        )

        return np.array([
            derivative,
            np.sign(self.m) * 1j * derivative,
            0.0
        ])


    def gradient(self, r):

        if self.m == 0:
            return self.gradient_m_0(r)
        elif abs(self.m) == 1:
            return self.gradient_m_mod_1(r)
        else:
            raise ValueError(f"Gradient not implemented for m = {self.m}")

def build_wavefunctions(params):

    return [Wavefunction(n,l,params) for (n,l) in params.basis]