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

        def integrand(r):

            W = mp.whitw(
                self.n_eff,
                self.l + 0.5,
                2 * self.rf * r / self.n_eff
            )
            return abs(W)**2
        
        integral = mp.quad(integrand, [0, mp.inf])
    
        N = 1 / mp.sqrt(integral)
    
        return N

    def calculate_k(self, r, params):

        k_squared = (- self.rf**2 / self.n_eff**2 
                     + 2 * self.rf / r )
        
        if k_squared < 0:
            k = None
            if not params.classically_forbidden:
                print(f"r = {r} classically forbidden region reached")
                params.classically_forbidden = True
        else:
            k = np.sqrt(k_squared)
    
        return k

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

    def radial_gradient(self, r):

        Y_lm = self.angular()
        dR_dr = self.radial_derivative(r)
        
        grad = dR_dr * Y_lm

        return grad

def build_wavefunctions(params):

    return [Wavefunction(n,l,params) for (n,l) in params.basis]