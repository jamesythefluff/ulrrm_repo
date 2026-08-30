"""
Implementation of Alisee's method to calculate ULRRM potentials for He*-He.

The Fermi-Omont pseudopotential (zero-range approx.) is used with scattering
phase shifts interpolated from literature data. Currently quantum defects are
approximated as negligable, but they should be added in later.

The Hamiltonian matrix is diagonalised in the basis of Rydberg eigenstates,
with states from nearby n manifolds included to account for
eigenstate mixing due to the perturbing atom, allowing for convergence.

Starting from large R, the overlap of the resulting eigenvectors with the 
initial Rydberg eigenvector is computed, and the eigenvector with the largest 
overlap is chosen.

R is decreased and the process is repeated, taking the overlap with the 
previously chosen eigenvector, in this way the state is "tracked" inwards.

The resulting ULRRM potential is labelled by the initial Rydberg eigenstate 
to which it converges in the long-range limit of R.
"""
# IMPORTS
import numpy as np
import matplotlib.pyplot as plt
from parameters import Params
import mpmath as mp
from scipy.special import sph_harm_y
import data.scattering_data_He as He
import time

class Wavefunction:

    def __init__(self,n,l,params):
        self.n = n
        self.l = l
        self.m = params.m

        self.quantum_defect = params.quantum_defect
        self.n_eff = self.n - self.quantum_defect

        self.theta = params.theta
        self.phi = params.phi

    @property
    def normalisation(self):

        def integrand(r):

            W = mp.whitw(
                self.n_eff,
                self.l + 0.5,
                2*r/self.n_eff
            )
            return abs(W)**2
        
        integral = mp.quad(integrand, [0, mp.inf])
    
        N = 1 / mp.sqrt(integral)
    
        return N

    def radial(self,r):

        W = mp.whitw(self.n_eff, self.l + 1/2, 2*r / self.n_eff)
        R_nl = self.normalisation * W / r

        return R_nl

    def angular(self):

        Y_lm = sph_harm_y(self.l, self.m, 0, 0)

        return Y_lm

    def psi(self,r):

        psi = self.radial(r) * self.angular()

        return psi

    def radial_derivative(self, r):

        n_eff = self.n_eff
        l = self.l
        z = 2 * r / n_eff

        W = mp.whitw(n_eff, l + 0.5, z)
        W_plus = mp.whitw(n_eff + 1, l + 0.5, z)

        dW_dz = (0.5 - n_eff / z) * W - W_plus / z
        dW_dR = (2 / n_eff) * dW_dz

        radial_derivative = self.normalisation * (dW_dR / r - W / r**2)

        return radial_derivative

    def gradient(self,r):

        Y_lm = self.angular()
        dR_dr = self.radial_derivative(r)
        
        grad = dR_dr * Y_lm

        return grad


def calculate_a_s(k):
    """
    Calculate the s-wave scattering length a_s for a given momentum k using:

        a_s = - tan(delta_0) / k

    """
    a_s = - np.tan(He.calculate_delta_0(k)) / k

    return a_s

def calculate_a_p_cubed(k):
    """
    Calculate the p-wave scattering volume a_p**3 for a given momentum k using:

        a_p**3 = - tan(delta_1) / k**3

    """
    a_p_cubed = - np.tan(He.calculate_delta_1(k)) / k**3

    return a_p_cubed

def calculate_k(r,params):
    """
    Calculates the electron semi-classical momentum at radial position r by
    approximating the total energy by the nth Hydrogen energy level.

    Arguments
    ----------------------------------------------------------------------
    r : float
        Radial position of the electron relative to the parent ion.

    params : Params
        Contains the calculation parameters such as the principal quantum 
        number.

    Returns
    ----------------------------------------------------------------------
    k : float
        Momentum at the given position r.
    """
    k_squared = -1 / params.n**2 + 2 / r

    if k_squared < 0:
        k = 0
        if not params.classically_forbidden:
            print(f"r = {r} classically forbidden region reached")
            params.classically_forbidden = True
    else:
        k = np.sqrt(k_squared)

    return k

def calculate_H_matrix(R,wavefunctions,params):
    """
    Calculates the matrix element of Rydberg states with quantum numbers
    n1, l1 and n2, l2 when the ground state He is at position R relative to 
    the Rydberg ion.

    Arguments
    ----------------------------------------------------------------------
    R : float
        Radial position of ground state He relative to the Rydberg He ion.
    """
    k = calculate_k(R,params)
    
    psi = np.array([wf.psi(R) for wf in wavefunctions], dtype=np.complex128)
    grad_psi = np.array([wf.gradient(R) for wf in wavefunctions], dtype=np.complex128)

    V_s = 2*np.pi * calculate_a_s(k) * np.outer(np.conj(psi), psi)
    V_p = 6*np.pi * calculate_a_p_cubed(k) * np.outer(np.conj(grad_psi), grad_psi)

    # unperturbed hamiltonian
    H_0 = np.zeros((len(params.basis),len(params.basis)), dtype=np.complex128)

    for i, (n,l) in enumerate(params.basis):
        H_0[i,i] = - 1 / (2*n**2)

    H_matrix = H_0 + V_s + V_p

    return H_matrix

def calculate_potential():

    params = Params()

    basis = params.basis
    R_arr = params.R_arr

    target_index = basis.index((params.n,params.l))

    potential_energy = []

    previous_eigenvector = None

    wavefunctions = [Wavefunction(n,l,params) for (n,l) in basis]

    progress = 0

    for R in reversed(R_arr):

        H_matrix = calculate_H_matrix(R,wavefunctions,params)

        eigenvalues, eigenvectors = np.linalg.eigh(H_matrix)

        if previous_eigenvector is None:

            # Identify target from params.n, params.l
            target_index = basis.index((params.n, params.l))
            overlaps = np.abs(eigenvectors[target_index, :])**2
            eigenvalue_index = np.argmax(overlaps)

        else:

            # Follow the same adiabatic state
            overlaps = np.abs(
                previous_eigenvector.conj() @ eigenvectors
            )**2
            eigenvalue_index = np.argmax(overlaps)

        selected_energy = eigenvalues[eigenvalue_index]
        potential_energy.append(selected_energy)
        previous_eigenvector = eigenvectors[:, eigenvalue_index]

        progress += 1

        print(f"{progress * 100 / len(R_arr)} %")

    potential_energy = potential_energy[::-1]

    return R_arr, potential_energy