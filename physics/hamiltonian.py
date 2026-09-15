# IMPORTS
import numpy as np

def calculate_k(R, params):

    k_2 = (- params.rydberg_factor**2 / params.reference_level**2 
           + 2 * params.rydberg_factor / R)

    if k_2 < 0:
        raise ValueError(f"k_2 = {k_2} must not be negative")
    
    return np.sqrt(k_2)

def build_H0(wavefunctions, params):

    H0 = np.zeros(
            (len(params.basis), len(params.basis)),
            dtype=np.complex128
        )
    
    for i, wf in enumerate(wavefunctions):
        H0[i, i] = - params.rydberg_factor / (2 * wf.n_eff**2)

    return H0

def build_H(R,wavefunctions, H0, scattering, params):

    k = calculate_k(R, params)
    a_s = scattering.a_s(k)
    a_p_cubed = scattering.a_p(k)
    
    psi = np.array([wf.psi(R) for wf in wavefunctions], dtype=np.complex128)
    grad_psi = np.array([wf.gradient(R) for wf in wavefunctions], dtype=np.complex128)

    V_s = 2*np.pi * a_s * np.outer(np.conj(psi), psi)
    V_p = 6*np.pi * a_p_cubed * np.conj(grad_psi) @ grad_psi.T

    H = H0 + V_s + V_p

    return H