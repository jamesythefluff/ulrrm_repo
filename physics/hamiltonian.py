# IMPORTS
import numpy as np
from physics.scattering import Scattering

def build_H0(wavefunctions, params):

    H0 = np.zeros(
            (len(params.basis), len(params.basis)),
            dtype=np.complex128
        )
    
    for i, wf in enumerate(wavefunctions):
        H0[i, i] = - params.rydberg_factor / (2 * wf.n_eff**2)

    return H0

def build_H(R,wavefunctions, H0, params):

    scattering = Scattering(params.perturber_element_data.scattering_data)

    k_arr = np.array([wf.calculate_k(R,params) for wf in wavefunctions])
    a_s = np.array([scattering.a_s(k) for k in k_arr])
    a_p_cubed = np.array([scattering.a_p(k) for k in k_arr])
    
    psi = np.array([wf.psi(R) for wf in wavefunctions], dtype=np.complex128)
    grad_psi = np.array([wf.radial_gradient(R) for wf in wavefunctions], dtype=np.complex128)

    V_s = 2*np.pi * a_s * np.outer(np.conj(psi), psi)
    V_p = 6*np.pi * a_p_cubed * np.outer(np.conj(grad_psi), grad_psi)

    H = H0 + V_s + V_p

    return H