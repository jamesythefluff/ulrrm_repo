import numpy as np
from scipy.optimize import linear_sum_assignment

from physics.wavefunctions import build_wavefunctions
from physics.hamiltonian import build_H0, build_H
from numerics.eigenvalues import diagonalise, label_eigenstates, adiabatic_reordering
from physics.scattering import Scattering

HARTREE_TO_GHZ = 6.5796839205e6

def calculate_pecs(params):

    pecs = {}
    labels = {}

    for m in params.m_arr:

        params.m = m

        wavefunctions = build_wavefunctions(params)

        H0 = build_H0(wavefunctions, params)

        pecs[m] = {}
        labels[m] = {}

        for channel in params.scattering_channels:

            channel_pecs, channel_labels = calculate_channel_pecs(channel, wavefunctions, H0, params)

            pecs[m][channel] = channel_pecs
            labels[m][channel] = channel_labels

    return pecs, labels

def calculate_channel_pecs(channel, wavefunctions, H0, params):

    channel_pecs = []
    
    large_R_eigenvectors = None

    prev_eigvecs = None

    data = params.perturber_element_data.scattering_data
    
    k_s_data = data.k_s_data[channel]
    k_p_data = data.k_p_data[channel]
    delta_0_data = data.delta_0_data[channel]
    delta_1_data = data.delta_1_data[channel]

    scattering = Scattering(k_s_data, delta_0_data, k_p_data, delta_1_data)

    for i, R in enumerate(reversed(params.R_arr)):

        H = build_H(R, wavefunctions, H0, scattering, params)

        assert np.allclose(H, H.conj().T) # Check Hermicity of Hamiltonian

        eigenvalues, eigenvectors = diagonalise(H)

        eigenvalues, eigenvectors, prev_eigvecs = adiabatic_reordering(prev_eigvecs, eigenvectors, eigenvalues)

        channel_pecs.append(eigenvalues)

        # Save eigenvectors at largest R for labelling purposes
        if large_R_eigenvectors is None:
            large_R_eigenvectors = eigenvectors

        print(f"{(i + 1) * 100 / len(params.R_arr):.1f} %")

    # R was calculated from large → small, so restore small → large
    channel_pecs = np.array(channel_pecs)
    channel_pecs = channel_pecs[::-1, :]

    # Shift PECs so reference level is zero
    channel_pecs = channel_pecs - params.reference_energy

    channel_labels = label_eigenstates(large_R_eigenvectors, params.basis)

    # Unit conversion to match literature (GHz)
    if params.energy_units == "GHz":
        channel_pecs *= HARTREE_TO_GHZ

    return channel_pecs, channel_labels