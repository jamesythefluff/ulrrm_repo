import numpy as np
from scipy.optimize import linear_sum_assignment

from physics.wavefunctions import build_wavefunctions
from physics.hamiltonian import build_H0, build_H
from numerics.eigenvalues import diagonalise, label_eigenstates

HARTREE_TO_GHZ = 6.5796839205e6

def calculate_pecs(params):

    pecs = {}
    labels = {}

    wavefunctions = build_wavefunctions(params)

    H0 = build_H0(wavefunctions, params)

    for channel in params.scattering_channels:

        channel_pecs = []

        large_R_eigenvectors = None

        prev_eigvecs = None

        for i, R in enumerate(reversed(params.R_arr)):

            H = build_H(R, wavefunctions, H0, channel, params)

            assert np.allclose(H, H.conj().T) # Check Hermicity of Hamiltonian

            eigenvalues, eigenvectors = diagonalise(H)

            # Reordering to track adiabatic states
            if prev_eigvecs is not None:
                overlap = np.abs(prev_eigvecs.conj().T @ eigenvectors)
                perm = linear_sum_assignment(-overlap)[1]   # scipy.optimize
                eigenvalues = eigenvalues[perm]
                eigenvectors = eigenvectors[:, perm]

            prev_eigvecs = eigenvectors

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

        if params.energy_units == "GHz":
            channel_pecs *= HARTREE_TO_GHZ

        pecs[channel] = channel_pecs
        labels[channel] = channel_labels

    return pecs, labels