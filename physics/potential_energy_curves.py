import numpy as np
from pathlib import Path

from physics.wavefunctions import build_wavefunctions
from physics.hamiltonian import build_H0, build_H
from numerics.eigenvalues import diagonalise, label_eigenstates, track_eigenstates
from scattering.scattering import Scattering
from scattering.scattering_noble_gases import MERTScattering

HARTREE_TO_GHZ = 6.5796839205e6
HARTREE_TO_MHZ = 6.5796839205e9

DIR = Path(__file__).parent / "data" / "calculated_pecs"


def calculate_pecs(system):

    config = system.config

    pecs = {}
    labels = {}

    for m in config.m_arr:

        config.m = m

        wavefunctions = build_wavefunctions(system)

        H0 = build_H0(wavefunctions, system)

        pecs[m] = {}
        labels[m] = {}

        for channel in config.scattering_channels:

            channel_pecs, channel_labels = calculate_channel_pecs(channel, wavefunctions, H0, system)

            pecs[m][channel] = channel_pecs
            labels[m][channel] = channel_labels

    save_pecs(config.filename, config.R_arr, pecs, labels)

    return pecs, labels


def calculate_channel_pecs(channel, wavefunctions, H0, system):

    config = system.config

    large_R_eigenvectors = None

    eigenvalues_list = []
    eigenvectors_list = []

    data = system.perturber_element_data.scattering_data

    if config.perturber_element == "Kr":
        scattering = MERTScattering(data, channel)
    else:
        alpha = system.perturber_element_data.ground_state_polarisability
        fix_zero_len = config.fix_zero_energy_scattering_length
        scattering = Scattering(data, channel, alpha, fix_zero_len)

    for i, R in enumerate(reversed(config.R_arr)):

        H = build_H(R, wavefunctions, H0, scattering, system)

        assert np.allclose(H, H.conj().T) # Check Hermicity of Hamiltonian

        eigenvalues, eigenvectors = diagonalise(H)

        eigenvalues_list.append(eigenvalues)
        eigenvectors_list.append(eigenvectors)

        # Save eigenvectors at largest R for labelling purposes
        if large_R_eigenvectors is None:
            large_R_eigenvectors = eigenvectors

        print(f"{(i + 1) * 100 / len(config.R_arr):.1f} %")

    # Match each R step to the one before it so curves follow adiabatic states
    # rather than the energy ordering eigh happens to return.
    channel_pecs, _ = track_eigenstates(eigenvalues_list, eigenvectors_list)

    # R was calculated from large → small, so restore small → large
    channel_pecs = channel_pecs[::-1, :]

    # Shift PECs so reference level is zero
    channel_pecs = channel_pecs - system.reference_energy

    channel_labels = label_eigenstates(large_R_eigenvectors, system.basis)

    # Unit conversion to match literature (GHz)
    if config.energy_units == "GHz":
        channel_pecs *= HARTREE_TO_GHZ
    elif config.energy_units == "MHz":
        channel_pecs *= HARTREE_TO_MHZ

    return channel_pecs, channel_labels


def save_pecs(filename, R_arr, pecs, labels):

    np.savez(
        DIR / filename,
        R = R_arr,
        energies = pecs,
        labels = labels
    )