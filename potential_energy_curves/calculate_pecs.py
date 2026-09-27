import numpy as np

from parameters import RydbergSystem
from physics.wavefunctions import build_wavefunctions
from physics.hamiltonian import build_H0, build_H
from numerics.eigenvalues import diagonalise, label_eigenstates, track_eigenstates
from scattering.scattering import Scattering
from scattering.scattering_noble_gases import MERTScattering
from potential_energy_curves.store_pecs import PECData, save_pecs

HARTREE_TO_GHZ = 6.5796839205e6
HARTREE_TO_MHZ = 6.5796839205e9


def calculate_pecs(system : RydbergSystem) -> None:

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

    pec_data = PECData(
        R = config.R_arr,
        energies = pecs,
        labels = labels,
        metadata={
            "rydberg_element" : config.rydberg_element,
            "perturber_element" : config.perturber_element,
            "n_arr" : system.n_arr,
            "m_arr" : config.m_arr,
            "rydberg_spin" : config.rydberg_spin,
            "reference_n" : config.reference_n,
            "reference_l" : config.reference_l,
            "reference_level" : system.reference_energy,
            "energy_units" : config.energy_units,
            "p_wave_scattering" : config.p_wave_scattering,
            "vibrational_reduced_mass" : system.vibrational_reduced_mass,
        }
    )

    print(config.pecs_file)

    save_pecs(pec_data, config.pecs_file)


def calculate_channel_pecs(
        channel : str, 
        wavefunctions : np.ndarray, 
        H0 : np.ndarray, 
        system : RydbergSystem
        ) -> tuple[np.ndarray, np.ndarray]:

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

    # Unit conversion to match literature
    if config.energy_units == "GHz":
        channel_pecs *= HARTREE_TO_GHZ
    elif config.energy_units == "MHz":
        channel_pecs *= HARTREE_TO_MHZ

    return channel_pecs, channel_labels