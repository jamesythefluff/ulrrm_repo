"""
PECData class for storing calculated PECs and relevant metadata
for subsequent plotting and vibrational state analysis.
"""
import numpy as np
from dataclasses import dataclass
from pathlib import Path


@dataclass
class PECData:
    R : np.ndarray
    energies : dict
    labels : dict
    metadata : dict


def save_pecs(pec : PECData, pecs_file : Path) -> None:

    np.savez(
        pecs_file,
        R = pec.R,
        energies = pec.energies,
        labels = pec.labels,
        metadata = pec.metadata
    )


def load_pecs(pecs_file : Path) -> PECData:
    data = np.load(pecs_file, allow_pickle=True)

    return PECData(
        R = data["R"],
        energies = data["energies"].item(),
        labels = data["labels"].item(),
        metadata = data["metadata"].item()
    )