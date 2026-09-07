"""
Convert NIST energy level data to a formatted data frame.
"""
from pathlib import Path
import pandas as pd
import re 

DATA_DIR = Path(__file__).parent

L_MAP = {
    "S": 0,
    "P": 1,
    "D": 2,
    "F": 3,
    "G": 4,
    "H": 5,
    "I": 6,
    "K": 7,
    "L": 8,
    "M": 9,
    "N": 10,
    "O": 11,
    "Q": 12,
    "R": 13,
    "T": 14,
    "U": 15,
    "V": 16,
    "W": 17,
    "X": 18,
    "Y": 19,
    "Z": 20
}

def get_n_l(configuration):
    """
    Helper function for converting NIST's formalism to a more useful form.
    """
    # Find something like 10s, 12p, 15d, etc.
    match = re.search(r'(\d+)([a-z])$', configuration.lower())

    if match is None:
        return None, None

    n = int(match.group(1))
    l = L_MAP[match.group(2).upper()]

    return n, l

def load_energy_levels(energy_levels):
    """
    Read NIST data out of the input element's file, convert to a data
    frame, clean and format.
    """
    filename = energy_levels

    levels = pd.read_csv(DATA_DIR / "energy_levels" / filename)

    # Remove NIST's ="..." formatting
    for col in levels.columns:
        levels[col] = (
            levels[col]
            .astype(str)
            .str.replace('="', '', regex=False)
            .str.replace('"', '', regex=False)
        )

    # Make energy numerical
    levels["Level (Hartree)"] = pd.to_numeric(levels["Level (Hartree)"])

    levels[["n", "l"]] = levels["Configuration"].apply(
        lambda x: pd.Series(get_n_l(x))
    )

    return levels

def get_energy(n, l, rydberg_spin, levels):

    # Singlet = 1, triplet = 3
    multiplicity = 1 if rydberg_spin == 0 else 3

    # Find matching states
    matches = levels[
        (levels["n"] == n) &
        (levels["l"] == l) &
        (levels["Term"].str.startswith(str(multiplicity)))
    ]

    if len(matches) == 0:
        return None

    # Energy should be the same for different J if fine structure
    # isn't being considered, so take the first one.
    return matches.iloc[0]["Level (Hartree)"]

    