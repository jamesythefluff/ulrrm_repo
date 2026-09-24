"""

"""
import numpy as np


def rydberg_ritz(n, mu_0, mu_0_prime):
    """
    Two-term Rydberg-Ritz expansion of a fine-structure-resolved quantum defect,

        mu(n) = mu_0 + mu_0_prime / (n - mu_0)**2 .
    """
    return mu_0 + mu_0_prime / (n - mu_0)**2


def j_values_Rb(l):
    """
    The j of the fine structure doublet for orbital angular momentum l.
    l = 0 has only j = 1/2.
    """
    if l == 0:
        return [0.5]

    return [l - 0.5, l + 0.5]


def j_values_He(l, s):

    if s == 1:
        return [l]

    elif s == 3:
        if l == 0:
            return [1]
        else:
            return [l-1, l, l+1]
        
    else:
        raise Exception(f"S={s} is not a valid Helium total spin state.")


def spin_free_quantum_defect_He(n, l, s, defect_data):

    if l > defect_data.L_MAX:
        return 0.0

    energy = 0.0
    weight_total = 0.0

    for j in j_values_He(l, s):

        mu = rydberg_ritz(n, *defect_data.QUANTUM_DEFECTS[(l, s, j)])

        weight = 2*j + 1

        # Rydberg factor omitted: it cancels when the average is inverted below.
        energy += weight * (-0.5 / (n - mu)**2)
        weight_total += weight

    energy /= weight_total

    n_eff = np.sqrt(-0.5 / energy)

    return n - n_eff

    
def spin_free_quantum_defect_Rb(n, l, defect_data):
    """
    A single quantum defect per (n, l), with the fine structure collapsed.

    This model carries no Rydberg electron spin, so it needs one defect per l,
    but every measured alkali quantum defect is fine-structure resolved. The
    two are reconciled here by taking the (2j + 1) weighted centre of gravity
    of the fine structure multiplet in ENERGY,

        E_l = sum_j (2j + 1) E_n,l,j / sum_j (2j + 1),

    and converting that back into an effective defect. This is the standard
    definition of a fine-structure-free level position, and it is exact in the
    sense that it preserves the multiplet's centre of gravity; what it discards
    is the splitting itself.

    The splitting discarded is small but not uniformly negligible. At n = 30 in
    Rb the j dependence of the defect is worth roughly 4 GHz for p states,
    0.5 GHz for d states and under 10 MHz for f states. The p states sit far
    from the high-l manifold that dominates the trilobite and butterfly wells,
    so this is not expected to matter there, but it is an approximation of the
    model and should be quoted as one.

    For the fully spin-resolved treatment see M. T. Eiles and C. H. Greene,
    "Hamiltonian for the inclusion of spin effects in long-range Rydberg
    molecules", Phys. Rev. A 95, 042515 (2017).
    """
    if l > defect_data.L_MAX:
        # High-l states are treated as hydrogenic and exactly degenerate,
        # following C. H. Greene, A. S. Dickinson and H. R. Sadeghpour,
        # Phys. Rev. Lett. 85, 2458 (2000).
        return 0.0

    energy = 0.0
    weight_total = 0.0

    for j in j_values_Rb(l):

        mu = rydberg_ritz(n, *defect_data.QUANTUM_DEFECTS[(l, j - l)])

        weight = 2*j + 1

        # Rydberg factor omitted: it cancels when the average is inverted below.
        energy += weight * (-0.5 / (n - mu)**2)
        weight_total += weight

    energy /= weight_total

    n_eff = np.sqrt(-0.5 / energy)

    return n - n_eff


def calculate_quantum_defect(n, l, s, defect_data, rydberg_element):
    """
    Spin-free quantum defect from tabulated Rydberg-Ritz parameters.
    """
    if rydberg_element == "Rb":
        defect = spin_free_quantum_defect_Rb(n, l, defect_data)
        return defect
    elif rydberg_element == "He":
        defect = spin_free_quantum_defect_He(n, l, s, defect_data)
        return defect
    else:
        raise Exception(f"Unable to calculate defect for Rydberg element {rydberg_element}")