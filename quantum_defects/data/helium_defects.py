"""
The first 2 Helium Rydberg-Ritz expansion parameters
mu_0, mu_2 keyed by:

    (L, S, J)

where L is the angular momentum of the Rydberg electron,
S is the total spin of the valence electrons (singlet
or triplet), and J is the fine structure quantum number.

J dependence will be removed using a centroid (2J+1)
average.

Data from 
https://iopscience.iop.org/article/10.1238/Physica.Topical.083a00083
table VII
"""

QUANTUM_DEFECTS = {
    (0, 1, 0): (0.13972, 0.027836),
    (0, 3, 1): (0.29666, 0.038297),

    (1, 1, 1): (-0.012142, 0.007519),
    (1, 3, 0): (0.068328, -0.018642),
    (1, 3, 1): (0.068358, -0.018630),
    (1, 3, 2): (0.068360, -0.018629),

    (2, 1, 2): (0.002113, -0.003090),
    (2, 3, 1): (0.002886, -0.006358),
    (2, 3, 2): (0.002891, -0.006357),
    (2, 3, 3): (0.002891, -0.006358),

    (3, 1, 3): (0.000440, -0.001689),
    (3, 3, 2): (0.000444, -0.001739),
    (3, 3, 3): (0.000449, -0.001727),
    (3, 3, 4): (0.000447, -0.001739)
}

L_MAX = 3