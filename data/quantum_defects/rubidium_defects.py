"""
Rydberg-Ritz quantum defect parameters for rubidium.

Values are Table 1 of

    M. T. Eiles, "Trilobites, butterflies, and other exotic specimens of
    long-range Rydberg molecules", J. Phys. B: At. Mol. Opt. Phys. 52,
    113001 (2019).

Eiles in turn takes them from the millimetre-wave measurements of

    W. Li, I. Mourachko, M. W. Noel and T. F. Gallagher, "Millimeter-wave
    spectroscopy of cold Rb Rydberg atoms in a magneto-optical trap: Quantum
    defects of the ns, np, and nd series", Phys. Rev. A 67, 052502 (2003)

for the s, p and d series, and

    J. Han, Y. Jamil, D. V. L. Norum, P. J. Tanner and T. F. Gallagher,
    "Rb nf quantum defects from millimeter-wave spectroscopy of cold 85Rb
    Rydberg atoms", Phys. Rev. A 74, 054502 (2006)

for the f series.

Each entry is keyed by (l, j - l) and holds (mu_0, mu_0_prime), the first two
coefficients of the Rydberg-Ritz expansion

    mu_l,j(n) = mu_0 + mu_0_prime / (n - mu_0)**2 .

Note that the tabulated defects are fine-structure resolved: there is no such
thing as a measured "spin-free" quantum defect. Collapsing the j dependence is
done in physics.quantum_defects and is a stated approximation of this model.
"""

QUANTUM_DEFECTS = {
    (0, 1/2): (3.1311804, 0.1784),

    (1, -1/2): (2.6548849, 0.2900),
    (1,  1/2): (2.6416737, 0.2950),

    (2, -1/2): (1.34809171, -0.60286),
    (2,  1/2): (1.34646572, -0.59600),

    (3, -1/2): (0.0165192, -0.085),
    (3,  1/2): (0.0165437, -0.086),
}

# Highest l for which the table above has entries. Above this the defect is
# taken to be zero, i.e. the high-l manifold is treated as hydrogenic.
L_MAX = 3
