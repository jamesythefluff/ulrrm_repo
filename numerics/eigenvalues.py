import numpy as np
from scipy.optimize import linear_sum_assignment

def diagonalise(H):
     return np.linalg.eigh(H)

def label_eigenstates(large_R_eigenvectors, basis):

    labels = []
     
    for i in range(len(basis)):

        eigenvector = large_R_eigenvectors[:, i]

        overlaps = np.abs(eigenvector)**2

        basis_index = np.argmax(overlaps)

        n, l = basis[basis_index]

        labels.append((n, l))

    return labels

def adiabatic_reordering(prev_eigvecs, eigenvectors, eigenvalues):

    # Reordering to track adiabatic states
    if prev_eigvecs is not None:
        overlap = np.abs(prev_eigvecs.conj().T @ eigenvectors)
        perm = linear_sum_assignment(-overlap)[1]   # scipy.optimize
        eigenvalues = eigenvalues[perm]
        eigenvectors = eigenvectors[:, perm]

    prev_eigvecs = eigenvectors

    return eigenvalues, eigenvectors, prev_eigvecs
