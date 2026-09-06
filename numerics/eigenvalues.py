import numpy as np

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
