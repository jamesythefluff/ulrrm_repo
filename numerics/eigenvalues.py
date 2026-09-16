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


def track_eigenstates(eigenvalues_list, eigenvectors_list):
    """
    Reorder a sequence of independently diagonalised (eigenvalues, eigenvectors)
    pairs so that column i follows the same adiabatic state along the whole
    sequence.

    numpy.linalg.eigh returns eigenvalues sorted by energy, so wherever two
    curves approach each other the columns swap and a curve plotted by column
    index jumps between physical states. Each step is therefore matched to the
    previous one by eigenvector overlap, and the pairing that maximises the
    total overlap is chosen with the Hungarian algorithm.

    Steps must be passed in the order they were computed, so that neighbouring
    entries correspond to neighbouring R.

    Arguments
    ----------
    eigenvalues_list : list
        list of 1D arrays, one per R step

    eigenvectors_list : list
        list of 2D arrays, one per R step, columns matching
        eigenvalues_list

    Returns
    -------
    tracked_eigenvalues : np.ndarray
        2D array, shape (n_steps, n_states)

    tracked_eigenvectors : list
        list of 2D arrays, reordered to match
    """
    n_steps = len(eigenvalues_list)
    n_states = len(eigenvalues_list[0])

    tracked_eigenvalues = np.zeros((n_steps, n_states))
    tracked_eigenvectors = [None] * n_steps

    # First step: nothing to compare against, so it sets the ordering.
    tracked_eigenvalues[0] = eigenvalues_list[0]
    tracked_eigenvectors[0] = eigenvectors_list[0]

    for step in range(1, n_steps):

        prev_vecs = tracked_eigenvectors[step - 1]
        curr_vecs = eigenvectors_list[step]
        curr_vals = eigenvalues_list[step]

        # overlap[i, j] = |<prev state i | current state j>|
        overlap = np.abs(prev_vecs.conj().T @ curr_vecs)

        # linear_sum_assignment minimises, so minimise -overlap to maximise
        # total overlap. perm[i] is the current column belonging in slot i.
        perm = linear_sum_assignment(-overlap)[1]

        tracked_eigenvalues[step] = curr_vals[perm]
        tracked_eigenvectors[step] = curr_vecs[:, perm]

    return tracked_eigenvalues, tracked_eigenvectors
