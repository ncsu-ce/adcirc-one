from math import sqrt
import numpy as np


def elevation_boundaries(boundaries, lhs, rhs):

    num_nodes = lhs.shape[0]
    dof = num_nodes - len(boundaries)

    lhs_mask = np.ones(lhs.shape, dtype=bool)
    rhs_mask = np.ones(rhs.shape, dtype=bool)

    for node, value in boundaries.items():

        lhs_mask[node, :] = False
        lhs_mask[:, node] = False
        rhs_mask[node] = False

        for row in range(num_nodes):

            lhs_val = lhs[row, node]
            rhs[row] -= lhs_val * value

    lhs = np.reshape(lhs[lhs_mask], (dof, dof))
    rhs = np.reshape(rhs[rhs_mask], (dof, 1))

    return lhs, rhs


def elevation_specified(boundaries, lhs, rhs):

    num_nodes = lhs.shape[0]
    rms = 0

    for node in boundaries:

        lhs[node].fill(0)

    for row in range(num_nodes):

        rms += lhs[row][row]

    rms = sqrt(rms / (num_nodes - len(boundaries)))

    for node, value in boundaries.items():

        rhs[node] = rms * value

        for row in range(num_nodes):

            if row not in boundaries:

                rhs[row] -= value * lhs[row][node]
