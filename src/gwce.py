import numpy as np

g = 9.81


def initialize_gwce(num_nodes):

    lhs = np.zeros((num_nodes, num_nodes), dtype=float)


def build_gwce_lhs(lhs, dt, nodes, alpha, h):

    num_nodes = nodes.shape[0]
    num_elements = num_nodes - 1
    alpha_1 = alpha[0]

    for element in range(num_elements):

        n1 = element
        n2 = element + 1

        le = nodes[n2] - nodes[n1]      # element length

        # dddz/dt contribution
        dddz = le / (6 * dt**2)
        lhs[n1][n1] += 2 * dddz
        lhs[n1][n2] += dddz
        lhs[n2][n1] += dddz
        lhs[n2][n2] += 2 * dddz

        # g*h*dz/dx contribution
        h_avg = (h[n1] + h[n2]) / 2
        dzdx = g * alpha_1 * h_avg / le
        lhs[n1][n1] += dzdx
        lhs[n1][n2] -= dzdx
        lhs[n2][n1] -= dzdx
        lhs[n2][n2] += dzdx

    return lhs
