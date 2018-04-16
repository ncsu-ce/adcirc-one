import numpy as np

g = 9.81


def initialize_gwce(num_nodes):

    lhs = np.zeros((num_nodes, num_nodes), dtype=float)


def build_gwce_rhs(rhs, dt, nodes, alpha, z_current, z_previous, tau, q, h):

    num_nodes = nodes.shape[0]
    num_elements = num_nodes - 1
    alpha_1 = alpha[0]
    alpha_2 = alpha[1]
    alpha_3 = alpha[2]

    for element in range(num_elements):

        n1 = element
        n2 = element + 1

        le = nodes[n2] - nodes[n1]
        dz1 = z_current[n1] - z_previous[n1]
        dz2 = z_current[n2] - z_previous[n2]

        # (5) contribution
        tau_avg = (tau[n1] + tau[n2]) / 2
        coeff = tau_avg / 3
        rhs[n1] += 2 * coeff * q[n1]
        rhs[n1] += coeff * q[n2]
        rhs[n2] += coeff * q[n1]
        rhs[n2] += 2 * coeff * q[n2]

        # (6) contribution
        coeff = le / (6 * dt)
        rhs[n1] += 2 * coeff * dz1
        rhs[n1] += coeff * dz2
        rhs[n2] += coeff * dz1
        rhs[n2] += 2 * coeff * dz2

        # (7) contribution
        coeff = 1 / (12 * dt)
        rhs[n1] -= 2 * coeff * dz1
        rhs[n1] -= coeff * dz2
        rhs[n2] -= coeff * dz1
        rhs[n2] -= 2 * coeff * dz2

        # (8) contribution
        h_avg = (h[n1] + h[n2]) / 2
        coeff = g * h_avg / le
        rhs[n1] -= coeff * ((alpha_1 + alpha_2) * z_current[n1] + alpha_3 * z_previous[n1])
        rhs[n1] += coeff * ((alpha_1 + alpha_2) * z_current[n2] + alpha_3 * z_previous[n2])
        rhs[n2] -= coeff * ((alpha_1 + alpha_2) * z_current[n1] + alpha_3 * z_previous[n1])
        rhs[n2] += coeff * ((alpha_1 + alpha_2) * z_current[n2] + alpha_3 * z_previous[n2])


def build_gwce_lhs(lhs, dt, nodes, alpha, h, tau):

    num_nodes = nodes.shape[0]
    num_elements = num_nodes - 1
    alpha_1 = alpha[0]

    for element in range(num_elements):

        n1 = element
        n2 = element + 1

        le = nodes[n2] - nodes[n1]      # element length

        # (1) contribution
        coeff = le / (6 * dt**2)
        lhs[n1][n1] += 2 * coeff
        lhs[n1][n2] += coeff
        lhs[n2][n1] += coeff
        lhs[n2][n2] += 2 * coeff

        # (2) contribution
        tau_avg = (tau[n1] + tau[n2]) / 2
        coeff = (tau_avg + le) / (12 * dt)
        lhs[n1][n1] += 2 * coeff
        lhs[n1][n2] += coeff
        lhs[n2][n1] += coeff
        lhs[n2][n2] += 2 * coeff

        # g*h*dz/dx contribution
        h_avg = (h[n1] + h[n2]) / 2
        coeff = g * alpha_1 * h_avg / le
        lhs[n1][n1] += coeff
        lhs[n1][n2] -= coeff
        lhs[n2][n1] -= coeff
        lhs[n2][n2] += coeff

    return lhs
