
g = 9.81


# (1)
def d2zdt2(lhs, rhs, nodes, z, dt):

    num_nodes = nodes.shape[0]
    num_elements = num_nodes - 1
    z_prev = z[0]
    z_curr = z[1]

    for element in range(num_elements):

        n1 = element
        n2 = element + 1
        le = nodes[n2] - nodes[n1]

        dz1 = z_curr[n1] - z_prev[n1]
        dz2 = z_curr[n2] - z_prev[n2]

        coeff = le / (6 * dt**2)
        lhs[n1][n1] += 2 * coeff
        lhs[n1][n2] += coeff
        lhs[n2][n1] += coeff
        lhs[n2][n2] += 2 * coeff

        rhs[n1] += 2 * coeff * dz1
        rhs[n1] += coeff * dz2
        rhs[n2] += coeff * dz1
        rhs[n2] += 2 * coeff * dz2


# (2)
def tau0dzdt(lhs, rhs, nodes, z, dt, tau):

    num_nodes = nodes.shape[0]
    num_elements = num_nodes - 1
    z_prev = z[0]
    z_curr = z[1]

    for element in range(num_elements):

        n1 = element
        n2 = element + 1
        le = nodes[n2] - nodes[n1]
        tau_avg = (tau[n1] + tau[n2]) / 2
        dz1 = z_curr[n1] - z_prev[n1]
        dz2 = z_curr[n2] - z_prev[n2]

        coeff = (tau_avg + le) / (12 * dt)
        lhs[n1][n1] += 2 * coeff
        lhs[n1][n2] += coeff
        lhs[n2][n1] += coeff
        lhs[n2][n2] += 2 * coeff

        rhs[n1] -= 2 * coeff * dz1
        rhs[n1] -= coeff * dz2
        rhs[n2] -= coeff * dz1
        rhs[n2] -= 2 * coeff * dz2


# (3)
def ghdzdx(lhs, rhs, nodes, z, alpha, h):

    num_nodes = nodes.shape[0]
    num_elements = num_nodes - 1
    z_prev = z[0]
    z_curr = z[1]
    h_curr = h[1]
    alpha_1 = alpha[0]
    alpha_2 = alpha[1]
    alpha_3 = alpha[2]

    for element in range(num_elements):

        n1 = element
        n2 = element + 1
        le = nodes[n2] - nodes[n1]
        h_avg = (h_curr[n1] + h_curr[n2]) / 2

        coeff = g * alpha_1 * h_avg / le
        lhs[n1][n1] += coeff
        lhs[n1][n2] -= coeff
        lhs[n2][n1] -= coeff
        lhs[n2][n2] += coeff

        rhs[n1] -= coeff * ((alpha_1 + alpha_2) * z_curr[n1] + alpha_3 * z_prev[n1])
        rhs[n1] += coeff * ((alpha_1 + alpha_2) * z_curr[n2] + alpha_3 * z_prev[n2])
        rhs[n2] += coeff * ((alpha_1 + alpha_2) * z_curr[n1] + alpha_3 * z_prev[n1])
        rhs[n2] -= coeff * ((alpha_1 + alpha_2) * z_curr[n2] + alpha_3 * z_prev[n2])


# (4)
def jx(lhs, rhs):

    pass


# (5)
# This is the ADCIRC version that pulls
# out Q as an elemental average
def qxdtaudx(rhs, nodes, q, tau):

    num_nodes = nodes.shape[0]
    num_elements = num_nodes - 1
    q_curr = q[1]

    for element in range(num_elements):

        n1 = element
        n2 = element + 1
        q_avg = (q_curr[n1] + q_curr[n2]) / 2

        coeff = q_avg / 2
        rhs[n1] -= coeff * tau[n1]
        rhs[n1] += coeff * tau[n2]
        rhs[n2] -= coeff * tau[n1]
        rhs[n2] += coeff * tau[n2]


# (6)
def jxboundary(lhs, rhs):

    pass
