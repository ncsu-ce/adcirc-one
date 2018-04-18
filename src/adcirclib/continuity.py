
g = 9.81
rho_0 = 1000


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
def tau0dzdt(lhs, rhs, nodes, z, dt, tau_0):

    num_nodes = nodes.shape[0]
    num_elements = num_nodes - 1
    z_prev = z[0]
    z_curr = z[1]

    for element in range(num_elements):

        n1 = element
        n2 = element + 1
        le = nodes[n2] - nodes[n1]
        tau_avg = (tau_0[n1] + tau_0[n2]) / 2
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
    alpha_1 = alpha[0]
    alpha_2 = alpha[1]
    alpha_3 = alpha[2]

    for element in range(num_elements):

        n1 = element
        n2 = element + 1
        le = nodes[n2] - nodes[n1]
        h_avg = (h[n1] + h[n2]) / 2

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
# Do these terms need to be divided by two??? (except the tau terms...)
#  - It would appear not, for version 1 (this version)
#  - Looks like version 2 will have the terms properly averaged
def jx_conservative(rhs, nodes, z, q, u, h, p, tau_0, tau_s, tau_b):

    num_nodes = nodes.shape[0]
    num_elements = num_nodes - 1

    z_curr = z[1]
    u_curr = u[1]
    q_curr = q[1]

    for element in range(num_elements):

        n1 = element
        n2 = element + 1
        le = nodes[n2] - nodes[n1]

        # Q * du/dx
        coeff = (q_curr[n2] * u_curr[n2] - q_curr[n1] * u_curr[n1]) / (2 * le)
        rhs[n1] += coeff
        rhs[n2] -= coeff

        # (g/2) * d(z^2)/dx
        coeff = g * (z_curr[n2]**2 - z_curr[n1]**2) / (4 * le)
        rhs[n1] += coeff
        rhs[n2] -= coeff

        # g * H * dXi/dx
        # Note: this is winds only. Tides to be added here.
        h1 = h[n1] + z_curr[n1]
        h2 = h[n2] + z_curr[n2]
        h_avg = (h1 + h2) / 2
        coeff = g * h_avg * (p[n2] - p[n1]) / (2 * le)
        rhs[n1] += coeff
        rhs[n2] -= coeff

        # tau_sx
        coeff = (tau_s[n1] + tau_s[n2]) / (2 * rho_0)
        rhs[n1] += coeff
        rhs[n2] -= coeff

        # tau_bx
        coeff = (tau_b[n1] + tau_b[n2]) / (2 * rho_0)
        rhs[n1] += coeff
        rhs[n2] -= coeff

        # tau_0 * q
        coeff = (tau_0[n1] * q_curr[n1] + tau_0[n2] * q_curr[n2]) / 2
        rhs[n1] += coeff
        rhs[n2] -= coeff


# (5)
# This is the ADCIRC version that pulls
# out Q as an elemental average
def qxdtaudx(rhs, nodes, q, tau_0):

    num_nodes = nodes.shape[0]
    num_elements = num_nodes - 1
    q_curr = q[1]

    for element in range(num_elements):

        n1 = element
        n2 = element + 1
        q_avg = (q_curr[n1] + q_curr[n2]) / 2

        coeff = q_avg / 2
        rhs[n1] -= coeff * tau_0[n1]
        rhs[n1] += coeff * tau_0[n2]
        rhs[n2] -= coeff * tau_0[n1]
        rhs[n2] += coeff * tau_0[n2]


# (6)
def jxboundary(rhs, dt, q, tau_0):

    l = 0
    r = rhs.shape[0] - 1

    # Left boundary
    q_prev = q[0][l]
    q_curr = q[1][l]
    q_next = q[2][l]
    dqdt = (q_next - q_prev) / (2 * dt)
    rhs[l] -= (dqdt + tau_0[l] * q_curr)

    # Right boundary
    q_prev = q[0][r]
    q_curr = q[1][r]
    q_next = q[2][r]
    dqdt = (q_next - q_prev) / (2 * dt)
    rhs[r] -= (dqdt + tau_0[r] * q_curr)
