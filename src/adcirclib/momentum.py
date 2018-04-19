import numpy as np


g = 9.81
rho_0 = 1000


# (1)
def dqdt(lhs, rhs, nodes, q):

    num_nodes = nodes.shape[0]

    q_curr = q[1]

    for node in range(num_nodes):

        lhs[node][node] += 1
        rhs[node] += q_curr[node]


# (2)
def dqudx(rhs, nodes, dt, q, u):

    num_nodes = nodes.shape[0]
    num_elements = num_nodes - 1
    q_curr = q[1]
    u_curr = u[1]

    rhs_temp = np.zeros_like(rhs)

    for element in range(num_elements):

        n1 = element
        n2 = element + 1

        q1 = q_curr[n1]
        q2 = q_curr[n2]
        u1 = u_curr[n1]
        u2 = u_curr[n2]

        coeff = q2*u2 - q1*u1
        rhs_temp[n1] += coeff
        rhs_temp[n2] += coeff

    for node in range(num_nodes):

        le_l = 0
        le_r = 0
        if node > 0:
            le_l = nodes[node] - nodes[node - 1]
        if node < num_nodes - 1:
            le_r = nodes[node + 1] - nodes[node]
        le_n = le_l + le_r

        rhs[node] -= (dt * le_n) * rhs_temp[node]


# (3)
def dzdx(rhs, nodes, dt, h, z, p):

    num_nodes = nodes.shape[0]
    num_elements = num_nodes - 1
    z_curr = z[1]
    z_next = z[2]

    rhs_temp = np.zeros_like(rhs)

    for element in range(num_elements):

        n1 = element
        n2 = element + 1

        # Current timestep contribution
        z1 = z_curr[n1]
        z2 = z_curr[n2]
        h1 = h[n1] + z1
        h2 = h[n2] + z2
        h_avg = (h1 + h2) / 2
        dz = z2 - z1
        dp = p[n2] - p[n1]

        coeff = h_avg * (dz + dp) / 2
        rhs_temp[n1] += coeff
        rhs_temp[n2] += coeff

        # Next timestep contribution
        z1 = z_next[n1]
        z2 = z_next[n2]
        h1 = h[n1] + z1
        h2 = h[n2] + z2
        h_avg = (h1 + h2) / 2
        dz = z2 - z1
        dp = p[n2] - p[n1]

        coeff = h_avg * (dz + dp) / 2
        rhs_temp[n1] += coeff
        rhs_temp[n2] += coeff

    for node in range(num_nodes):

        le_l = 0
        le_r = 0
        if node > 0:
            le_l = nodes[node] - nodes[node - 1]
        if node < num_nodes - 1:
            le_r = nodes[node + 1] - nodes[node]
        le_n = le_l + le_r

        rhs[node] -= (dt / le_n) * rhs_temp[node]


# (4)
def tausx(rhs, nodes, dt, tau_s):

    num_nodes = nodes.shape[0]
    tau_s_curr = tau_s
    tau_s_next = tau_s

    for node in range(num_nodes):

        rhs[node] += (dt / (2 * rho_0)) * (tau_s_curr[node] + tau_s_next[node])


# (5)
def taubx(lhs, rhs, nodes, dt, q, u, h, z, n):

    num_nodes = nodes.shape[0]
    q_curr = q[1]
    u_curr = u[1]
    z_curr = z[1]
    z_next = z[2]
    h_curr = h + z_curr
    h_next = h + z_next
    cd_curr = (g * n**2) / (h + z_curr)**(1.0/3.0)
    k_slip = cd_curr * u_curr

    for node in range(num_nodes):

        lhs[node][node] += (dt * k_slip[node]) / (2 * h_next[node])
        rhs[node] -= (dt * k_slip[node] * q_curr[node]) / (2 * h_curr[node])
