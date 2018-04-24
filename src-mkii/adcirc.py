import numpy as np
from scipy.integrate import quadrature as integrate


# Problem setup
width = 1000
depth = 5
wind = 1
atm_pressure = 1013
mannings_n = 0.035
g = 9.81

# Mesh setup
num_nodes = 100
num_ts = 7200
dt = 0.5

# Solver setup
tau_0 = 1
alpha = [0.35, 0.3, 0.35]

# Run ADCIRC
nodes = np.linspace(0, width, num=num_nodes, endpoint=True)
h = depth * np.ones(num_nodes)
z = np.zeros((3, num_nodes))
u = np.zeros((3, num_nodes))

wind = np.full(num_nodes, wind)
atm_pressure = np.full(num_nodes, atm_pressure)
mannings_n = np.full(num_nodes, mannings_n)


def z_prev(n):

    return z[0][n]


def z_curr(n):

    return z[1][n]


def phi_l(xl, xr):

    return lambda x: (xr - x) / (xr - xl)


def phi_r(xl, xr):

    return lambda x: (x - xl) / (xr - xl)


def dphi_l(xl, xr):

    return lambda x: -1 / (xr - xl)


def dphi_r(xl, xr):

    return lambda x: 1 / (xr - xl)


lhs = np.zeros((num_nodes, num_nodes))
rhs = np.zeros(num_nodes)
for ts in range(2, num_ts):

    time = ts * dt
    lhs.fill(0)
    rhs.fill(0)

    q = (h + z) * u
    # Calculate tau_b
    # Calculate tau_s

    for element in range(num_nodes - 1):

        n1 = element
        n2 = element + 1
        x1 = nodes[n1]
        x2 = nodes[n2]
        N1 = phi_l(x1, x2)
        N2 = phi_r(x1, x2)
        dN1 = dphi_l(x1, x2)
        dN2 = dphi_r(x1, x2)

        z1_curr = z_curr(n1)
        z1_prev = z_prev(n1)
        z2_curr = z_curr(n2)
        z2_prev = z_prev(n2)
        dz1 = z1_curr - z1_prev
        dz2 = z2_curr - z2_prev

        N11 = integrate(lambda x: N1(x) * N1(x), x1, x2)[0]
        N12 = integrate(lambda x: N1(x) * N2(x), x1, x2)[0]
        N21 = integrate(lambda x: N2(x) * N1(x), x1, x2)[0]
        N22 = integrate(lambda x: N2(x) * N2(x), x1, x2)[0]

        dN11 = integrate(lambda x: dN1(x) * dN1(x), x1, x2)[0]
        dN12 = integrate(lambda x: dN1(x) * dN2(x), x1, x2)[0]
        dN21 = integrate(lambda x: dN2(x) * dN1(x), x1, x2)[0]
        dN22 = integrate(lambda x: dN2(x) * dN2(x), x1, x2)[0]

        # (0)
        coeff = 1 / dt**2
        lhs[n1][n1] += coeff * N11
        lhs[n1][n2] += coeff * N12
        lhs[n2][n1] += coeff * N21
        lhs[n2][n2] += coeff * N22
        rhs[n1][n1] += coeff * dz1 * N11
        rhs[n1][n2] += coeff * dz1 * N12
        rhs[n2][n1] += coeff * dz2 * N21
        rhs[n2][n2] += coeff * dz2 * N22

        # (1)
        coeff = tau_0 / (2 * dt)
        lhs[n1][n1] += coeff * N11
        lhs[n1][n2] += coeff * N12
        lhs[n2][n1] += coeff * N21
        lhs[n2][n2] += coeff * N22

        # (2)
        h_avg = (h[n1] + h[n2]) / 2
        coeff = g * h_avg * alpha[0]
        lhs[n1][n1] += coeff * dN11
        lhs[n1][n2] += coeff * dN12
        lhs[n2][n1] += coeff * dN21
        lhs[n2][n2] += coeff * dN22
