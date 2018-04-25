from math import sqrt
import numpy as np
import numpy.linalg as solver
import matplotlib.pyplot as plt
# from scipy.integrate import quadrature as integrate
import util_output as out
import util_stability as stability


# Mesh setup
num_nodes = 100
num_ts = 7200
dt = 10

# Problem setup
width = 100000  # 100km
depth = 100  # 100m
wind = 0
atm_pressure = 1013
mannings_n = 0.035
g = 9.80665
rho_0 = 1000
rho_air = 0.001293
flow_boundary_node = num_nodes - 1
flow_boundary_val = 0
elev_boundary_node = 0
elev_boundary_val = depth
wave_loc = int(num_nodes / 2)
out.center_column = wave_loc

# Solver setup
equations = 'nonconservative'
tau_0 = 1
alpha = [0.35, 0.3, 0.35]
instability_z_cutoff = 10


#
# Run ADCIRC
#

nodes = np.linspace(0, width, num=num_nodes, endpoint=True)
h = depth * np.ones(num_nodes)
z = np.zeros((3, num_nodes))
u = np.zeros((3, num_nodes))

wind = np.full(num_nodes, wind)
atm_pressure = np.full(num_nodes, atm_pressure)
mannings_n = np.full(num_nodes, mannings_n)
tau_s = np.zeros(num_nodes)
tau_b = np.zeros(num_nodes)

elev_boundary_val -= depth

# Make a tiny wave
z[0][wave_loc] += 0.001
z[1][wave_loc] += 0.001


def z_prev(n=None):

    if n is None:
        return z[0]
    return z[0][n]


def z_curr(n=None):

    if n is None:
        return z[1]
    return z[1][n]


def z_next(n=None):

    if n is None:
        return z[2]
    return z[2][n]


def q_prev(n=None):

    if n is None:
        return q[0]
    return q[0][n]


def q_curr(n=None):

    if n is None:
        return q[1]
    return q[1][n]


def u_prev(n=None):

    if n is None:
        return u[0]
    return u[0][n]


def u_curr(n=None):

    if n is None:
        return u[1]
    return u[1][n]


def phi_l(xl, xr):

    return lambda x: (xr - x) / (xr - xl)


def phi_r(xl, xr):

    return lambda x: (x - xl) / (xr - xl)


def dphi_l(xl, xr):

    return lambda x: -1 / (xr - xl)


def dphi_r(xl, xr):

    return lambda x: 1 / (xr - xl)


fig, ax = plt.subplots()
lhs = np.zeros((num_nodes, num_nodes))
rhs = np.zeros(num_nodes)
for ts in range(1, num_ts):

    time = ts * dt
    lhs.fill(0)
    rhs.fill(0)

    out.print_vector(h + z[1], time, 'H')

    # Check stability before advancing
    stability.check_unstable(z, 10)

    # Calculate Q
    q = (h + z) * u

    # Calculate tau_b
    cd = (g * mannings_n**2) / h**(1/3)
    tau_b = cd * np.fabs(u_curr()) * u_curr()

    # Calculate tau_s
    cd = 0.001 * (0.75 * 0.067 * np.fabs(wind))
    tau_s = rho_air * cd * np.fabs(wind) * wind

    # Assemble matrices and vectors for continuity solution
    for element in range(num_nodes - 1):

        n1 = element
        n2 = element + 1
        x1 = nodes[n1]
        x2 = nodes[n2]
        # N1 = phi_l(x1, x2)
        # N2 = phi_r(x1, x2)
        # dN1 = dphi_l(x1, x2)
        # dN2 = dphi_r(x1, x2)

        z1_curr = z_curr(n1)
        z1_prev = z_prev(n1)
        z2_curr = z_curr(n2)
        z2_prev = z_prev(n2)
        dz1 = z1_curr - z1_prev
        dz2 = z2_curr - z2_prev

        q1 = q_curr(n1)
        q2 = q_curr(n2)
        u1 = u_curr(n1)
        u2 = u_curr(n2)

        tau_s1 = tau_s[n1]
        tau_s2 = tau_s[n2]
        tau_b1 = tau_b[n1]
        tau_b2 = tau_b[n2]

        N11 = (x2 - x1) / 3  # integrate(lambda x: N1(x) * N1(x), x1, x2)[0]
        N12 = (x2 - x1) / 6  # integrate(lambda x: N1(x) * N2(x), x1, x2)[0]
        N21 = (x2 - x1) / 6  # integrate(lambda x: N2(x) * N1(x), x1, x2)[0]
        N22 = (x2 - x1) / 3  # integrate(lambda x: N2(x) * N2(x), x1, x2)[0]

        dN11 = 1 / (x2 - x1)  # integrate(lambda x: dN1(x) * dN1(x), x1, x2)[0]
        dN12 = -1 / (x2 - x1)  # integrate(lambda x: dN1(x) * dN2(x), x1, x2)[0]
        dN21 = -1 / (x2 - x1)  # integrate(lambda x: dN2(x) * dN1(x), x1, x2)[0]
        dN22 = 1 / (x2 - x1)  # integrate(lambda x: dN2(x) * dN2(x), x1, x2)[0]

        NdN11 = -0.5  # integrate(lambda x: N1(x) * dN1(x), x1, x2)[0]
        NdN12 = 0.5  # integrate(lambda x: N1(x) * dN2(x), x1, x2)[0]
        NdN21 = -0.5  # integrate(lambda x: N2(x) * dN1(x), x1, x2)[0]
        NdN22 = 0.5  # integrate(lambda x: N2(x) * dN2(x), x1, x2)[0]

        # (0)
        coeff = 1 / dt**2
        lhs[n1][n1] += coeff * N11
        lhs[n1][n2] += coeff * N21
        lhs[n2][n1] += coeff * N12
        lhs[n2][n2] += coeff * N22
        rhs[n1] += coeff * dz1 * N11
        rhs[n1] += coeff * dz2 * N21
        rhs[n2] += coeff * dz1 * N12
        rhs[n2] += coeff * dz2 * N22

        # (1)
        coeff = tau_0 / (2 * dt)
        lhs[n1][n1] += coeff * N11
        lhs[n1][n2] += coeff * N21
        lhs[n2][n1] += coeff * N12
        lhs[n2][n2] += coeff * N22
        rhs[n1] -= coeff * dz1 * dN11
        rhs[n1] -= coeff * dz2 * dN21
        rhs[n2] -= coeff * dz1 * dN12
        rhs[n2] -= coeff * dz2 * dN22

        # (2)
        h_avg = (h[n1] + h[n2]) / 2
        coeff = g * h_avg * alpha[0]
        lhs[n1][n1] += coeff * dN11
        lhs[n1][n2] += coeff * dN21
        lhs[n2][n1] += coeff * dN12
        lhs[n2][n2] += coeff * dN22
        coeff = g * h_avg * (alpha[1] + alpha[2])
        rhs[n1] -= coeff * dz1 * dN11
        rhs[n1] -= coeff * dz2 * dN21
        rhs[n2] -= coeff * dz1 * dN12
        rhs[n2] -= coeff * dz2 * dN22
        coeff = g * h_avg
        rhs[n1] += coeff * z1_prev * dN11
        rhs[n1] += coeff * z2_prev * dN21
        rhs[n2] += coeff * z1_prev * dN12
        rhs[n2] += coeff * z2_prev * dN22

        # (3)
        if equations == 'nonconservative':

            # (a)
            q_avg = (q1 + q2) / 2
            rhs[n1] -= q_avg * u1 * dN11
            rhs[n1] -= q_avg * u2 * dN21
            rhs[n2] -= q_avg * u1 * dN12
            rhs[n2] -= q_avg * u2 * dN22

            # (b)
            coeff = g / 2
            z1sq = z1_curr**2
            z2sq = z2_curr**2
            rhs[n1] -= coeff * z1sq * dN11
            rhs[n1] -= coeff * z2sq * dN21
            rhs[n2] -= coeff * z1sq * dN12
            rhs[n2] -= coeff * z2sq * dN22

            # (c)
            coeff = 1 / rho_0
            rhs[n1] += coeff * tau_s1 * NdN11
            rhs[n1] += coeff * tau_s2 * NdN21
            rhs[n2] += coeff * tau_s1 * NdN12
            rhs[n2] += coeff * tau_s2 * NdN22

            # (d)
            rhs[n1] -= tau_b1 * NdN11
            rhs[n1] -= tau_b2 * NdN21
            rhs[n2] -= tau_b1 * NdN12
            rhs[n2] -= tau_b2 * NdN22

            # (e)
            u_avg = (u1 + u2) / 2
            dzdt1 = dz1 / dt
            dzdt2 = dz2 / dt
            rhs[n1] += u_avg * dzdt1 * NdN11
            rhs[n1] += u_avg * dzdt2 * NdN21
            rhs[n2] += u_avg * dzdt1 * NdN12
            rhs[n2] += u_avg * dzdt2 * NdN22

            # (f)
            rhs[n1] += tau_0 * q1 * NdN11
            rhs[n1] += tau_0 * q2 * NdN21
            rhs[n2] += tau_0 * q1 * NdN12
            rhs[n2] += tau_0 * q2 * NdN22

        # (4)
        # Because we are currently assuming tau_0 is constant,
        # dtau_0/dx = 0, so this term does not contribute.

        # (5)
        if n1 == flow_boundary_node or n2 == flow_boundary_node:
            qb_next = flow_boundary_val
            qb_curr = flow_boundary_val
            qb_prev = flow_boundary_val
            dqdt = (qb_next - qb_prev) / (2 * dt)  # Note: This will currently always be zero
            rhs[flow_boundary_node] -= dqdt
            rhs[flow_boundary_node] -= tau_0 * qb_curr

    # Apply elevation specified boundary condition to continuity equation
    rms = 0
    lhs[elev_boundary_node].fill(0)

    for row in range(num_nodes):

        rms += lhs[row][row]**2
        rhs[row] -= lhs[row][elev_boundary_node] * elev_boundary_val
        lhs[row][elev_boundary_node] = 0

    rms = sqrt(rms / (num_nodes - 1))
    lhs[elev_boundary_node][elev_boundary_node] = rms
    rhs[elev_boundary_node] = rms * elev_boundary_val

    # Solve GWCE
    # out.print_matrix(lhs, time, 'LHS')
    # out.print_vector(rhs, time, 'rhs')
    dz = solver.solve(lhs, rhs)
    z[2] = z[1] + dz

    # Set up matrices and vectors for momentum
    lhs.fill(0)
    rhs.fill(0)

    for element in range(num_nodes - 1):

        n1 = element
        n2 = element + 1
        x1 = nodes[n1]
        x2 = nodes[n2]
        # N1 = phi_l(x1, x2)
        # N2 = phi_r(x1, x2)
        # dN1 = dphi_l(x1, x2)
        # dN2 = dphi_r(x1, x2)

        u1 = u_curr(n1)
        u2 = u_curr(n2)
        h1 = h[n1]
        h2 = h[n2]
        z1_curr = z_curr(n1)
        z1_next = z_next(n1)
        z2_curr = z_curr(n2)
        z2_next = z_next(n2)

        tau_s1 = tau_s[n1]
        tau_s2 = tau_s[n2]

        N11 = (x2 - x1) / 3  # integrate(lambda x: N1(x) * N1(x), x1, x2)[0]
        N12 = (x2 - x1) / 6  # integrate(lambda x: N1(x) * N2(x), x1, x2)[0]
        N21 = (x2 - x1) / 6  # integrate(lambda x: N2(x) * N1(x), x1, x2)[0]
        N22 = (x2 - x1) / 3  # integrate(lambda x: N2(x) * N2(x), x1, x2)[0]

        dNN11 = -0.5  # integrate(lambda x: dN1(x) * N1(x), x1, x2)[0]
        dNN12 = -0.5  # integrate(lambda x: dN1(x) * N2(x), x1, x2)[0]
        dNN21 = 0.5  # integrate(lambda x: dN2(x) * N1(x), x1, x2)[0]
        dNN22 = 0.5  # integrate(lambda x: dN2(x) * N2(x), x1, x2)[0]

        # (0)
        coeff = 1 / dt
        lhs[n1][n1] += coeff * N11
        lhs[n1][n2] += coeff * N21
        lhs[n2][n1] += coeff * N12
        lhs[n2][n2] += coeff * N22
        rhs[n1] += coeff * u1 * N11
        rhs[n1] += coeff * u2 * N21
        rhs[n2] += coeff * u1 * N12
        rhs[n2] += coeff * u2 * N22

        # (1)
        u_avg = (u1 + u2) / 2
        rhs[n1] -= u_avg * u1 * dNN11
        rhs[n1] -= u_avg * u2 * dNN21
        rhs[n2] -= u_avg * u1 * dNN12
        rhs[n2] -= u_avg * u2 * dNN22

        # (2)
        coeff = g / 2
        rhs[n1] -= coeff * (z1_curr + z1_next) * dNN11
        rhs[n1] -= coeff * (z2_curr + z2_next) * dNN21
        rhs[n2] -= coeff * (z1_curr + z1_next) * dNN12
        rhs[n2] -= coeff * (z2_curr + z2_next) * dNN22

        # (3) IFNLFA = 0, so H1, H2 = h1, h2, i.e. not dependent on dz
        #     but keep these separate for when we change this
        coeff = 1 / 2
        rhs[n1] += coeff * (tau_s1 / h1) * N11
        rhs[n1] += coeff * (tau_s2 / h2) * N21
        rhs[n2] += coeff * (tau_s1 / h1) * N12
        rhs[n2] += coeff * (tau_s2 / h2) * N22

        rhs[n1] += coeff * (tau_s1 / h1) * N11
        rhs[n1] += coeff * (tau_s2 / h2) * N21
        rhs[n2] += coeff * (tau_s1 / h1) * N12
        rhs[n2] += coeff * (tau_s2 / h2) * N22

        # (4)
        coeff = 1 / 2
        cd_1 = cd[n1]
        cd_2 = cd[n2]
        k_slip_1 = cd_1 * abs(u1)
        k_slip_2 = cd_2 * abs(u2)
        lhs[n1][n1] += coeff * (k_slip_1 / h1) * N11
        lhs[n1][n2] += coeff * (k_slip_2 / h2) * N21
        lhs[n2][n1] += coeff * (k_slip_1 / h1) * N12
        lhs[n2][n2] += coeff * (k_slip_2 / h2) * N22
        rhs[n1] -= coeff * (k_slip_1 * u1 / h1) * N11
        rhs[n1] -= coeff * (k_slip_2 * u2 / h2) * N21
        rhs[n2] -= coeff * (k_slip_1 * u1 / h1) * N12
        rhs[n2] -= coeff * (k_slip_2 * u2 / h2) * N22

    # Solve momentum
    if equations == 'nonconservative':
        u[2] = solver.solve(lhs, rhs)

    # Advance timestep
    z = np.roll(z, -1, axis=0)
    u = np.roll(u, -1, axis=0)

    z[2].fill(0)
    u[2].fill(0)

    plt.cla()
    ax.set_ylim(depth-0.1, depth+0.1)
    ax.plot(h+z[1])
    plt.pause(0.05)
