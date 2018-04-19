import numpy as np
import scipy as sp
import numpy.linalg as solver
from adcirclib.continuity import *
from adcirclib.momentum import *
from adcirclib.boundaries import *
from adcirclib.stress_bottom import calculate_bottom_stress_mannings_n
from adcirclib.stress_surface import calculate_wind_stress

#
# User defined variables
#
num_nodes = 1000
num_ts = 7200
dt = 0.5
width = 1000
depth = 1
wind = 1
atm_pressure = 1013
manning = 0.035
tau_0 = 1
alpha = [0.35, 0.3, 0.35]
# alpha = [0, 1, 0]
np.set_printoptions(precision=3, linewidth=175)

essential_boundaries = {
    0: 0
}

#
# Initialization
#
nodes = np.linspace(0, width, num=num_nodes, endpoint=True)

h = depth * np.ones(num_nodes)
z = np.zeros((3, num_nodes))
u = np.zeros((3, num_nodes))

# Wind, atmospheric pressure, bottom stress, and tau all constant over whole domain
wind = np.full(num_nodes, wind)
atm_pressure = np.full(num_nodes, atm_pressure)
manning = np.full(num_nodes, manning)
tau_0 = np.full(num_nodes, tau_0)


#
# Timestep loop
#
for ts in range(num_ts):

    time = ts * dt
    lhs = np.zeros((num_nodes, num_nodes))
    rhs = np.zeros(num_nodes)

    q = (h + z) * u
    tau_b = calculate_bottom_stress_mannings_n(q[1], h, z[1], manning)
    tau_s = calculate_wind_stress(wind)

    # Solve continuity to get dz
    d2zdt2(lhs, rhs, nodes, z, dt)
    tau0dzdt(lhs, rhs, nodes, z, dt, tau_0)
    ghdzdx(lhs, rhs, nodes, z, alpha, h)
    jx_conservative(rhs, nodes, z, q, u, h, atm_pressure, tau_0, tau_s, tau_b)
    qxdtaudx(rhs, nodes, q, tau_0)
    jxboundary(rhs, dt, q, tau_0)

    # Doesn't work
    # elevation_specified(essential_boundaries, lhs, rhs)
    # dz = solver.solve(lhs, rhs)

    # At least does something
    lhs, rhs = elevation_boundaries(essential_boundaries, lhs, rhs)
    dz = solver.solve(lhs, rhs)
    for node, value in essential_boundaries.items():
        dz = np.insert(dz, node, value)

    z[2] = z[1] + dz

    # Solve momentum to get u
    lhs = np.zeros((num_nodes, num_nodes))
    rhs = np.zeros(num_nodes)

    dqdt(lhs, rhs, nodes, q)
    dqudx(rhs, nodes, dt, q, u)
    dzdx(rhs, nodes, dt, h, z, atm_pressure)
    tausx(rhs, nodes, dt, tau_s)
    taubx(lhs, rhs, nodes, dt, q, u, h, z, manning)

    q[2] = solver.solve(lhs, rhs)
    u[2] = q[2] / h

    # Advance the timestep
    z = np.roll(z, -1, axis=0)
    u = np.roll(u, -1, axis=0)

    print((h + z[1])[-25:])
    if np.any(np.isnan(z[1])):
        break

print('Done.')
