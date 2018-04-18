import numpy as np
import numpy.linalg as solver
from adcirclib.continuity import *
from adcirclib.stress_bottom import calculate_bottom_stress_mannings_n
from adcirclib.stress_surface import calculate_wind_stress

#
# User defined variables
#
num_nodes = 5
num_ts = 100
dt = 0.5
width = 100
wind = 1
atm_pressure = 1013
manning = 0.035
tau_0 = 1
alpha = [0, 1, 0]


#
# Initialization
#
nodes = np.linspace(0, width, num=num_nodes, endpoint=True)

h = np.ones(num_nodes)
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
    tau_b = calculate_bottom_stress_mannings_n(u[1], h, z[1], manning)
    tau_s = calculate_wind_stress(wind)

    # Solve continuity to get dz
    d2zdt2(lhs, rhs, nodes, z, dt)
    tau0dzdt(lhs, rhs, nodes, z, dt, tau_0)
    ghdzdx(lhs, rhs, nodes, z, alpha, h)
    jx_conservative(rhs, nodes, z, q, u, h, atm_pressure, tau_0, tau_s, tau_b)
    qxdtaudx(rhs, nodes, q, tau_0)
    jxboundary(rhs, dt, q, tau_0)

    dz = solver.solve(lhs, rhs)
    print(lhs, '*', dz, '=', rhs, sep='\t')

    # Solve momentum to get u


print('Done.')
