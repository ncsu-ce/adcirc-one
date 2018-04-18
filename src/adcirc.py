import numpy as np
import numpy.linalg as solver
from adcirclib.continuity import *


# User defined variables
num_nodes = 5
num_ts = 100
dt = 0.5
width = 100
tau_0 = 1
alpha = [0, 1, 0]


# Run ADCIRC
nodes = np.linspace(0, width, num=num_nodes, endpoint=True)

h = np.ones((3, num_nodes))
z = np.zeros((3, num_nodes))
u = np.zeros((3, num_nodes))
tau_0 = np.full(num_nodes, tau_0)


for ts in range(num_ts):

    q = (h + z) * u

    time = ts * dt
    lhs = np.zeros((num_nodes, num_nodes))
    rhs = np.zeros(num_nodes)

    # Solve continuity to get dz
    d2zdt2(lhs, rhs, nodes, z, dt)
    tau0dzdt(lhs, rhs, nodes, z, dt, tau_0)
    ghdzdx(lhs, rhs, nodes, z, alpha, h)
    jx_conservative(lhs, rhs)
    qxdtaudx(rhs, nodes, q, tau_0)
    jxboundary(lhs, rhs)

    dz = solver.solve(lhs, rhs)
    print(lhs, '*', dz, '=', rhs, sep='\t')

    # Solve momentum to get u


print('Done.')
