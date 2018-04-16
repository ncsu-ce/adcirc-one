import numpy as np

# User defined variables
num_elements = 3
width = 3
tau_0 = 1
g = 9.81
dt = 0.5

# Constants
num_nodes = num_elements + 1
element_width = width / num_elements

# Solve

# x-coordinates of nodes
x = np.linspace(0, width, num=num_nodes, endpoint=True)

# Matrices
#
# Each will be (n x 3) for past, present, and future timesteps
#
# Initialize depth to one meter
h = np.ones((num_nodes, 3))     # bathymetric depth
z = np.zeros((num_nodes, 3))    # free-surface departure from geoid
H = h + z                       # total water-column thickness
U = np.zeros((num_nodes, 3))    # x-velocity
Q = H * U                       # x-flux

print('h =')
print(h)
print('z =')
print(z)
print('H =')
print(H)
print('U =')
print(U)
print('Q =')
print(Q)


# zeta = np.zeros(num_nodes)      # free-surface departure from geoid
# H = h + zeta                    # total water column thickness

# K = np.zeros((num_nodes, num_nodes))
# M = np.zeros((num_nodes, num_nodes))
# F = np.zeros((num_nodes, num_nodes))


