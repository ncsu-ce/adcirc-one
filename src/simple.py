import numpy as np

# User defined variables
num_elements = 3
width = 3
tau_0 = 1
g = 9.81

# Constants
num_nodes = num_elements + 1
element_width = width / num_elements

# Solve
x = np.linspace(0, width, num=num_nodes, endpoint=True)
h = np.ones(num_nodes)
H = np.ones(num_nodes)
K = np.zeros((num_nodes, num_nodes))
M = np.zeros((num_nodes, num_nodes))
F = np.zeros((num_nodes, num_nodes))

for node in range(num_nodes):

    # Term (1) contribution
    if node - 1 >= 0:
        K[node][node-1] += 0.5 * tau_0
    if node + 1 < num_nodes:
        K[node][node+1] -= 0.5 * tau_0

    # Term (2) contribution
    if node - 1 >= 0:

        h_bar_l = (h[node-1] + h[node]) / 2.0
        K[node][node-1] -= g * h_bar_l
        K[node][node] += g * h_bar_l

    if node + 1 < num_nodes:

        h_bar_r = (h[node] + h[node+1]) / 2.0
        K[node][node] -= g * h_bar_r
        K[node][node+1] += g * h_bar_r


print(K)
