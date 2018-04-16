import numpy as np


g = 9.81                    # Gravity
cf = 1                      # Bottom friction coefficient used unless otherwise specified TODO: Set this value
friction_type = 'linear'    # Can be 'linear', 'quadratic', or 'hybrid' TODO: Implement quadratic and hybrid


def calculate_bottom_friction_mannings_n(u, h, z, n):

    # Determine number of nodes
    num_nodes = u.shape[1]

    # Create iterator that will calculate Cd for each node
    iterator = ((g * n[i]**2) / ((h[i] + z[i])**(1.0/3.0)) for i in range(num_nodes))

    # Create new array using iterator
    cd = np.fromiter(iterator, float, count=num_nodes)

    # Apply lower limit
    np.clip(cd, cf, None, out=cd)

    # Calculate bottom friction based on Cd
    return calculate_bottom_friction_cd(u, h, z, cd)


def calculate_bottom_friction_chezy(u, h, z, c):

    # Determine number of nodes
    num_nodes = u.shape[1]

    # Create iterator that will calculate Cd for each node
    iterator = ((g / c[i]**2) for i in range(num_nodes))

    # Create new array using iterator
    cd = np.fromiter(iterator, float, count=num_nodes)

    # Calculate bottom friction based on Cd
    return calculate_bottom_friction_cd(u, h, z, cd)


def calculate_bottom_friction_cd(u, h, z, cd):

    if friction_type == 'linear':
        return cd

    if friction_type == 'quadratic':
        print('Quadtratic bottom friction not yet implemented. Falling back on linear.')
        return cd

    if friction_type == 'hybrid':
        print('Hybrid bottom friction not yet implemented. Falling back on linear.')
        return cd
