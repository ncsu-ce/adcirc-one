import numpy as np


g = 9.81                    # Gravity
cf = 1                      # Bottom friction coefficient used unless otherwise specified TODO: Set this value
friction_type = 'linear'    # Can be 'linear', 'quadratic', or 'hybrid' TODO: Implement quadratic and hybrid


def calculate_bottom_stress_mannings_n(u, h, z, n):

    # Create iterator that will calculate Cd for each node
    cd = (g * n**2) / (h + z)**(1.0/3.0)

    # Apply lower limit
    np.clip(cd, cf, None, out=cd)

    # Calculate bottom friction based on Cd
    return calculate_bottom_stress_cd(u, cd)


def calculate_bottom_stress_cd(u, cd):

    if friction_type == 'linear':

        return cd * u

    if friction_type == 'quadratic':
        print('Quadtratic bottom friction not yet implemented. Falling back on linear.')
        return cd

    if friction_type == 'hybrid':
        print('Hybrid bottom friction not yet implemented. Falling back on linear.')
        return cd
