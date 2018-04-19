import numpy as np


g = 9.81                    # Gravity
cf = 1                      # Bottom friction coefficient used unless otherwise specified TODO: Set this value
friction_type = 'quadratic'    # Can be 'linear', 'quadratic', or 'hybrid' TODO: Implement quadratic and hybrid


def calculate_bottom_stress_mannings_n(q, h, z, n):

    # Calculate Cd for each node
    cd = (g * n**2) / (h + z)**(1.0/3.0)

    # Apply lower limit
    np.clip(cd, cf, None, out=cd)

    # Calculate bottom friction based on Cd
    return calculate_bottom_stress_cd(q, h, cd)


def calculate_bottom_stress_cd(q, h, cd):

    if friction_type == 'linear':

        return cd

    if friction_type == 'quadratic':

        return cd * (q / h)

    if friction_type == 'hybrid':
        print('Hybrid bottom friction not yet implemented. Falling back on linear.')
        return cd

    print('ERROR: Bottom friction method not specified.')
