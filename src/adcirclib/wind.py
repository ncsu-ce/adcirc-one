import numpy as np


wind_drag_limit = 0.0035


def drag_coefficient(u):

        cd = 0.001 * (0.75 + 0.067 * u)
        np.clip(cd, None, wind_drag_limit, out=cd)

        return cd
