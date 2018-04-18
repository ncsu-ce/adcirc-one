from .wind import drag_coefficient


rho_air = 0.001293


def calculate_wind_stress(u_wind):

    cd = drag_coefficient(u_wind)
    return rho_air * cd * u_wind ** 2
