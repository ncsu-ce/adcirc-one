from scipy.integrate import quadrature as integrate


def phi_l(xl, xr):

    return lambda x: (xr - x) / (xr - xl)


def phi_r(xl, xr):

    return lambda x: (x - xl) / (xr - xl)


def dphi_l(xl, xr):

    return lambda x: -1 / (xr - xl)


def dphi_r(xl, xr):

    return lambda x: 1 / (xr - xl)


l = phi_l(0, 1)
r = phi_r(0, 1)
dl = dphi_l(0, 1)
dr = dphi_r(0, 1)

print(integrate(lambda x: dr(x)*dl(x), 0, 1))
