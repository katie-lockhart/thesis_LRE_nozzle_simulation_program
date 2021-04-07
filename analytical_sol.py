from const import *
from scipy import optimize


def area_mach_relation(x, area_ratio):
    mach = x
    f_mach = (1 / mach ** 2) * ((G5 * (1 + (G7 * mach ** 2))) ** G9) - area_ratio ** 2
    return f_mach


def get_mach(area_ratio, mach_guess):
    mach = optimize.newton(area_mach_relation, mach_guess, args=(area_ratio,), maxiter=200)
    return mach


def get_vars(local_mach, p_stag, T_stag, rho_stag):
    p_analytical_quasi = p_stag / ((1 + (G7 * local_mach**2)) ** G10)
    T_analytical_quasi = T_stag / (1 + (G7 * local_mach**2))
    rho_analytical_quasi = rho_stag / ((1 + (G7 * local_mach**2)) ** G11)
    return [p_analytical_quasi, T_analytical_quasi, rho_analytical_quasi]


def solve_analytical(area_ratio, p_stag, T_stag, rho_stag, mach_guess):
    local_mach = get_mach(area_ratio, mach_guess)
    [p_analytical_quasi, T_analytical_quasi, rho_analytical_quasi] = get_vars(local_mach, p_stag, T_stag, rho_stag)
    return [local_mach, p_analytical_quasi, T_analytical_quasi, rho_analytical_quasi]


if __name__ == "__main__":

    mach_guess = 0.00001
    area_ratio = 7.2551
    p_stag = 1e+6
    T_stag = 2500.0
    rho_stag = 1.393

    local_mach = get_mach(area_ratio, mach_guess)

    [p_analytical_quasi, T_analytical_quasi, rho_analytical_quasi] = get_vars(local_mach, p_stag, T_stag, rho_stag)


    print('local mach = ', local_mach)
    print('analytical p = ', p_analytical_quasi)
    print('analytical T = ', T_analytical_quasi)
    print('analytical rho = ', rho_analytical_quasi)
