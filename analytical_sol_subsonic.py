from const import *
from scipy import optimize
from sympy import symbols, Eq, solve

def get_outlet_mach(p_stag, p_outlet):
    p_ratio = p_stag/p_outlet
    p_ratio_exp = p_ratio ** (G12)
    mach_squared = (p_ratio_exp - 1) / G7
    mach = np.sqrt(mach_squared)
    return mach


def get_area_star(mach_outlet, area_outlet):
    RHS = np.sqrt((1 / (mach_outlet ** 2)) * ((G5 * (1 + (G7 * (mach_outlet ** 2)))) ** G9))
    a_star = area_outlet / RHS

    return a_star

def area_mach_relation(x, area_ratio):
    mach = x
    f_mach = (1 / mach ** 2) * ((G5 * (1 + (G7 * mach ** 2))) ** G9) - area_ratio ** 2
    return f_mach

def get_mach(area_ratio, mach_guess):
    mach = optimize.newton(area_mach_relation, mach_guess, args=(area_ratio,), maxiter=2000)
    return mach


def get_vars(local_mach, p_stag, T_stag, rho_stag):
    p_analytical_quasi = p_stag / ((1 + (G7 * local_mach**2)) ** G10)
    T_analytical_quasi = T_stag / (1 + (G7 * local_mach**2))
    rho_analytical_quasi = rho_stag / ((1 + (G7 * local_mach**2)) ** G11)
    return [p_analytical_quasi, T_analytical_quasi, rho_analytical_quasi]


def solve_subsonic_analytical(p_stag, p_outlet, area_outlet, area, T_stag, rho_stag, mach_guess):
    mach_outlet = get_outlet_mach(p_stag, p_outlet)
    area_star = get_area_star(mach_outlet, area_outlet)

    area_ratio = (area / area_star)

    local_mach_sub = get_mach(area_ratio, mach_guess)
    [p_analytical_sub, T_analytical_sub, rho_analytical_sub] = get_vars(local_mach_sub, p_stag, T_stag, rho_stag)

    return [local_mach_sub, p_analytical_sub, T_analytical_sub, rho_analytical_sub]




if __name__ == "__main__":

    p_outlet = 101235
    p_stag = p_outlet*1.01
    area_outlet = 0.22
    area = 0.075
    T_stag = 2000.0
    rho_stag = 1.393
    mach_guess = 0.001

    mach_outlet = get_outlet_mach(p_stag, p_outlet)

    area_star = get_area_star(mach_outlet, area_outlet)

    print('mach_outlet = ', mach_outlet)
    print('area_star', area_star)

    area_ratio = (area / area_star)

    local_mach_sub = get_mach(area_ratio, mach_guess)

    print('local mach = ', local_mach_sub)

    [p_analytical_sub, T_analytical_sub, rho_analytical_sub] = get_vars(local_mach_sub, p_stag, T_stag, rho_stag)


    print('analytical p = ', p_analytical_sub)
    print('analytical T = ', T_analytical_sub)
    print('analytical rho = ', rho_analytical_sub)