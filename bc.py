import numpy as np
import math
from scipy import optimize

from const import *

# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------
# SUBSONIC INLET
def get_OS_inlet_subsonic(p_s,T_O, p_O):
    rho_O = p_O / (R_gas_const * T_O)
    rho_O_s = rho_O * (p_s/p_O)**(1/GAMMA)
    T_O_s = p_s/(R_gas_const * rho_O_s)
    return [rho_O_s,T_O_s]

def get_IS_inlet_subsonic(p_s,rho_I, p_I, u_I):
    rho_I_s = rho_I * (p_s/p_I)**(1/GAMMA)
    a_I_s = (GAMMA*p_s/rho_I_s)**0.5
    a_I   = (GAMMA*p_I/rho_I  )**0.5
    u_s = u_I + (a_I_s-a_I)/delta
    return [rho_I_s,a_I_s,a_I,u_s]


def f_p_inlet_subsonic(x, T_O, p_O, rho_I, p_I, u_I):
    p_s = x
    [rho_O_s,T_O_s] = get_OS_inlet_subsonic(p_s,T_O, p_O)
    [rho_I_s,a_I_s,a_I,u_s] = get_IS_inlet_subsonic(p_s,rho_I, p_I, u_I)
    f_i_s =  cp * (T_O_s - T_O) + u_s**2.0/2.0
    return f_i_s


def get_p_inlet_subsonic(T_O, p_O, rho_I, p_I, u_I):
    p_guess = 0.5*(p_O + p_I)
    p_star = optimize.newton(f_p_inlet_subsonic, p_guess, args=(T_O, p_O, rho_I, p_I, u_I))
    return p_star


def get_flux_inlet_subsonic(T_O, p_O, rho_I, p_I, u_I):
    p_s = get_p_inlet_subsonic(T_O, p_O, rho_I, p_I, u_I)
    [rho_O_s,T_O_s] = get_OS_inlet_subsonic(p_s,T_O, p_O)
    [rho_I_s,a_I_s,a_I,u_s] = get_IS_inlet_subsonic(p_s,rho_I, p_I, u_I)
    f1 = rho_O_s*u_s
    f2 = rho_O_s*u_s**2.0+p_s
    E_f = p_s / (GAMMA-1.0) + 0.5*rho_O_s*u_s**2.0
    f3 =  u_s * (E_f + p_s)
    print('bc subsonic inlet: u= ',u_s,'p= ',p_s )
    return [f1, f2, f3]


def get_flux_outlet_subsonic(p_R_out, rho_L_out, p_L_out, u_L_out):
    p_s_out = p_R_out
    rho_s_L_out = rho_L_out * ((p_s_out / p_L_out) ** (1/GAMMA))
    a_L_out = ((GAMMA * p_L_out) / rho_L_out) ** 0.5
    a_s_L_out = a_L_out * ((p_s_out / p_L_out) ** (G1))
    u_s_out = u_L_out + ((a_L_out - a_s_L_out) / delta)

    f1_out = rho_s_L_out * u_s_out
    f2_out = rho_s_L_out * u_s_out ** 2.0 + p_s_out
    E_f_out = p_s_out / (GAMMA-1.0) + 0.5 * rho_s_L_out * u_s_out ** 2.0
    f3_out =  u_s_out * (E_f_out + p_s_out)
    print('bc subsonic outlet: u_st = ', u_s_out, 'p_st = ', p_s_out, 'rho_st = ', rho_s_L_out,
          'a_L = ', a_L_out, 'a_st = ', a_s_L_out)
    return [f1_out, f2_out, f3_out]


# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------

 
if __name__ == "__main__":
    # this is a test program, not exectuted when the module is imported
    # in another script
    
    rho_O = 1
    p_O = 2
    T_O = p_O/(R_gas_const * rho_O)
    rho_I = 1
    p_I = 2
    u_I = 1
    p_s = get_p_inlet_subsonic(T_O, p_O, rho_I, p_I, u_I)
    [rho_O_s,T_O_s] = get_OS_inlet_subsonic(p_s,T_O, p_O)
    [rho_I_s,a_I_s,a_I,u_s] = get_IS_inlet_subsonic(p_s,rho_I, p_I, u_I)

    [f1, f2, f3] = get_flux_inlet_subsonic(T_O, p_O, rho_I, p_I, u_I)


    p_R_out = 101325
    rho_L_out = 0.3
    p_L_out = 200000
    u_L_out = 200

    [f1_out, f2_out, f3_out] = get_flux_outlet_subsonic(p_R_out, rho_L_out, p_L_out, u_L_out)


 


