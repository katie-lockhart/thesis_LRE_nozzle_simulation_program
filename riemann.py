import numpy as np
import math
from scipy import optimize

from const import *


def flux_interpolate(rho_L,u_L,p_L,rho_R,u_R,p_R,gamma_L,gamma_R):
    f1 = 0.5*(rho_L * u_L + rho_R * u_R)
    f2 = 0.5*(rho_L * u_L**2.0 +p_L + rho_R * u_R**2.0 + p_R)
    E_L = p_L / (gamma_L-1.0) + 0.5*rho_L*u_L**2.0
    E_R = p_R / (gamma_R-1.0) + 0.5*rho_R*u_R**2.0
    f3 = 0.5*( u_L * (E_L + p_L) + u_R * (E_R + p_R) )
    return [f1, f2, f3]



def prefun(P,DK,PK,CK):
    if P <= PK:
        #  Rarefaction wave
        PRATIO = P/PK
        F    = G4*CK*(PRATIO**G1 - 1.0)
        FD   = (1.0/(DK*CK))*PRATIO**(-G2)
    else:
        # Shock wave
        AK  = G5/DK
        BK  = G6*PK
        QRT = (AK/(BK + P))**0.5
        F   = (P - PK)*QRT
        FD  = (1.0 - 0.5*(P - PK)/(BK + P))*QRT
    return [F,FD]


def func_p(x,rho_L,u_L,p_L,rho_R,u_R,p_R):
    P = x
    cc_L = (GAMMA*p_L/rho_L)**0.5
    cc_R = (GAMMA*p_R/rho_R)**0.5
    [FL, FD_L] = prefun(P,rho_L,p_L,cc_L)
    [FR, FD_R] = prefun(P,rho_R,p_R,cc_R)
    p_f = FL + FR + u_R - u_L
    return p_f


def ex_solver(rho_L,u_L,p_L,rho_R,u_R,p_R):
    # Speed of sound
    cc_L = (GAMMA*p_L/rho_L)**0.5
    cc_R = (GAMMA*p_R/rho_R)**0.5
    # assign first guess
    #p_guess = 0.5 * (p_R + p_L)
    CUP = 0.25*(rho_L + rho_R)*(cc_L + cc_R)
    PPV = 0.5*(p_L + p_R) + 0.5*(u_L - u_R)*CUP
    p_guess = max([1e-10,PPV])
    # solve
    p_star = optimize.newton(func_p, p_guess, args=(rho_L,u_L,p_L,rho_R,u_R,p_R), maxiter=2000)
    # compute velocity     
    [FL, FD_L] = prefun(p_star,rho_L,p_L,cc_L)
    [FR, FD_R] = prefun(p_star,rho_R,p_R,cc_R)
    u_star = 0.5*(u_L + u_R + FR - FL)
    return [p_star, u_star]

def sample(rho_L,u_L,p_L,rho_R,u_R,p_R,PM,UM,xxx,ttt):
    # Samples the solution at a given position relative to the
    # position of the initial discontinuity 
    # note that ttt (the time at which is sampled) should be positive (ttt>0)

    S = xxx/ttt # sampling at s=x/t

    DL = rho_L
    UL = u_L
    PL = p_L
    CL = (GAMMA*p_L/rho_L)**0.5
    
    DR = rho_R
    UR = u_R
    PR = p_R
    CR = (GAMMA*p_R/rho_R)**0.5

    if (S<=UM):
        # Sampling point lies to the left of the contact
        # discontinuity
        if (PM<=PL):
           # Left rarefaction
           SHL = UL - CL
           if (S<=SHL):
               # Sampled point is left data state
               D = DL
               U = UL
               P = PL
           else:
               CML = CL*(PM/PL)**G1
               STL = UM - CML
               if (S>STL):
                   # Sampled point is Star Left state
                   D = DL*(PM/PL)**(1.0/GAMMA)
                   U = UM
                   P = PM
               else:
                   # Sampled point is inside left fan
                   U = G5*(CL + G7*UL + S)
                   C = G5*(CL + G7*(UL - S))
                   D = DL*(C/CL)**G4
                   P = PL*(C/CL)**G3
        else:
            # Left shock
            PML = PM/PL
            SL  = UL - CL*(G2*PML + G1)**0.5
            if (S<=SL):
                # Sampled point is left data state
                D = DL
                U = UL
                P = PL
            else:
                # Sampled point is Star Left state
                D = DL*(PML + G6)/(PML*G6 + 1.0)
                U = UM
                P = PM
    else:
        # Sampling point lies to the right of the contact
        # discontinuity
        if (PM>PR):
            # Right shock
            PMR = PM/PR
            SR  = UR + CR*(G2*PMR + G1)**0.5
            if (S>=SR):
                # Sampled point is right data state
                D = DR
                U = UR
                P = PR
            else:
                # Sampled point is Star Right state
                D = DR*(PMR + G6)/(PMR*G6 + 1.0)
                U = UM
                P = PM
        else:
            # Right rarefaction
            SHR = UR + CR
            if (S>=SHR):
                # Sampled point is right data state
                D = DR
                U = UR
                P = PR
            else:
                CMR = CR*(PM/PR)**G1
                STR = UM + CMR
                if (S<=STR):
                    #  Sampled point is Star Right state
                    D = DR*(PM/PR)**(1.0/GAMMA)
                    U = UM
                    P = PM
                else:
                    # Sampled point is inside left fan
                    U = G5*(-CR + G7*UR + S)
                    C = G5*(CR - G7*(UR - S))
                    D = DR*(C/CR)**G4
                    P = PR*(C/CR)**G3


    return [U,D,P]


def get_flux(u_f,rho_f,p_f):
    f1 = rho_f*u_f
    f2 = rho_f*u_f**2.0+p_f 
    E_f = p_f / (GAMMA-1.0) + 0.5*rho_f*u_f**2.0
    f3 =  u_f * (E_f + p_f)
    return [f1, f2, f3]


def flux_RM_exact(rho_L,u_L,p_L,rho_R,u_R,p_R):
    [p_st, u_st] = ex_solver(rho_L,u_L,p_L,rho_R,u_R,p_R)
    [u_flux,rho_flux,p_flux] = sample(rho_L,u_L,p_L,rho_R,u_R,p_R,p_st, u_st,0.0,1.0)
    [f1, f2, f3] =  get_flux(u_flux,rho_flux,p_flux)
    return [f1, f2, f3]


 
if __name__ == "__main__":
    
    rho_L = 1.0
    u_L   = 0.0
    p_L   = 1.0
    
    rho_R = 0.125
    u_R   = 0.0
    p_R   = 0.1

    
    [p_st, u_st] = ex_solver(rho_L,u_L,p_L,rho_R,u_R,p_R)

    [u_flux,rho_flux,p_flux] = sample(rho_L,u_L,p_L,rho_R,u_R,p_R,p_st, u_st,0.0,1.0)

    [f1, f2, f3] =  get_flux(u_flux,rho_flux,p_flux)

    [f1_, f2_, f3_] = flux_RM_exact(rho_L,u_L,p_L,rho_R,u_R,p_R)

    print(p_st,u_st)
    print(u_flux,rho_flux,p_flux)
    print(f1, f2, f3)
    print(f1_, f2_, f3_)


 


