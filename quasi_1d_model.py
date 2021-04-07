import numpy as np
import matplotlib.pyplot as plt
import riemann
import math
import bc
from const import *
import analytical_sol
import analytical_sol_subsonic
import parabolic_nozzle
from scipy.interpolate import interp1d
import os

def solve_flow(inlet_pressure, inlet_temperature, outlet_pressure, nozzle_geometry, boundaries, sample_exact, heat,
               coolant_wall_temp, T_w_hg, viscosity_fluid, cp_fluid, k_fluid, t_max, t_steady, init, init_rho_iter,
               init_u_iter, init_p_iter, file_no, x_domain, x_min, x_max):

    # useful function to pause the code
    def pause():
        programPause = input("Press the <ENTER> key to continue...")


    # ------------------------------------------
    # Create output folder
    folder_out_0 = 'output_full_solution_iter'      #full solution at each global iteration
    if not os.path.exists(folder_out_0):
        os.makedirs(folder_out_0)

    folder_out = 'output_subsonic_51'

    if not os.path.exists(folder_out):
        os.makedirs(folder_out)
    # ------------------------------------------

    # ------------------------------------------
    # geometry definition and domain discretization
    nx = 51


    if nozzle_geometry == 0: #rao nozzle
        throat_radius = 0.15
        ratio_area = 3
        bell_percent = 0.8

        [nozzle_length, nozzle_x_coord, nozzle_y_coord] = parabolic_nozzle.rao_nozzle(throat_radius, ratio_area, bell_percent)

        L = np.abs(max(nozzle_x_coord) - min(nozzle_x_coord))
        x = np.linspace(min(nozzle_x_coord), max(nozzle_x_coord), nx)
        x_min = min(nozzle_x_coord)
        x_max = max(nozzle_x_coord)


    if nozzle_geometry == 1: #define any radius function
        L = 1
        x = np.linspace(x_min, x_max, nx)

    if nozzle_geometry == 2: #constant area
        L = x_domain
        x = np.linspace(x_min, x_max, nx)
        print('L = ', L)

    # ------------------------------------------
    # define constants
    # numerical

    dt_out = 0.02

    sigma = 1
    CFL = 0.5

    dx = L / (nx - 1)
    dt = sigma * dx

    # ------------------------------------------

    # ------------------------------------------
    # set up arrays
    rho = np.ones(nx)
    u = np.zeros(nx)
    p = np.ones(nx)
    E = np.ones(nx)
    T = np.ones(nx)

    # old variables
    u_old = np.zeros(nx)
    p_old = np.zeros(nx)
    rho_old = np.zeros(nx)


    radius_nozzle = np.ones(nx)
    Area_nozzle_centred = np.ones(nx)
    dA_dx_nozzle_centred = np.ones(nx) #centred means at cell centre
    area_ratio = np.ones(nx)
    A_throat = np.ones(nx)
    radius_throat = np.ones(nx)
    dia_hydraulic = np.ones(nx)
    dia_nozzle = np.ones(nx)
    dia_throat = np.ones(nx)

    #print('x = ', x)

    W_exhaust = np.zeros((nx, 3))
    F_exhaust = np.zeros((nx + 1, 3))
    Q_exhaust = np.zeros((nx, 3))

    rho_exact = np.ones(nx)
    u_exact = np.zeros(nx)
    p_exact = np.ones(nx)

    rho_analytical_quasi = np.ones(nx)
    u_analytical_quasi = np.zeros(nx)
    p_analytical_quasi = np.ones(nx)
    T_analytical_quasi = np.ones(nx)
    mach_analytical_quasi = np.ones(nx)

    mach_analytical_sub = np.ones(nx)
    p_analytical_sub = np.ones(nx)
    T_analytical_sub = np.ones(nx)
    rho_analytical_sub = np.ones(nx)
    u_analytical_sub = np.ones(nx)

    u_sol_diff = np.ones(nx)
    p_sol_diff = np.ones(nx)
    mach_sol_diff = np.ones(nx)
    T_sol_diff = np.ones(nx)
    rho_sol_diff = np.ones(nx)

    heat_flux = np.ones(nx)
    h_coeff_exhaust = np.ones(nx)
    T_aw_exhaust = np.ones(nx)
    T_w_exhaust = np.ones(nx)
    recover_fac = np.ones(nx)

    T_w_coolant = np.ones(nx)
    wall_conductivity = np.ones(nx)
    # ------------------------------------------
    # compute geometry parameters

    if nozzle_geometry == 0:

        radius_function = interp1d(nozzle_x_coord, nozzle_y_coord)
        for i in range(0, nx):


            radius_nozzle[i] = radius_function(x[i])
            dia_hydraulic[i] = 2 * radius_nozzle[i]
            dia_nozzle[i] = 2 * radius_nozzle[i]
            Area_nozzle_centred[i] = math.pi * radius_nozzle[i]**2

            file_name_geo = folder_out_0 + '/data_' + 'geometry' + '.dat'  # file for residuals
            f_geo = open(file_name_geo, 'w')

            for i in range(0, nx):
                s_tot = str(x[i]) + '  ' + str(Area_nozzle_centred[i]) + '  ' + str(radius_nozzle[i])
                f_geo.write(s_tot + '\n')

            f_geo.close()

        A_throat = math.pi * throat_radius **2
        radius_throat = np.min(radius_nozzle)
        dia_throat = radius_throat * 2

        for i in range(0, nx):
            area_ratio[i] = Area_nozzle_centred[i] / A_throat

        print('area ratio = ', area_ratio)

        for i in range(1, nx - 1):
            dA_dx_nozzle_centred[i] = (Area_nozzle_centred[i + 1] - Area_nozzle_centred[i - 1]) / (2 * dx)

        print('dA_dx_nozzle_centred = ', dA_dx_nozzle_centred)


    if nozzle_geometry == 1:
        for i in range(0, nx):

            radius_nozzle[i] = 0.2 + 0.01 * math.cos(2 * math.pi * x[i])
            dia_hydraulic[i] = 2 * radius_nozzle[i]
            dia_nozzle[i] = 2 * radius_nozzle[i]
            Area_nozzle_centred[i] = math.pi * radius_nozzle[i]**2

        A_throat = np.min(Area_nozzle_centred)
        radius_throat = np.min(radius_nozzle)
        dia_throat = radius_throat * 2

        for i in range(0, nx):
            area_ratio[i] = Area_nozzle_centred[i] / A_throat

        print('area ratio = ', area_ratio)

        for i in range(1, nx - 1):
            dA_dx_nozzle_centred[i] = (Area_nozzle_centred[i + 1] - Area_nozzle_centred[i - 1]) / (2 * dx)

        print('dA_dx_nozzle_centred = ', dA_dx_nozzle_centred)


    if nozzle_geometry == 2:
        for i in range(0, nx):
            Area_nozzle_centred[i] = 0.000127
            radius_nozzle[i] = np.sqrt(Area_nozzle_centred[i] / math.pi)
            dia_hydraulic[i] = radius_nozzle[i] * 2
            dia_nozzle[i] = 2 * radius_nozzle[i]

        A_throat = np.min(Area_nozzle_centred)

        for i in range(0, nx):
            area_ratio[i] = Area_nozzle_centred[i] / A_throat

        print('area ratio = ', area_ratio)

        for i in range(0, nx):
            dA_dx_nozzle_centred[i] = 0

        print('dA_dx_nozzle_centred = ', dA_dx_nozzle_centred)


    throat_pos = np.argmin(Area_nozzle_centred)

    print('throat position = ', throat_pos)


    # ------------------------------------------
    # Initial Conditions

    if init == 0:
        init_rho_L = 1
        init_u_L = 0.2
        init_p_L = 101325

        init_rho_R = 1
        init_u_R = 0.2
        init_p_R = 101325

        for i in range(0, nx):

            if (x[i] < 0.5):
                rho[i] = init_rho_L
                u[i] = init_u_L
                p[i] = init_p_L
            else:
                rho[i] = init_rho_R
                u[i] = init_u_R
                p[i] = init_p_R

            E[i] = p[i] / (gamma - 1) + 0.5 * rho[i] * u[i] ** 2
            # E = total energy per unit volume

        a_squared = gamma * p / rho
        a = np.sqrt(a_squared)

    if init == 1:

        init_rho_cool = init_rho_iter
        init_u_cool = init_u_iter
        init_p_cool = init_p_iter

        for i in range(0, nx):

            rho[i] = init_rho_cool[i]
            u[i] = init_u_cool[i]
            p[i] = init_p_cool[i]

            E[i] = p[i] / (gamma - 1) + 0.5 * rho[i] * u[i] ** 2
            # E = total energy per unit volume

        a_squared = gamma * p / rho
        a = np.sqrt(a_squared)

    # ------------------------------------------


    # ------------------------------------------
    # Plot initial conditions
    plt.figure(0)
    plt.plot(x, rho, '-')
    plt.plot(x, u, '-')
    plt.plot(x, p, '-')
    plt.plot(x, E, '-')
    plt.xlabel('x')
    plt.ylabel('vars')

   # plt.show()

    # ------------------------------------------

    # ------------------------------------------
    # Boundary Conditions

    # LEFT
    # bc_left = 'zero-grad'

    # bc_left = 'supersonic'
    p_in = 1
    rho_in = 1.0
    u_in = 2

    if boundaries == 0:  # subsonic in, supersonic out
        bc_left = 'subsonic'
        p_O_inlet = inlet_pressure
        # p_O_inlet = 1.01*101325
        T_O_inlet = inlet_temperature
        # rho_O_inlet = 1
        rho_O_inlet = p_O_inlet / (R_gas_const * T_O_inlet)
        # T_O_inlet = p_O_inlet / (rho_O_inlet * R_gas_const)

        print('rho_inlet = ', rho_O_inlet)

        # RIGHT
        bc_right = 'zero-grad'

    if boundaries == 1:  # subsonic in, subsonic out
        #bc_left = 'zero_grad'
        bc_left = 'subsonic'
        p_O_inlet = inlet_pressure
        # p_O_inlet = 1.01*101325
        T_O_inlet = inlet_temperature
        # rho_O_inlet = 1
        rho_O_inlet = p_O_inlet / (R_gas_const * T_O_inlet)
        # T_O_inlet = p_O_inlet / (rho_O_inlet * R_gas_const)

        bc_right = 'subsonic'
        #bc_right = 'zero-grad'
        p_out = outlet_pressure
    # ------------------------------------------

    p_throat = p_O_inlet * (G5 ** G10)
    rho_throat = rho_O_inlet * (G5 ** G11)
    T_throat = p_throat / (rho_throat * R_gas_const)
    # ------------------------------------------


    # ------------------------------------------
    # Plot Area
    plt.figure(1)
    plt.plot(x, Area_nozzle_centred, '-',x,radius_nozzle , '.')
    plt.xlabel('x')
    plt.ylabel('A, r')
    # ------------------------------------------

    # ------------------------------------------
    # open file for global results
    file_name = folder_out_0 + '/data_' + 'residual' + '_' + str(file_no) + '.dat' #file for residuals
    f_io = open(file_name,'w')
    # ------------------------------------------


    # ------------------------------------------
    # MAIN TIME LOOP
    t = 0
    t_out = 0
    t_file_res = 0

    while (t < t_max):

        # ----------------------------------------------------
        # Store old variables
        u_old = np.copy(u)
        p_old = np.copy(p)
        rho_old = np.copy(rho)
        # ----------------------------------------------------
        # Heat source term

        Mach = u / a
        perimeter_wall = (4 * Area_nozzle_centred) / dia_hydraulic

        if heat == 0:

            mass_flux = rho * u * Area_nozzle_centred

            characteristic_velocity = (inlet_pressure * A_throat) / mass_flux

            Re_exhaust = rho * u * dia_nozzle / viscosity_fluid
            Pr_exhaust = viscosity_fluid * cp_fluid / k_fluid  # stagnation conds
            nusselt_number_exhaust = 0.026 * (Re_exhaust ** 0.8) * (Pr_exhaust ** 0.4)

            recover_fac = Pr_exhaust ** 0.3

            bl_corrector_exponent = 0.6 / 5 #from Bartz

            boundary_layer_correction = 1 / ((((0.5 * (T_w_hg / inlet_temperature) * (1 * delta * (Mach**2))) + 0.5)
                                              ** (0.8-bl_corrector_exponent)) * ((1 + (delta * (Mach**2)))**bl_corrector_exponent))

            h_coeff_exhaust = (0.026 / (dia_throat**0.2)) * ((viscosity_fluid**0.2 * cp_fluid) / (Pr_exhaust**0.6)) * \
                ((inlet_pressure / characteristic_velocity)**0.8) * ((A_throat / Area_nozzle_centred)**0.9) * boundary_layer_correction


            #h_coeff_exhaust = (nusselt_number_exhaust * k_fluid) / dia_nozzle
            #h_coeff_exhaust = (rho * u) ** 0.8
            #h_coeff_exhaust = 800

            #T_aw = 3000
            T_aw = T * (1 + (recover_fac * ((gamma - 1) / 2) * Mach ** 2))
            T_w = T_w_hg

            heat_flux = h_coeff_exhaust * (T_w - T_aw) # negative source for heat leaving exhaust

            heat_flux_return = heat_flux * -1 # return positive heat flux


        if heat == 1:

            Re_exhaust = (rho * u * dia_nozzle) / viscosity_fluid
            Pr_exhaust = (viscosity_fluid * cp_fluid) / k_fluid

            nusselt_number_coolant = 0.023 * (Re_exhaust**0.8) * (Pr_exhaust**0.4)

            h_coeff_exhaust = (nusselt_number_coolant * k_fluid) / dia_nozzle

            T_aw = T * (1 + (recover_fac * ((gamma - 1) / 2) * Mach ** 2))

            T_w = coolant_wall_temp

            heat_flux = h_coeff_exhaust * (T_w - T) # positive for heat entering coolant

            #print('convective coeff coolant = ', h_coeff_exhaust)

            heat_flux_return = heat_flux

        if heat == 2:

            T_aw = T * (1 + (recover_fac * ((gamma - 1) / 2) * Mach ** 2))
            heat_flux = np.zeros(nx)

            #print('convective coeff coolant = ', h_coeff_exhaust)

            heat_flux_return = heat_flux

        # ----------------------------------------------------

        # ----------------------------------------------------

        darcy_friction = 0.005

        shear_stress_wall = (1 / 8) * rho * (u ** 2) * darcy_friction
        friction_source_momentum = ((1 / Area_nozzle_centred) * perimeter_wall * shear_stress_wall)

        # ----------------------------------------------------
        W_exhaust[:, 0] = rho
        W_exhaust[:, 1] = rho * u
        W_exhaust[:, 2] = E
        # ----------------------------------------------------

        Q_exhaust[:, 0] = (-1 / Area_nozzle_centred) * (dA_dx_nozzle_centred * rho * u)
        Q_exhaust[:, 1] = (-1 / Area_nozzle_centred) * (dA_dx_nozzle_centred * rho * u ** 2) - \
                          ((1 / Area_nozzle_centred) * perimeter_wall * shear_stress_wall)
        Q_exhaust[:, 2] = (-1 / Area_nozzle_centred) * (dA_dx_nozzle_centred * u * (E + p)) + \
                          ((1 / Area_nozzle_centred) * perimeter_wall * heat_flux)

        # ----------------------------------------------------
        # Compute dt according to eigenvalues and CFL
        max_u = np.max(np.abs(u) + a)
        dt = (CFL * dx) / max_u
        t = t + dt
        t_out = t_out + dt

        print('----------------------------------------------')
        print('t= ', t, '    dt= ', dt, '    t_out=', t_out)
        # ----------------------------------------------------


        # ----------------------------------------------------
        # Riemann problem at each interface and flux computation (upwind)
        for i in range(1, nx):

            rho_L = rho[i - 1]
            u_L = u[i - 1]
            p_L = p[i - 1]
            rho_R = rho[i]
            u_R = u[i]
            p_R = p[i]
            [f1, f2, f3] = riemann.flux_RM_exact(rho_L, u_L, p_L, rho_R, u_R,p_R)

            F_exhaust[i, 0] = f1
            F_exhaust[i, 1] = f2
            F_exhaust[i, 2] = f3
        # ----------------------------------------------------


        # ----------------------------------------------------
        #Left BCs
        dA_dx_nozzle_centred[0] = dA_dx_nozzle_centred[1]

        if bc_left == 'zero-grad':
            F_exhaust[0, 0] = rho[0] * u[0]
            F_exhaust[0, 1] = rho[0] * u[0] ** 2 + p[0]
            F_exhaust[0, 2] = u[0] * (E[0] + p[0])

        if bc_left == 'supersonic':
            E_in = p_in / (gamma - 1) + 0.5 * rho_in * u_in ** 2
            F_exhaust[0, 0] = rho_in * u_in
            F_exhaust[0, 1] = rho_in * u_in ** 2 + p_in
            F_exhaust[0, 2] = u_in * (E_in + p_in)

        if bc_left == 'subsonic':
            [f1, f2, f3] = bc.get_flux_inlet_subsonic(T_O_inlet, p_O_inlet, rho[0], p[0], u[0])
            F_exhaust[0, 0] = f1
            F_exhaust[0, 1] = f2
            F_exhaust[0, 2] = f3
        # ----------------------------------------------------


        # ----------------------------------------------------
        #Right BCs
        dA_dx_nozzle_centred[nx-1] = dA_dx_nozzle_centred[nx-2]

        if bc_right == 'zero-grad':
            F_exhaust[nx, 0] = rho[nx-1] * u[nx-1]
            F_exhaust[nx, 1] = rho[nx-1] * u[nx-1] ** 2 + p[nx-1]
            F_exhaust[nx, 2] = u[nx-1] * (E[nx-1] + p[nx-1])

        if bc_right == 'subsonic':
            [f1_out, f2_out, f3_out] = bc.get_flux_outlet_subsonic(p_out, rho[nx-1], p[nx-1], u[nx-1])
            F_exhaust[nx, 0] = f1_out
            F_exhaust[nx, 1] = f2_out
            F_exhaust[nx, 2] = f3_out
        # ----------------------------------------------------


        # ----------------------------------------------------
        # Conserved variable update with Godunov method
        for i in range(0, nx):
            W_exhaust[i, 0] = W_exhaust[i, 0] - dt / dx * (F_exhaust[i + 1, 0] - F_exhaust[i, 0]) + dt * Q_exhaust[i, 0]
            W_exhaust[i, 1] = W_exhaust[i, 1] - dt / dx * (F_exhaust[i + 1, 1] - F_exhaust[i, 1]) + dt * Q_exhaust[i, 1]
            W_exhaust[i, 2] = W_exhaust[i, 2] - dt / dx * (F_exhaust[i + 1, 2] - F_exhaust[i, 2]) + dt * Q_exhaust[i, 2]
        # ----------------------------------------------------

        # ----------------------------------------------------
        # Computation of primitive and additional vars
        rho = W_exhaust[:, 0]
        u = W_exhaust[:, 1] / rho
        E = W_exhaust[:, 2]
        p = (gamma - 1.0) * (E - 0.5 * rho * u ** 2)
        T = p / (rho * R_gas_const)
        E_specific = p / ((gamma - 1) * rho)

        a_squared = gamma * p / rho
        a = np.sqrt(a_squared)

        Mach = u / a
        print('Max Mach= ',max(Mach))
        # ----------------------------------------------------

        # ----------------------------------------------------
        # Old variables to be returned. Store as full solution for this iteration
        rho_return = rho
        u_return = u
        E_return = E
        p_return = p
        T_return = T
        mass_flow = rho*u*Area_nozzle_centred

        file_name= folder_out_0 + '/data_steady_state' + str(nx) + '_' + str(file_no) + '.dat'
        f_fullsol = open(file_name,'w')

        for i in range(0, nx):
            s_tot = str(x[i]) + '  ' + str(u[i]) + '  ' + str(p[i]) + '  ' + str(rho[i]) + '  ' + str(T[i]) + '  ' + str(Mach[i]) + '  ' + str(mass_flow[i])
            f_fullsol.write(s_tot + '\n')

        f_fullsol.close()
        # -------------------

        # ----------------------------------------------------
        # Compute residual
        u_res = np.sum(np.abs(u-u_old))
        p_res = np.sum(np.abs(p-p_old))
        rho_res = np.sum(np.abs(rho-rho_old))
        s1 = str(t)
        s2 = str(u_res)
        s3 = str(p_res)
        s4 = str(rho_res)
        s_tot = s1 + ' ' + s2 + ' ' + s3  + ' ' + s4  + ' '
        f_io.write(s_tot + '\n')   #into output.dat file
        print('Residual u, p, rho= ',u_res, p_res, rho_res)

        t_file_res = 0
        #if residuals all come within some tolerance, t_steady = t

        # ----------------------------------------------------

        # ----------------------------------------------------
        # Dump solution on file
        if t_out >= dt_out:
            file_name= folder_out + '/data_' + str(nx) + '_' + str(t) + '.dat'
            f_sol = open(file_name,'w')
            for i in range(0, nx):
                s_tot = str(x[i]) + '  ' + str(u[i]) + '  ' + str(p[i]) + '  ' + str(rho[i]) + '  ' + str(T[i]) + '  ' + str(Mach[i])
                f_sol.write(s_tot + '\n')
            f_sol.close()
        # ----------------------------------------------------

        # ----------------------------------------------------
        # Plot output
        # Analytical solutions
        if t_out >= dt_out:
            t_out = 0

            # ----------------------------------------------------
            if sample_exact == 1:  # riemann exact solution - steady or unsteady 1D flow
                [p_st, u_st] = riemann.ex_solver(init_rho_L, init_u_L, init_p_L, init_rho_R, init_u_R, init_p_R)
                for i in range(0, nx):
                    [u_exact[i], rho_exact[i], p_exact[i]] = riemann.sample(init_rho_L, init_u_L, init_p_L, init_rho_R,
                                                                            init_u_R, init_p_R, p_st, u_st, x[i] - 0.5, t)

            # ----------------------------------------------------
            if sample_exact == 2: #subsonic/supersonic transitional flow analytical solution quasi 1D
                if t > t_steady:
                    for i in range(0, nx):

                        if i < throat_pos: #POSITION of throat
                            mach_guess = 0.1
                            [mach_analytical_quasi[i], p_analytical_quasi[i], T_analytical_quasi[i], rho_analytical_quasi[i]] = \
                                analytical_sol.solve_analytical(area_ratio[i], p_O_inlet, T_O_inlet, rho_O_inlet, mach_guess)

                            u_analytical_quasi = mach_analytical_quasi * a

                        else:
                            mach_guess = 20
                            [mach_analytical_quasi[i], p_analytical_quasi[i], T_analytical_quasi[i], rho_analytical_quasi[i]] = \
                                analytical_sol.solve_analytical(area_ratio[i], p_O_inlet, T_O_inlet, rho_O_inlet, mach_guess)

                            u_analytical_quasi = mach_analytical_quasi * a

                        u_sol_diff[i] = abs(u[i] - u_analytical_quasi[i])
                        p_sol_diff[i] = abs(p[i] - p_analytical_quasi[i])
                        T_sol_diff[i] = abs(T[i] - T_analytical_quasi[i])
                        mach_sol_diff[i] = abs(Mach[i] - mach_analytical_quasi[i])
                        rho_sol_diff[i] = abs(rho[i] - rho_analytical_quasi[i])

                    # Dump analytical solution on file
                    file_name= folder_out + '/exact_sol_quasi1d' + str(nx) + '.dat'
                    f_ana = open(file_name,'w')
                    for i in range(0, nx):
                        s_tot = str(x[i]) + '  ' + str(u_analytical_quasi[i]) + '  ' + str(p_analytical_quasi[i]) + '  ' + \
                                str(rho_analytical_quasi[i]) + '  ' + str(T_analytical_quasi[i]) + '  ' + str(mach_analytical_quasi[i])
                        f_ana.write(s_tot + '\n')
                    f_ana.close()

            # ----------------------------------------------------
            if sample_exact == 3: #subsonic inlet and outlet exact solution quasi 1D
                if t > t_steady:
                    for i in range(0, nx):
                        mach_guess = 0.001
                        [mach_analytical_sub[i], p_analytical_sub[i], T_analytical_sub[i], rho_analytical_sub[i]] = \
                                analytical_sol_subsonic.solve_subsonic_analytical(p_O_inlet, p_out, Area_nozzle_centred[nx-1],
                                        Area_nozzle_centred[i], T_O_inlet, rho_O_inlet, mach_guess)

                        u_analytical_sub = mach_analytical_sub * a

                    # Dump analytical solution on file
                    file_name= folder_out + '/exact_sol_1d' + str(nx) + '.dat'
                    f_ana = open(file_name,'w')
                    for i in range(0, nx):
                        s_tot = str(x[i]) + '  ' + str(u_analytical_sub[i]) + '  ' + str(p_analytical_sub[i]) + '  ' + \
                                str(rho_analytical_sub[i]) + '  ' + str(T_analytical_sub[i]) + '  ' + str(mach_analytical_sub[i])
                        f_ana.write(s_tot + '\n')
                    f_ana.close()

                    print('mach subsonic ana = ', mach_analytical_sub)
            # ----------------------------------------------------
        # ----------------------------------------------------

            plt.figure(11)
            plt.plot(x, p, '.')
            if sample_exact == 1:
                plt.plot(x, p_exact, '-')
            plt.xlabel('x')
            plt.ylabel('p')
            if sample_exact == 2:
                if t > t_steady:
                    plt.plot(x, p_analytical_quasi, '-')
            if sample_exact == 3:
                if t > t_steady:
                    plt.plot(x, p_analytical_sub, '-')
            plt.xlabel('x')
            plt.ylabel('p')

            plt.figure(12)
            plt.plot(x, u, '.')
            if sample_exact == 1:
                plt.plot(x, u_exact, '-')
            plt.xlabel('x')
            plt.ylabel('u')
            if sample_exact == 2:
                if t > t_steady:
                    plt.plot(x, u_analytical_quasi, '-')
            if sample_exact == 3:
                if t > t_steady:
                    plt.plot(x, u_analytical_sub, '-')
            plt.xlabel('x')
            plt.ylabel('u')

            plt.figure(13)
            plt.plot(x, rho, '.')
            if sample_exact == 1:
                plt.plot(x, rho_exact, '-')
            plt.xlabel('x')
            plt.ylabel('rho')
            if sample_exact == 2:
                if t > t_steady:
                    plt.plot(x, rho_analytical_quasi, '-')
            if sample_exact == 3:
                if t > t_steady:
                    plt.plot(x, rho_analytical_sub, '-')
            plt.xlabel('x')
            plt.ylabel('rho')

            plt.figure(14)
            plt.plot(x, T, '.')
            if sample_exact == 2:
                if t > t_steady:
                    plt.plot(x, T_analytical_quasi, '-')
            if sample_exact == 3:
                if t > t_steady:
                    plt.plot(x, T_analytical_sub, '-')
            plt.xlabel('x')
            plt.ylabel('T')

            plt.figure(15)
            plt.plot(x, Mach, '.')
            if sample_exact == 2:
                if t > t_steady:
                    plt.plot(x, mach_analytical_quasi, '-')
            if sample_exact == 3:
                if t > t_steady:
                    plt.plot(x, mach_analytical_sub, '-')
            plt.xlabel('x')
            plt.ylabel('Mach')

            plt.figure(16)
            plt.plot(x, E_specific, '-')
            plt.xlabel('x')
            plt.ylabel('Energy')

            plt.figure(21)
            plt.plot(x, Q_exhaust[:,0], '-r')
            plt.xlabel('x')
            plt.ylabel('Q1')

            plt.figure(22)
            plt.plot(x, Q_exhaust[:, 1], '-g')
            plt.xlabel('x')
            plt.ylabel('Q2')

            plt.figure(23)
            plt.plot(x, Q_exhaust[:, 2], '-b')
            plt.xlabel('x')
            plt.ylabel('Q3')

            plt.figure(24)
            plt.plot(x, rho*u*Area_nozzle_centred, '-')
            plt.xlabel('x')
            plt.ylabel('Mass flux')

        # ----------------------------------------------------

    #plt.show()


    # ------------------------------------------
    # Close file
    f_io.close()
    # ------------------------------------------
    return [heat_flux_return, T, T_aw, h_coeff_exhaust, rho_return, u_return, p_return, E_return, L, x_min, x_max]
