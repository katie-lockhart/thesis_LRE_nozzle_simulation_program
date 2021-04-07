import quasi_1d_model
import numpy as np
import os
import matplotlib.pyplot as plt


def solve_exhaust(chamber_pressure, chamber_temperature, outlet_pressure, nozzle_geometry,
                                        exhaust_boundaries, exhaust_sample_exact, exhaust_heat, coolant_wall_temp_exhaust,
                                T_w_hg, visc_exhaust, cp_exhaust, k_exhaust, t_max_exhaust, t_steady_exhaust, init,
                                init_rho, init_u, init_p, file_no_ex, L, x_min, x_max):

    [heat_flux_exhaust, T_exhaust_free_stream, T_adiabatic_wall, h_exhaust, rho_ex_return, u_ex_return, p_ex_return, E_ex_return,
                                        x_domain, x_min, x_max] = quasi_1d_model.solve_flow(chamber_pressure,
                                            chamber_temperature, outlet_pressure, nozzle_geometry, exhaust_boundaries,
                                                    exhaust_sample_exact, exhaust_heat, coolant_wall_temp_exhaust, T_w_hg,
                                                            visc_exhaust, cp_exhaust, k_exhaust, t_max_exhaust, t_steady_exhaust, init,
                                                                     init_rho, init_u, init_p, file_no_ex, L, x_min, x_max)

    return [heat_flux_exhaust, T_exhaust_free_stream, T_adiabatic_wall, h_exhaust, rho_ex_return, u_ex_return, p_ex_return, E_ex_return, x_domain, x_min, x_max]
    # returns heat flux out across ALL channels

def wall_conduction(heat_flux_exhaust, conductivity_k, T_w_hg, wall_thickness, no_channels):

    coolant_wall_temp = T_w_hg - (heat_flux_exhaust * wall_thickness / conductivity_k)
    flux_out_per_channel = (heat_flux_exhaust / no_channels)

    return [coolant_wall_temp, flux_out_per_channel]


def solve_coolant(coolant_inlet_pressure, coolant_inlet_temperature, coolant_outlet_pressure, channel_geometry,
                  coolant_boundaries, coolant_sample_exact, coolant_heat, coolant_wall_temp, T_w_hg,
                  visc_coolant, cp_coolant, k_coolant, t_max_coolant, t_steady_coolant, init, init_rho, init_u, init_p,
                  file_no_cool, x_domain, x_min, x_max):

    [heat_flux_coolant, T_fs_c, T_aw_c, h_coeff_coolant, rho_cool_return, u_cool_return, p_cool_return, E_cool_return, x_domain,
                        x_min, x_max] = quasi_1d_model.solve_flow(coolant_inlet_pressure, coolant_inlet_temperature,
                                coolant_outlet_pressure, channel_geometry, coolant_boundaries, coolant_sample_exact,
                                            coolant_heat, coolant_wall_temp, T_w_hg, visc_coolant, cp_coolant, k_coolant, t_max_coolant,
                                                            t_steady_coolant, init, init_rho, init_u, init_p, file_no_cool, x_domain, x_min, x_max)

    return [heat_flux_coolant, T_fs_c, T_aw_c, h_coeff_coolant, rho_cool_return, u_cool_return, p_cool_return, E_cool_return, x_domain, x_min, x_max]
    # returns heat flux in to coolant PER CHANNEL

def normalize(heat_flux_exhaust, new_lower_bound, new_upper_bound):

    min_flux = np.min(heat_flux_exhaust)
    max_flux = np.max(heat_flux_exhaust)
    range = max_flux - min_flux
    new_range = new_upper_bound - new_lower_bound

    return [((i - min_flux) / range) * new_range + new_lower_bound for i in heat_flux_exhaust]


def solve_coupled_flow(chamber_pressure, chamber_temperature, outlet_pressure, nozzle_geometry,
                            exhaust_boundaries, exhaust_sample_exact, exhaust_heat, coolant_wall_temp_exhaust,
                            T_w_hg, conductivity_k, wall_thickness, coolant_inlet_pressure, coolant_inlet_temperature,
                            coolant_outlet_pressure, channel_geometry, coolant_boundaries, coolant_sample_exact, coolant_heat,
                            T_w_hg_coolant, tol, increment, visc_exhaust, cp_exhaust, k_exhaust, visc_coolant, cp_coolant, k_coolant,
                            t_max_exhaust, t_steady_exhaust, t_max_coolant, t_steady_coolant, init_exhaust, init_rho_ex, init_u_ex,
                            init_p_ex, init_coolant, init_rho_cool, init_u_cool, init_p_cool, file_no_ex, file_no_cool, L, x_min, x_max, no_channels):

    # ----------------------------------------------------
    # STEP 2
    # get frozen exhaust state and initial values for coolant

    [heat_flux_exhaust, T_exhaust_free_stream, T_adiabatic_wall, h_exhaust, rho_ex_return, u_ex_return, p_ex_return, E_ex_return, x_domain,
                                        x_min, x_max] = solve_exhaust(chamber_pressure, chamber_temperature,
                                        outlet_pressure, nozzle_geometry, exhaust_boundaries, exhaust_sample_exact, exhaust_heat,
                                        coolant_wall_temp_exhaust, T_w_hg, visc_exhaust, cp_exhaust, k_exhaust, t_max_exhaust, t_steady_exhaust, init_exhaust,
                                        init_rho_ex, init_u_ex, init_p_ex, file_no_ex, L, x_min, x_max)


    #STEP 3
    # get exhaust side wall temp distribution and first guess

    T_w_hg_old = normalize(heat_flux_exhaust, 500, 900)

    size_flux = heat_flux_exhaust.size
    x = np.linspace(x_min, x_max, 51)

    # STEP 4
    # wall conduction, get coolant wall temperature distribution
    [coolant_wall_temp, flux_out_per_channel] = wall_conduction(heat_flux_exhaust, conductivity_k, T_w_hg_old, wall_thickness, no_channels)

    plt.figure(31)
    plt.plot(x, coolant_wall_temp, '.')
    plt.xlabel('x')
    plt.ylabel('coolant wall temp')


    plt.figure(33)
    plt.plot(x, (T_w_hg_old), '.')
    plt.xlabel('x')
    plt.ylabel('Twh')

  #  plt.show()

    # STEP 5
    # solve coolant flow

    [heat_flux_coolant, T_fs_c, T_aw_c, h_coeff_coolant, rho_cool_return, u_cool_return, p_cool_return, E_cool_return, x_domain,
                                    x_min, x_max] = solve_coolant(coolant_inlet_pressure, coolant_inlet_temperature,
                                    coolant_outlet_pressure, channel_geometry, coolant_boundaries, coolant_sample_exact, 2,
                                    coolant_wall_temp, T_w_hg_coolant, visc_coolant, cp_coolant, k_coolant, t_max_coolant, t_steady_coolant,
                                    init_coolant, init_rho_cool, init_u_cool, init_p_cool, file_no_cool, x_domain, x_min, x_max)


    file_no_cool_old = file_no_cool
    itr_number_old = 0
    rho_cool_return_old = rho_cool_return
    u_cool_return_old = u_cool_return
    p_cool_return_old = p_cool_return

    heat_flux_in_per_channel = heat_flux_coolant / no_channels

    heat_flux_difference = (flux_out_per_channel - heat_flux_in_per_channel) #per channel

    domain = np.linspace(x_min, x_max, 51)

    file_name = folder_out_2 + '/coupled_output.iter' + str(itr_number_old) + '.dat'  # file
    f_itr = open(file_name, 'w')

    size_arr = heat_flux_exhaust.size
    for i in range(0, size_arr):
        s_tot = str(domain[i]) + '         ' + str(flux_out_per_channel[i]) + '         ' + str(heat_flux_in_per_channel[i]) + '         ' + str(heat_flux_difference[i]) + '         ' \
                + str(T_fs_c[i]) + '         ' + str(coolant_wall_temp[i]) + '         ' + str(T_w_hg) + '         ' + str(h_exhaust[i]) + '         ' + str(h_coeff_coolant[i])
        f_itr.write(s_tot + '\n')

    f_itr.close()

    x = np.linspace(0, size_arr, 51)

    # ----------------------------------------------------

    # ----------------------------------------------------
    # iterate on the coolant wall temp

    #STEP 6
    while np.max(np.absolute(heat_flux_difference)) > tol:

        itr_number_new = itr_number_old + 1
        file_no_cool_new = file_no_cool_old + 1

        init_rho_iter_new = rho_cool_return_old
        init_p_iter_new = p_cool_return_old
        init_u_iter_new = u_cool_return_old

        size_temp = size_flux
        T_w_hg_new = np.ones(size_temp)

        #STEP 7
        # get new exhaust-side wall temperature guess
        for i in range(size_temp):

            T_w_hg_new[i] = T_w_hg_old[i] + increment * (heat_flux_difference[i]*6e-6)

        print('new wall temp', T_w_hg_new)

        heat_flux_exhaust_new = h_exhaust * (T_adiabatic_wall - T_w_hg_new)     #new flux across ALL CHANNELS

        [coolant_wall_temp, flux_out_per_channel_new] = wall_conduction(heat_flux_exhaust_new, conductivity_k, T_w_hg_new, wall_thickness, no_channels)

        [heat_flux_coolant, T_fs_c, T_aw_c, h_coeff_coolant, rho_cool_return, u_cool_return, p_cool_return, E_cool_return, x_domain,
                                                                x_min, x_max] = solve_coolant(coolant_inlet_pressure, coolant_inlet_temperature,
                                                                       coolant_outlet_pressure, channel_geometry,
                                                                       coolant_boundaries, coolant_sample_exact,
                                                                       coolant_heat, coolant_wall_temp, T_w_hg_coolant, visc_coolant, cp_coolant, k_coolant,
                                                                        t_max_coolant, t_steady_coolant, 1,
                                                                        init_rho_iter_new, init_u_iter_new, init_p_iter_new, file_no_cool_new, x_domain, x_min, x_max)     #new coolant solution

        itr_number_old = itr_number_new
        file_no_cool_old = file_no_cool_new
        T_w_hg_old = T_w_hg_new
        heat_flux_in_per_channel_new = heat_flux_coolant
        heat_flux_difference = (flux_out_per_channel_new - heat_flux_in_per_channel_new) #out of exhaust per channel, in to each cooling channel

        rho_cool_return_old = rho_cool_return
        u_cool_return_old = u_cool_return
        p_cool_return_old = p_cool_return


        file_name = folder_out_2 + '/coupled_output_iter' + str(itr_number_new) + '.dat'  # file to store parameters at every iteration
        f_itr = open(file_name, 'w')

        size_arr = heat_flux_exhaust.size
        for i in range(0, size_arr):
            s_tot = str(domain[i]) + '         ' + str(flux_out_per_channel_new[i]) + '         ' + str(heat_flux_in_per_channel_new[i]) + '         ' + str(heat_flux_difference[i]) \
                    + '         ' + str(T_fs_c[i]) + '         ' + str(coolant_wall_temp[i]) + '         ' + str(T_w_hg_new[i]) + '         ' + str(h_exhaust[i]) + '         ' + str(h_coeff_coolant[i])
            f_itr.write(s_tot + '\n')

        f_itr.close()

    file_name = 'coupled_output_final_inner_iter' + '.dat'  # file to give final exhaust wall temp to outer iteration
    f_final = open(file_name, 'w')

    for i in range(0, size_arr):
        s_tot = str(T_w_hg_new[i])
        f_final.write(s_tot + '\n')

    f_final.close()

    print('fluxes converged within tol')
    return [flux_out_per_channel, heat_flux_coolant, heat_flux_difference, coolant_wall_temp, T_w_hg_new, file_no_cool_old]


if __name__ == "__main__":

    #-------------------------------------------------

    nx = 51

    chamber_pressure = 1000000
    chamber_temperature = 2600
    outlet_pressure = 101325    #doesn't matter for exhaust as zero grad outlet boundary
    nozzle_geometry = 0         #0 for Rao nozzle, 1 for radius defined, 2 for constant area
    exhaust_boundaries = 0      #0 for subsonic inlet and supersonic outlet, 1 for subsonic inlet subsonic outlet
    exhaust_sample_exact = 2    #1 for riemann, 2 for choked exhaust flow, 3 for fully subsonic
    exhaust_heat = 0            #0 for exhaust flow, 1 for coolant flow
    coolant_wall_temp_exhaust = 1
    t_max_exhaust = 0.05
    t_steady_exhaust = 0.018

    L = 1
    x_min = 0
    x_max = 1

    init_exhaust = 0
    init_rho_ex = 0
    init_u_ex = 0
    init_p_ex = 0

    file_no_ex = 0

    T_w_hg = 1200               #hot gas wall temp initial (STEP 1) guess
    conductivity_k = 20         #inconel, W/mK
    wall_thickness = 0.0015     #1.5 mm wall
    no_channels = 100

    coolant_inlet_pressure = 101325 * 1.0075
    coolant_inlet_temperature = 273
    coolant_outlet_pressure = 101325
    channel_geometry = 2
    coolant_boundaries = 1
    coolant_sample_exact = 3
    coolant_heat = 1
    T_w_hg_coolant = 1
    t_max_coolant = 0.1
    t_steady_coolant = 0.08


    init_coolant = 1
    init_rho_cool = np.ones(nx)
    init_u_cool = np.ones(nx)
    init_p_cool = np.ones(nx)

    file_no_cool = 1

    for i in range(0, nx):
        init_rho_cool[i] = 1.0
        init_u_cool[i] = 2.0
        init_p_cool[i] = 101325


    flux_tol = 15000
    T_w_hg_inc = 20             #incremental increase of T_w_hg


    visc_exhaust = 0.0001
    cp_exhaust = 2519.0
    k_exhaust = 0.3231

    visc_coolant = 0.00116
    cp_coolant = 2700
    k_coolant = 142


    folder_out_2 = 'output_coupled_model_iterations'

    if not os.path.exists(folder_out_2):
        os.makedirs(folder_out_2)


    [heat_flux_exhaust, heat_flux_coolant, heat_flux_difference, coolant_wall_temp, exhaust_wall_temp, last_file_no] = solve_coupled_flow(chamber_pressure,
                        chamber_temperature, outlet_pressure, nozzle_geometry, exhaust_boundaries,
                        exhaust_sample_exact, exhaust_heat, coolant_wall_temp_exhaust, T_w_hg, conductivity_k, wall_thickness,
                        coolant_inlet_pressure, coolant_inlet_temperature, coolant_outlet_pressure, channel_geometry,
                        coolant_boundaries, coolant_sample_exact,coolant_heat, T_w_hg_coolant, flux_tol, T_w_hg_inc,
                        visc_exhaust, cp_exhaust, k_exhaust, visc_coolant, cp_coolant, k_coolant, t_max_exhaust, t_steady_exhaust, t_max_coolant, t_steady_coolant,
                        init_exhaust, init_rho_ex, init_u_ex, init_p_ex, init_coolant, init_rho_cool, init_u_cool, init_p_cool,
                        file_no_ex, file_no_cool, L, x_min, x_max, no_channels)


    size_arr = heat_flux_exhaust.size
    exhaust_wall_temp_final = np.ones(size_arr)
    exhaust_wall_temp_final = np.loadtxt('./coupled_output_final_inner_iter.dat', usecols=(0))
    print('walltemp = ', exhaust_wall_temp_final)

    last_file_no_new = last_file_no + 1

    # STEP 8
    [heat_flux_exhaust_final, T_exhaust_free_stream_final, T_adiabatic_wall_final, h_exhaust_final, rho_ex_return, u_ex_return,
                                        p_ex_return, E_ex_return, x_domain, x_min, x_max] = solve_exhaust(chamber_pressure,
                                        chamber_temperature, outlet_pressure, nozzle_geometry, exhaust_boundaries,
                                        exhaust_sample_exact, exhaust_heat, coolant_wall_temp_exhaust, exhaust_wall_temp_final,
                                        visc_exhaust, cp_exhaust, k_exhaust, t_max_exhaust, t_steady_exhaust, init_exhaust,
                                        init_rho_ex, init_u_ex, init_p_ex, last_file_no_new, L, x_min, x_max) #all channels

    flux_out_per_channel_final = heat_flux_exhaust_final / no_channels #per channel
    difference_final = flux_out_per_channel_final - heat_flux_coolant
    # ----------------------------------------------------
    # Outer iteration output file
    domain = np.linspace(x_min, x_max, 51)
    file_name = folder_out_2 + '/exhaust_final_output.dat' #file
    f_ex = open(file_name,'w')

    size_arr = heat_flux_exhaust.size
    for i in range(0, size_arr):
        s_tot = str(domain[i]) + '  ' + str(flux_out_per_channel_final[i]) + '  ' + str(heat_flux_coolant[i]) + '  ' + str(difference_final[i]) + '  ' + \
                str(exhaust_wall_temp_final[i]) + '  ' + str(coolant_wall_temp[i]) + '  ' + str(h_exhaust_final[i])
        f_ex.write(s_tot + '\n')

    f_ex.close()
    #-----------------------------------------------------