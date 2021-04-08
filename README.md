# thesis_LRE_nozzle_simulation_program

USER INSTRUCTIONS

The user should download all .py filed into one folder. 

coupled_model provides all the user inputs, which are detailed below. The user should set inputs to their desired simulation conditions, \
and run the program through the coupled_model script

USER INPUTS

    nx = 51
    
    EXHAUST PARAMETERS
    
    chamber_pressure = 1000000            # combustion chamber pressure (stagnation pressure)
    chamber_temperature = 2600            # combustion chamber temperature (stagnation temperature)
    outlet_pressure = 101325              # exhaust outlet pressure. Doesn't matter for exhaust flow boundary condition is 'zero grad' and transmissive
    nozzle_geometry = 0                   # 0 for Rao nozzle, 1 for user-defined radius function, 2 for constant area
    exhaust_boundaries = 0                # 0 for subsonic inlet and supersonic outlet, 1 for subsonic inlet subsonic outlet
    exhaust_sample_exact = 2              # for computation of analytical solution. 1 for riemann, 2 for steady exhaust flow, 3 for fully subsonic cooling channel flow
    exhaust_heat = 0                      # 0 for exhaust flow, 1 for coolant flow
    coolant_wall_temp_exhaust = 1         # dummy variable. Advised to keep this as 1
    t_max_exhaust = 0.1                   # simulation time for exhaust flow
    t_steady_exhaust = 0.06               # look at residuals and define this at a time where resuduals reach an appropriately \
                                          low value. 0.06 is appropriate for these conditions
    T_w_hg = 1200                         #hot gas wall temp initial (STEP 1) guess
    
    DOMAIN GEOMETRY # user advised keep these as they are
    
    L = 1
    x_min = 0
    x_max = 1

    INITIAL EXHAUST FLOW CONDITIONS # dummy vars, advised to keep these as they are
    
    init_exhaust = 0
    init_rho_ex = 0
    init_u_ex = 0
    init_p_ex = 0

    file_no_ex = 0 # initial file number counter

    WALL PARAMETERS
    
    conductivity_k = 20                   # inconel, W/mK
    wall_thickness = 0.0015               # 1.5 mm wall
    no_channels = 100                     # number of cooling channels

    COOLANT PARAMETERS
    
    coolant_inlet_pressure = 101325 * 1.0075      # coolant inlet pressure (stagnation pressure)
    coolant_inlet_temperature = 273               # coolant inlet temperature (stagnation pressure)
    coolant_outlet_pressure = 101325              # coolant outlet pressure
    channel_geometry = 2                          # 0 for Rao nozzle, 1 for user-defined radius function, 2 for constant area
    coolant_boundaries = 1                        # 0 for subsonic inlet and supersonic outlet, 1 for subsonic inlet subsonic outlet
    coolant_sample_exact = 3                      # for computation of analytical solution. 1 for riemann, 2 for steady exhaust flow, 3 for fully subsonic cooling channel flow
    coolant_heat = 1                              # 0 for exhaust flow, 1 for coolant flow
    T_w_hg_coolant = 1                            # dummy variable. Advised to keep this as 1
    t_max_coolant = 0.1                           # simulation time for coolant flow
    t_steady_coolant = 0.08                       look at residuals and define this at a time where resuduals reach an appropriately \
                                                  low value. 0.08 is appropriate for these conditions


    INITIAL COOLANT FLOW CONDITIONS # assign coolant initial conditions at each new iteration as the solution of the last iteration 
     
    init_coolant = 1
    init_rho_cool = np.ones(nx)
    init_u_cool = np.ones(nx)
    init_p_cool = np.ones(nx)

    file_no_cool = 1

    for i in range(0, nx):                      # initial conditions for first coolant iteration 
        init_rho_cool[i] = 1.0
        init_u_cool[i] = 2.0
        init_p_cool[i] = 101325


    flux_tol = 15000                            # flux tolerance. 15,000 is advised to give trade-off between run-time and convergence
    T_w_hg_inc = 20                             #incremental increase of T_w_hg. Multiplied by convergence coefficient in script

    THERMODYNAMIC PROPERTIES

    visc_exhaust = 0.0001                       # exhaust viscosity
    cp_exhaust = 2519.0                         # exhaust specific heat
    k_exhaust = 0.3231                          # exhaust thermal conductivity

    visc_coolant = 0.00116                      # coolant viscosity
    cp_coolant = 2700                           # coolant specific heat
    k_coolant = 142                             # coolant thermal conductivity
