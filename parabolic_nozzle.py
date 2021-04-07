import math
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.interpolate import griddata
from scipy.interpolate import interp1d




def bell_coefficients(theta_n, theta_e, radius_throat, a_ratio, l_percent):
    N_x = (0.382 * radius_throat * math.cos(theta_n - math.pi/2))
    N_y = (0.382 * radius_throat * math.sin(theta_n - math.pi/2) + (0.382 * radius_throat) + radius_throat)

    E_x = l_percent * ((math.sqrt(a_ratio) - 1) * radius_throat) / math.tan(math.radians(15))
    E_y = math.sqrt(a_ratio) * radius_throat

    grad_1 = math.tan(theta_n)
    grad_2 = math.tan(theta_e)

    int_1 = N_y - (grad_1 * N_x)
    int_2 = E_y - (grad_2 * E_x)

    Q_x = (int_2 - int_1) / (grad_1 - grad_2)
    Q_y = ((grad_1 * int_2) - (grad_2 * int_1)) / (grad_1 - grad_2)

    return[N_x, N_y, E_x, E_y, Q_x, Q_y]



def get_contour(theta_n, theta_e, radius_throat, a_ratio, l_percent):

    #converging throat section
    converging_entrance_angle = -135
    converging_entrance_rad = math.radians(converging_entrance_angle)

    grid_points = 100

    conv_start = converging_entrance_rad
    conv_end = - math.pi / 2

    throat_contour_list = np.linspace(conv_start, conv_end, grid_points)
    x_point_conv = []
    y_point_conv = []

    for i in throat_contour_list:
        x_point_conv.append(1.5 * radius_throat * math.cos(i))
        y_point_conv.append(1.5 * radius_throat * math.sin(i) + (1.5 * radius_throat) + radius_throat)

    #diverging throat section
    diverg_start = - math.pi / 2
    diverg_end = theta_n - (math.pi / 2)
    throat_contour_list = np.linspace(diverg_start, diverg_end, grid_points)
    x_point_diverg = []
    y_point_diverg = []

    for i in throat_contour_list:
        x_point_diverg.append(0.382 * radius_throat * math.cos(i))
        y_point_diverg.append(0.382 * radius_throat * math.sin(i) + (0.382 * radius_throat) + radius_throat)

    #print('throat contour = ', throat_contour_list)
    #print('x throat contour conv = ', x_point_conv)
    #print('y throat contour conv = ', y_point_conv)

    #bell
    bell_start = 0
    bell_end = 1
    bell_contour_list = np.linspace(bell_start, bell_end, grid_points)
    x_point_bell = []
    y_point_bell = []

    [N_x, N_y, E_x, E_y, Q_x, Q_y] = bell_coefficients(theta_n, theta_e, radius_throat, a_ratio, l_percent)

    for t in bell_contour_list:
        x_point_bell.append(((1 - t) ** 2 * N_x) + (2 * (1 - t) * t * Q_x) + (t ** 2 * E_x))
        y_point_bell.append(((1 - t) ** 2 * N_y) + (2 * (1 - t) * t * Q_y) + (t ** 2 * E_y))


    return[x_point_conv, y_point_conv, x_point_diverg, y_point_diverg, x_point_bell, y_point_bell]



def get_angles(a_ratio, l_percent):

    if l_percent == 0.8:
        theta_n_80 = [19.5, 23.0, 23.9, 24.7, 25.3, 25.9, 26.2, 32.6]
        theta_e_80 = [16.0, 12.9, 12.0, 11.4, 11.1, 10.6, 10.3, 7.0]

        if a_ratio == 3:
            theta_n = theta_n_80[0]
            theta_e = theta_e_80[0]
        elif a_ratio == 5:
            theta_n = theta_n_80[1]
            theta_e = theta_e_80[1]
        elif a_ratio == 6:
            theta_n = theta_n_80[2]
            theta_e = theta_e_80[2]
        elif a_ratio == 7:
            theta_n = theta_n_80[3]
            theta_e = theta_e_80[3]
        elif a_ratio == 8:
            theta_n = theta_n_80[4]
            theta_e = theta_e_80[4]
        elif a_ratio == 9:
            theta_n = theta_n_80[5]
            theta_e = theta_e_80[5]
        elif a_ratio == 10:
            theta_n = theta_n_80[6]
            theta_e = theta_e_80[6]
        elif a_ratio == 70:
            theta_n = theta_n_80[7]
            theta_e = theta_e_80[7]


    return[math.radians(theta_n), math.radians(theta_e)]



def get_length(area_ratio, radius_throat, l_percent):

    length_func = ((math.sqrt(area_ratio) - 1) * radius_throat) / math.tan(math.radians(15))

    if l_percent == 0.8:
        l_nozzle = 0.8 * length_func

    return l_nozzle



def get_output(x_point_conv, y_point_conv, x_point_diverg, y_point_diverg, x_point_bell, y_point_bell):

    x_array_out = x_point_conv + x_point_diverg + x_point_bell

    y_array_out = y_point_conv + y_point_diverg + y_point_bell

    return[x_array_out, y_array_out]


def rao_nozzle(radius_throat, a_ratio, length_percent):

    nozzle_length = get_length(a_ratio, radius_throat, length_percent)

    [theta_n, theta_e] = get_angles(a_ratio, length_percent)

    [x_point_conv, y_point_conv, x_point_diverg, y_point_diverg, x_point_bell, y_point_bell] = \
        get_contour(theta_n, theta_e, radius_throat, a_ratio, length_percent)

    [x_array_out, y_array_out] = get_output(x_point_conv, y_point_conv, x_point_diverg, y_point_diverg, x_point_bell, y_point_bell)

    return[nozzle_length, x_array_out, y_array_out]



if __name__ == "__main__":

    radius_throat = 0.5
    a_ratio = 5
    length_percent = 0.8
    nx = 51
    nr = nx

    [nozzle_length, x_array_out, y_array_out] = rao_nozzle(radius_throat, a_ratio, length_percent)

    print('x out = ', x_array_out)
    print('y out = ', y_array_out)
    print('nozzle len = ', nozzle_length)

    L = nozzle_length

    x = np.linspace(0, L, nx)

    y_size = np.size(y_array_out)

    x_size = np.size(x_array_out)

    domain_size = np.size(x)

    print('x nozzle coord size = ', x_size)
    print('y nozzle coord size = ', y_size)
    print('main script domain size = ', domain_size)


    plt.figure(15)
    plt.plot(x_array_out, y_array_out, '-b')

    plt.xlabel('x')
    plt.ylabel('y')

    plt.figure(16)
    plt.plot(x_array_out, y_array_out)
    plt.xlabel('x')
    plt.ylabel('y')


plt.show()
