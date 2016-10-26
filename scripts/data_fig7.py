# -*- coding: utf-8 -*-
# Author: Vera Heck <heck@pik-potsdam.de>
# Script generates data for Fig. 7 of Heck et al. 2016 (ESD)

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy import integrate

array = np.array

# import model code
import main_model as main


def state_space_size(evol, tmax, n_tm, b_a, b_m, b_l):
    # loop over initial conditions
    AB_a = np.linspace(0.0, 0.4, n_tm)
    AB_t = np.linspace(0.4, 0.8, n_tm)
    ts = np.linspace(0, tmax, tmax)  # time grid
    c_r0 = 5e-4
    c_g0 = 0.00
    s = alm = am = al = lm = l = a = m = 0.0
    for c_a0 in AB_a:
        for c_t0 in AB_t:
            c_m0 = main.C_0 - c_a0 - c_t0
            if c_m0 >= 0.04999 and c_m0 <= 0.45011:  # solve the DEs 
                y0 = [c_r0, c_m0, c_t0, c_g0]  # initial condition vector
                soln = integrate.odeint(evol, y0, ts)
                C_R = soln[:, 0]
                C_M = soln[:, 1]
                C_T = soln[:, 2]
                C_G = soln[:, 3]
                C = main.C_0 + C_R - C_G
                C_A = C - C_M - C_T
                if np.all(np.round(C_T, 5) >= b_l) & np.all(np.round(C_A, 5) <= b_a) & np.all(
                                np.round(C_M, 5) <= b_m):  # SOS
                    s += 1
                elif np.any(np.round(C_T, 5) < b_l) & np.any(np.round(C_A, 5) > b_a) & np.any(
                                np.round(C_M, 5) > b_m):  # atm, ocean and land transgressed
                    alm += 1
                elif np.all(np.round(C_T, 5) > b_l) & np.any(np.round(C_A, 5) > b_a) & np.any(
                                np.round(C_M, 5) > b_m):  # atm and ocean transgressed
                    am += 1
                elif np.any(np.round(C_T, 5) < b_l) & np.any(np.round(C_A, 5) > b_a) & np.all(
                                np.round(C_M, 5) <= b_m):  # atm & land transgressed
                    al += 1
                elif np.any(np.round(C_T, 5) < b_l) & np.all(np.round(C_A, 5) <= b_a) & np.any(
                                np.round(C_M, 5) > b_m):  # ocean and land transgressed
                    lm += 1
                elif np.any(np.round(C_A, 5) > b_a) & np.all(np.round(C_T, 5) >= b_l) & np.all(
                                np.round(C_M, 5) <= b_m):  # atm transgressed
                    a += 1
                elif np.any(np.round(C_T, 5) < b_l) & np.all(np.round(C_A, 5) <= b_a) & np.all(
                                np.round(C_M, 5) <= b_m):  # land transgressed
                    l += 1
                elif np.any(np.round(C_M, 5) > b_m) & np.all(np.round(C_A, 5) <= b_a) & np.all(
                                np.round(C_T, 5) >= b_l):  # ocean transgressed
                    m += 1

    A = s + alm + am + al + lm + l + a + m  # size of possible space 
    l_safe = am + a + m + s
    m_safe = al + a + l + s
    a_safe = lm + l + m + s
    return np.array([s/A, a_safe/A, m_safe/A, l_safe/A, alm/A])

if __name__ == "__main__":
 
    var_thresh = True
    var_alpha_max = True

    b_a = 0.21
    b_m = 0.31
    b_l = 0.59
    nstep = 150  # number of variation points (points form line)
    tmax = 600  # differential equ ation time frame
    n_tm = 100  # number of startimg conditions for one (3d) 
    r_i = 0.03
    alpha = 0.0004
    s_geo = 200.
    c_max=0.4
    
    #data for Fig 7 a:
    alpha_max= 0.02
                  
    if var_thresh:
        main.alpha_max = alpha_max  
        main.alpha = alpha
        main.c_max = c_max
        main.r_i = r_i
        i = 0
        areas_thresh = np.empty((nstep, 5), dtype=object)
        
        for main.thresh_geo in np.linspace(0, 0.3, nstep):
            a = state_space_size(main.deriv, tmax, n_tm, b_a, b_m, b_l)
            areas_thresh[i,] = a
            i += 1
        s_file = '/save/var_thresh_alpha_max_' + str(alpha_max)+ '_c_max_'+str(c_max) + '.dat'
        np.save(s_file, areas_thresh)
                       
                                 
    #data for Fig 7 b:
    thresh_geo=0.2
    
    if var_alpha_max:
        main.thresh_geo = thresh_geo
        main.alpha = alpha
        main.c_max = c_max
        main.r_i = r_i

        i = 0
        areas_alpha_max = np.empty((nstep, 5), dtype=object)

        for main.alpha_max in np.linspace(0.000001, 0.03, nstep):  
            a = state_space_size(main.deriv, tmax, n_tm, b_a, b_m, b_l)
            areas_alpha_max[i,] = a
            i = i + 1
        s_file = '/save/var_alpha_max_thresh_' + str(thresh_geo)+ '_c_max_'+str(c_max) + '.dat'
        np.save(s_file, areas_alpha_max)
                            
