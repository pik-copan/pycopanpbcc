# -*- coding: utf-8 -*-
# Author: Vera Heck <heck@pik-potsdam.de>
# Script generates data and plots Fig. 6 of Heck et al. 2016 (ESD)

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib
import sys
import brewer2mpl
from scipy import integrate

array = np.array
import main_model as main


def plot_state_space_dots(fig,evol, tmax, n_tm, b_a, b_m, b_l, ii):
    AB_a = np.linspace(0.00, 0.4, n_tm)
    AB_t = np.linspace(0.4, 0.80, n_tm)
    
    ts = np.linspace(0, tmax, tmax)  # time grid 

    if ii==1: 
        ax = fig.add_subplot(131, projection='3d')
        main.thresh_geo = 0.16
        main.alpha_max = 0.000
        main.s_geo = 200.
        main.c_max = 0.25
        ax.text(0.26, 0.2,0.48, 'a) no tCDR',  fontsize=8, va='top')
        
    if ii==2:
        ax = fig.add_subplot(132, projection='3d')
        main.thresh_geo = 0.16
        main.alpha_max = 0.004
        main.s_geo = 200.
        main.c_max = 0.2
        print(b_a, b_l,b_m, main.thresh_geo, main.alpha_max, main.c_max, main.s_geo) 
        ax.text(0.26, 0.2,0.48, 'b) medium tCDR rate',  fontsize=8, va='top')
        
    if ii==3: 
        ax = fig.add_subplot(133, projection='3d')
        main.thresh_geo = 0.16
        main.alpha_max = 0.04
        main.s_geo = 200.
        main.c_max = 0.2
        ax.text(0.24, 0.2,0.48, 'c) high tCDR rate',  fontsize=8, va='top')
   
    ax.set_xlim(0.4, 0.815) 
    ax.set_ylim(0.45, 0.05) 
    ax.set_zlim(0.0, 0.4) 

    ax.set_xticklabels([0.4,"",0.5,"",0.6,"",0.7,"",0.8], fontsize=6,va='center')
    ax.set_yticklabels( ["",0.1,"",0.2,"",0.3,"",0.4,""], rotation = -1,fontsize=6,va='center', ha='left')
    ax.set_zticklabels(["","",0.1,"",0.2,"",0.3, "", 0.4], fontsize=6,va='center')
    ax.view_init(*(20, -68))
    
    ax.tick_params(axis='x', pad=-3)
    ax.tick_params(axis='y', pad=-5)
    ax.tick_params(axis = 'z', pad=-2)

    ax.text(.63, 0.465, 0, 'land carbon',(1,-0.02,0) , fontsize=6)
    ax.text(0.698, 0.19, 0, 'ocean carbon', (0.01,1,0), fontsize=6)
    ax.text(0.776, .0, .31, 'atmosphere carbon', 'z', fontsize=6)   

    set2 = brewer2mpl.get_map('Set3', 'Qualitative', 12).mpl_colors
    c_r0 = 5e-4
    c_g0 = 0.00
    s = alm = am = al = lm = l = a = m = c_g = c_r = 0.0
    i=0
    for c_a0 in AB_a:
        for c_t0 in AB_t:
            c_m0 = main.C_0 - c_a0 - c_t0
            if c_m0 >= 0.04999999999 and c_m0 <= 0.45000001:  # solve the DEs
                i=i+1
                y0 = [c_r0, c_m0, c_t0, c_g0]  # initial condition vector
                soln = integrate.odeint(evol, y0, ts)
                C_R = soln[:, 0]
                C_M = soln[:, 1]
                C_T = soln[:, 2]
                C_G = soln[:, 3]
                C = main.C_0 + C_R - C_G
                C_A = C - C_M - C_T
                if np.all(np.round(C_T, 5) > b_l) & np.all(np.round(C_A, 5) < b_a) & np.all(
                                np.round(C_M, 5) < b_m):  # SOS
                    if s == 0:
                        ax.scatter(c_t0, c_m0, c_a0, alpha=1, c=set2[6], edgecolors='None', lw = 0, label="SOS", s=0.5)  
                    else:
                        ax.scatter(c_t0, c_m0, c_a0, alpha=1, c=set2[6], edgecolors='None', lw = 0, s=0.5)
                    s += 1
                elif np.any(np.round(C_T, 5) <= b_l) & np.any(np.round(C_A, 5) >= b_a) & np.any(
                                np.round(C_M, 5) >= b_m):  #atm, ocean and land transgressed
                    if alm == 0:
                        ax.scatter(c_t0, c_m0, c_a0, alpha=1, c=set2[3], edgecolors='None', lw = 0, label="outside SOS", s=0.5)
                    else:
                        ax.scatter(c_t0, c_m0, c_a0, alpha=1, c=set2[3],edgecolors='None',  lw = 0,s=0.5)
                    alm += 1
                elif np.all(np.round(C_T, 5) > b_l) & np.any(np.round(C_A, 5) >= b_a) & np.any(
                                np.round(C_M, 5) >= b_m):  #atm and ocean transgressed
                    if am == 0:
                        ax.scatter(c_t0, c_m0, c_a0, alpha=1, c=set2[11], edgecolors='None', lw = 0, label="land SOS", s=0.5)
                    else:
                        ax.scatter(c_t0, c_m0, c_a0, alpha=1, c=set2[11], edgecolors='None',  lw = 0,s=0.5)
                    am += 1
                elif np.any(np.round(C_T, 5) <= b_l) & np.any(np.round(C_A, 5) >= b_a) & np.all(
                                np.round(C_M, 5) < b_m):  #atm & land transgressed
                    if al == 0:
                        ax.scatter(c_t0, c_m0, c_a0, alpha=1, c=set2[4], edgecolors='None', lw = 0, label="ocean SOS", s=0.5)
                    else:
                        ax.scatter(c_t0, c_m0, c_a0, alpha=1, c=set2[4], edgecolors='None', lw = 0, s=0.5)
                    al += 1
                elif np.any(np.round(C_T, 5) <= b_l) & np.all(np.round(C_A, 5) < b_a) & np.any(
                                np.round(C_M, 5) >= b_m):  #ocean and land transgressed
                    if lm == 0:
                        ax.scatter(c_t0, c_m0, c_a0, alpha=1, c='LightGrey', edgecolors='None', lw = 0, label="atm SOS", s=0.5)
                    else:
                        ax.scatter(c_t0, c_m0, c_a0, alpha=1, c='LightGrey', edgecolors='None', lw = 0, s=0.5)
                    lm += 1
                elif np.any(np.round(C_A, 5) >= b_a) & np.all(np.round(C_T, 5) > b_l) & np.all(
                                np.round(C_M, 5) < b_m):  # atm transgressed
                    if a == 0:
                        ax.scatter(c_t0, c_m0, c_a0, alpha=1, c=set2[2], edgecolors='None', lw = 0, label="land & ocean SOS", s=0.5)
                    else:
                        ax.scatter(c_t0, c_m0, c_a0, alpha=1, c=set2[2], edgecolors='None', lw = 0, s=0.5)
                    a += 1
                elif np.any(np.round(C_T, 5) <= b_l) & np.all(np.round(C_A, 5) < b_a) & np.all(
                                np.round(C_M, 5) < b_m):  # land transgressed
                    if l == 0:
                        ax.scatter(c_t0, c_m0, c_a0, alpha=1, c=set2[9], edgecolors='None', lw = 0, label="atm & ocean SOS", s=0.5)
                    else:
                        ax.scatter(c_t0, c_m0, c_a0, alpha=1, c=set2[9], edgecolors='None', lw = 0, s=0.5)
                    l += 1
                elif np.any(np.round(C_M, 5) >= b_m) & np.all(np.round(C_A, 5) < b_a) & np.all(
                                np.round(C_T, 5) > b_l):  # oecan transgressed
                    if m == 0:
                        ax.scatter(c_t0, c_m0, c_a0, alpha=1, c=set2[10], edgecolors='None', lw = 0, label="atm & land SOS", s=0.5)
                    else:
                        ax.scatter(c_t0, c_m0, c_a0, alpha=1, c=set2[10], edgecolors='None', lw = 0, s=0.5)
                    m += 1
    # plot current and preind state (normalized to 3989 PgC , IPCC 2014)
    ax.scatter(0.629, 0.228, 0.15, c='Black', edgecolors='None',marker='*', s=30)  # preindustrial

    ax.scatter(0.619, 0.264, 0.208, c='red', edgecolors='None',marker='*', s=30)  # current
    
    return c_g, c_r,i

if __name__ == "__main__":
    main.alpha_max = 0.00
    main.s_geo = 200.

    tmax = 600 
    n_tm = 300
    
    # boundaries
    b_a = 0.21 
    b_m = 0.31
    b_l = 0.59#43
      
    fig = plt.figure(figsize=(18/2.54, 6.5/2.54))
    plt.subplots_adjust(left=0, bottom=0.0, right=0.999, top=None, wspace=-0.03, hspace=-0.1)
    for ii in [1,2,3]: 
        [c_g, c_r,i] = plot_state_space_dots(fig, main.deriv, tmax, n_tm, b_a, b_m, b_l,ii)
    set2 = brewer2mpl.get_map('Set3', 'Qualitative', 12).mpl_colors
    proxy1 = matplotlib.lines.Line2D([1], [0], linestyle="none",  markeredgewidth=0, c=set2[6], marker='o',  alpha=1) #sos
    proxy2 = matplotlib.lines.Line2D([0], [0], linestyle="none", markeredgewidth=0,c=set2[3], marker='o', alpha=1)#no sos
    proxy3 = matplotlib.lines.Line2D([0], [0], linestyle="none", markeredgewidth=0,c=set2[11], marker='o', alpha=1) #land
    proxy4 = matplotlib.lines.Line2D([0], [0], linestyle="none", markeredgewidth=0,c=set2[4], marker='o', alpha=1) #ocean
    proxy5 = matplotlib.lines.Line2D([0], [0], linestyle="none", markeredgewidth=0,c='LightGrey', marker='o', alpha=1) #atm
    proxy6 = matplotlib.lines.Line2D([0], [0], linestyle="none", markeredgewidth=0,c=set2[2], marker='o', alpha=1) #land-ocean 
    proxy7 = matplotlib.lines.Line2D([0], [0], linestyle="none", markeredgewidth=0,c=set2[9], marker='o', alpha=1) #atm-ocean
    proxy8 = matplotlib.lines.Line2D([0], [0], linestyle="none", markeredgewidth=0,c=set2[10], marker='o', alpha=1) #atm-land
    proxy9 = matplotlib.lines.Line2D([0], [0], linestyle="none", markeredgewidth=0,c='Black', marker='*', alpha=1)
    proxy10 = matplotlib.lines.Line2D([0], [0], linestyle="none", markeredgewidth=0,c='red', marker='*', alpha=1)
 

    fig.legend([proxy1, proxy2, 
               proxy6, proxy7, proxy8,
               proxy3, 
               proxy4, proxy5,
               proxy9, proxy10],
            ['MCSOS', 'not safe',
            'land-ocean MD', 'atmosphere-ocean MD', 'atmosphere-land MD',
             'land MD',
            'ocean MD', 'atmosphere MD',
            'preindustrial state', 'current state'],
            'upper center',numpoints=1, fontsize=6, ncol=5,frameon=False)
            
    print i    # number of iterations    

    saveToFile = "FIG6.png" 

    plt.savefig(saveToFile, format='png', dpi=1200)

