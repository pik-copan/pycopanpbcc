# -*- coding: utf-8 -*-
# Author: Vera Heck <heck@pik-potsdam.de>
# Script generates data for Fig. 5 of Heck et al. 2016 (ESD)

# Model main

# This file is 
import numpy as np
from math import exp
import matplotlib
import matplotlib.pyplot as plt
from scipy import integrate

plt.ion()

# Definde parameters

# fixed parameters
C_0 = 1.0  # initial carbon stock
b_T = 0.227 # intercept of T - c_a relationship
a_T = 1.06 # slope of T - c_a relationship
beta = 0.654 # carbon solubility in sea water factor 
alpha = 0.0004 # human terrestrial carbon offtake 
a_m = 0.0166  # atmosphere ocean diffusion coefficient
r_tc = 2.5  # ecosystem dependent conversion factor 

a_p = 0.48 # scaling factor for T  
a_r = 0.4 # scaling factor for T  
b_r = 0.5 # exponent for increase in r(T) for low T
b_p = b_r  # exponent for increase in p(T) for low T
c_p = 0.833  # exponent for decrease in p(T) for high T
c_r = 0.556 # exponent for decrease in r(T) for high T

# varied parameters
a_k = -0.6      # terrestrial carbon carrying capacity 
b_k = -13.0     # terrestrial carbon carrying capacity 
c_k = 0.75       # terrestrial carbon carrying capacity 

global s_geo, thresh_geo, alpha_max, r_i, c_max

# geoengineering parameters
s_geo = 200
thresh_geo = 0.18
alpha_max = 0.005

# industrialization parameters
r_i = 0.03
c_max = 0.51

# define equations
def C_a(C, c_m, c_t):  
    return C - c_t - c_m

def T(c_a):
    return a_T * c_a + b_T

def P(t, c_a):
    return  a_p * np.power(t, b_p) * np.exp(-c_p * t)
  
def R(t):
    return a_r * np.power(t, b_r) * np.exp(-c_r * t)

def K(c_a):
    return a_k*np.exp(b_k*c_a)+c_k
    
def NEP(t, c_a, c_t):
    return r_tc * c_t *(1 - (c_t / K(c_a))) * (P(t, c_a) - R(t))

def H(c_t):
    return alpha * c_t

def H_geo(c_t, c_a):
    return c_t * alpha_max / (1 + np.exp(-s_geo * (c_a - thresh_geo)))

def I(c_r):
    return r_i * c_r * (1 - c_r / c_max)


# solve the system dy/dt = f(y, t)
def deriv(y, t):
    c_r = y[0]
    c_m = y[1]
    c_t = y[2]
    c_g = y[3]
    C = C_0 + c_r - c_g
    c_a = C_a(C, c_m, c_t)
    t = T(c_a)
    f0 = I(c_r)
    f1 = a_m * (C_a(C, c_m, c_t) - beta * c_m)
    f2 = NEP(t, c_a, c_t) - H(c_t) - H_geo(c_t, c_a)
    f3 = H_geo(c_t, c_a)
    return np.array([f0, f1, f2, f3])

# initial conditions
c_m0 = 0.226  # initial ocean carbon (preind)
c_t0 = 0.627 # initial terrestrial carbon (preind)
c_r0 = 5e-4
c_g0 = 0.00
y0 = np.array([c_r0, c_m0, c_t0, c_g0])  # initial condition vector
ts = np.linspace(0, 600, 600)  # time grid: define integration timeframe

if __name__ == "__main__":
    alpha1=0.00
    alpha2= 0.0025
    alpha3=0.025

    saveToFile = "FIG5.pdf"   
    f, (ax1, ax2,ax3) = plt.subplots(1, 3, frameon=False,sharey=True, figsize=(13/2.54, 7.5/2.54))
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.05, hspace=None)
    
    alpha_max = alpha1
    #  solve the DEs
    soln = integrate.odeint(deriv, y0, ts)  # ODE solver
    C_R = soln[:, 0]
    C_M = soln[:, 1]
    C_T = soln[:, 2]
    C_G = soln[:, 3]
    C = C_0 + C_R - C_G
    C_A = C - C_M - C_T
  
    ax1.plot(ts, C, color='red', label='total', linewidth=1)
    ax1.plot(ts, C_M, color='blue', label='ocean', linewidth=1)
    ax1.plot(ts, C_T, color='green', label='land', linewidth=1)
    ax1.plot(ts, C_A, color='grey', label='atmosphere', linewidth=1)
    
    ax1.text(400, 0.19, 'atm. PB',color= "grey",fontsize=6)
    ax1.text(10, 0.29, 'ocean PB', color= 'blue',fontsize=6)
    ax1.text(400, 0.57, 'land PB', color='green',fontsize=6)     
    ax1.fill_between(ts,0.59, 1.6, facecolor='green',alpha=0.15, edgecolor="")
    ax1.fill_between(ts,0.31, 0, facecolor='blue',alpha=0.1, edgecolor="")
    ax1.fill_between(ts,0.21, 0, facecolor='grey',alpha=0.3, edgecolor="")
    ax1.axhline(y=0.21,xmin=0,xmax=500,c="grey",linestyle=":",linewidth=1, zorder=0, alpha=0.4)
    ax1.axhline(y=0.31,xmin=0,xmax=500,c="blue",linestyle=":",linewidth=1, zorder=0, alpha=0.4)
    ax1.axhline(y=0.59,xmin=0,xmax=500,c="green",linestyle=":",linewidth=1, zorder=0, alpha=0.4)
 
    ax1.spines["top"].set_linewidth(0.5)     
    ax1.spines["bottom"].set_linewidth(0.5)     
    ax1.spines["right"].set_linewidth(0.5)    
    ax1.spines["left"].set_linewidth(0.5)  
    ax1.get_xaxis().tick_bottom()    
    ax1.get_yaxis().tick_right() 

    ax1.set_ylabel('carbon pools', fontsize=8)
    ax1.set_xlim([0, 600])
    ax1.set_ylim([0, c_max+1.05])
    ax1.set_xticklabels(["",1900,"",2100,"",2300,""], fontsize=6)
    ax1.set_yticklabels(["",0.2,"",0.6,"",1.0,"",1.4], fontsize=6)
    ax1.text(0.03, 0.98, "a) no tCDR", transform=ax1.transAxes, fontsize=6, verticalalignment='top')
    
    ###################################################################################################
    alpha_max =alpha2
    
    #  solve the DEs
    soln = integrate.odeint(deriv, y0, ts)  # ODE solver
    C_R = soln[:, 0]
    C_M = soln[:, 1]
    C_T = soln[:, 2]
    C_G = soln[:, 3]
    C = C_0 + C_R - C_G
    C_A = C - C_M - C_T
    
    print 'H_geo2', max(H_geo(C_T, C_A)*3989)
    print 'H_geo2', max(C_G*3989)
    ax2.plot(ts, C, color='red', label='total', linewidth=1)
    ax2.plot(ts, C_M, color='blue', label='ocean', linewidth=1)
    ax2.plot(ts, C_T, color='green', label='land', linewidth=1)
    ax2.plot(ts, C_A, color='grey', label='atmosphere', linewidth=1)
    ax2.plot(ts, C_G, color='Orchid', label='GE-pool', linewidth=1)

    ax2.text(400, 0.19, 'atm. PB',color= "grey",fontsize=6)
    ax2.text(10, 0.29, 'ocean PB', color= 'blue',fontsize=6)
    ax2.text(400, 0.57, 'land PB', color='green',fontsize=6)     
    ax2.fill_between(ts,0.59, 1.6, facecolor='green',alpha=0.15, edgecolor="")
    ax2.fill_between(ts,0.31, 0, facecolor='blue',alpha=0.1, edgecolor="")
    ax2.fill_between(ts,0.21, 0, facecolor='grey',alpha=0.3, edgecolor="")
    ax2.axhline(y=0.21,xmin=0,xmax=500,c="grey",linestyle=":",linewidth=1, zorder=0, alpha=0.4)
    ax2.axhline(y=0.31,xmin=0,xmax=500,c="blue",linestyle=":",linewidth=1, zorder=0, alpha=0.4)
    ax2.axhline(y=0.59,xmin=0,xmax=500,c="green",linestyle=":",linewidth=1, zorder=0, alpha=0.4)
 
    ax2.spines["top"].set_linewidth(0.5)     
    ax2.spines["bottom"].set_linewidth(0.5)     
    ax2.spines["right"].set_linewidth(0.5)    
    ax2.spines["left"].set_linewidth(0.5)  
    ax2.get_xaxis().tick_bottom()    
    ax2.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off
  
    ax2.set_xlabel('time', fontsize=8)
    ax2.set_xlim([0, 600])
    ax2.set_xticklabels(["",1900,"",2100,"",2300,""], fontsize=6)
    ax2.text(0.03, 0.98, "b) intermediate tCDR rate", transform=ax2.transAxes, fontsize=6, verticalalignment='top')
    
    ###################################################################################################
    alpha_max = alpha3 
    #  solve the DEs
    soln = integrate.odeint(deriv, y0, ts)  # ODE solver
    C_R = soln[:, 0]
    C_M = soln[:, 1]
    C_T = soln[:, 2]
    C_G = soln[:, 3]
    C = C_0 + C_R - C_G
    C_A = C - C_M - C_T
    print 'H_geo3', max(H_geo(C_T, C_A)*3989)
    print 'H_geo3', max(C_G*3989)
    ax3.plot(ts, C, color='red', label='total', linewidth=1)
    ax3.plot(ts, C_M, color='blue', label='ocean', linewidth=1)
    ax3.plot(ts, C_T, color='green', label='land', linewidth=1)
    ax3.plot(ts, C_A, color='grey', label='atmosphere', linewidth=1)
    ax3.plot(ts, C_G, color='Orchid', label='CE sink', linewidth=1)
    
    ax3.text(400, 0.19, 'atm. PB',color= "grey",fontsize=6)
    ax3.text(10, 0.29, 'ocean PB', color= 'blue',fontsize=6)
    ax3.text(400, 0.57, 'land PB', color='green',fontsize=6)     
    ax3.fill_between(ts,0.59, 1.6, facecolor='green',alpha=0.15, edgecolor="")
    ax3.fill_between(ts,0.31, 0, facecolor='blue',alpha=0.1, edgecolor="")
    ax3.fill_between(ts,0.21, 0, facecolor='grey',alpha=0.3, edgecolor="")
    ax3.axhline(y=0.21,xmin=0,xmax=500,c="grey",linestyle=":",linewidth=1, zorder=0, alpha=0.4)
    ax3.axhline(y=0.31,xmin=0,xmax=500,c="blue",linestyle=":",linewidth=1, zorder=0, alpha=0.4)
    ax3.axhline(y=0.59,xmin=0,xmax=500,c="green",linestyle=":",linewidth=1, zorder=0, alpha=0.4)
 
    ax3.spines["top"].set_linewidth(0.5)     
    ax3.spines["bottom"].set_linewidth(0.5)     
    ax3.spines["right"].set_linewidth(0.5)    
    ax3.spines["left"].set_linewidth(0.5)  
    ax3.get_xaxis().tick_bottom()    
    ax3.get_yaxis().tick_right() 

    ax3.set_xlim([0, 600])
    ax3.set_xticklabels(["",1900,"",2100,"",2300,""], fontsize=6)
    ax3.set_yticklabels(["",0.2,"",0.6,"",1.0,"",1.4], fontsize=6)
    ax3.text(0.03, 0.98, "c) high tCDR rate", transform=ax3.transAxes, fontsize=6, verticalalignment='top')

    proxy1 = matplotlib.lines.Line2D([1], [0], linestyle="-", c='red')
    proxy2 = matplotlib.lines.Line2D([0], [0], linestyle="-", c='green')
    proxy3 = matplotlib.lines.Line2D([0], [0], linestyle="-", c='blue')
    proxy4 = matplotlib.lines.Line2D([0], [0], linestyle="-", c='grey')
    proxy5 = matplotlib.lines.Line2D([0], [0], linestyle="-", c='Orchid')
    f.legend((proxy1, proxy2,proxy3, proxy4,proxy5), 
        ('active carbon','land','upper ocean','atmosphere', 'CE sink'), 
        'upper center',fontsize=6,ncol=5,frameon=False)

    if saveToFile:
        f.savefig(saveToFile)
    plt.show()



