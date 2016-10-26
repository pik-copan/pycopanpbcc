# -*- coding: utf-8 -*-
# Author: Vera Heck <heck@pik-potsdam.de>
# script designed to run parallel (step size given by number of nodes)


import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import matplotlib.cm as cm
import mpi 

array = np.array
# import code
import data_fig7 as areas
import main_model as main


# Define a function that may run on slaves:
def loop(par1, pars2):
    main.alpha_max = par1
    ii = 0
    a_vec = np.empty((nstep, 2 + 5))  # stores values for one alpha_max value and all alpha values
    for main.thresh_geo in pars2:
        a = areas.phase_space_area_3(main.deriv, tmax, n_tm, b_a, b_m, b_l)  
        a_vec[ii, ] = np.hstack([[main.alpha_max, main.thresh_geo], a])
        ii += 1
    return a_vec


# Define the master function that will run on the master:
def master():
    pars1 = np.linspace(0.0, 0.03, nstep)  # alpha_max
    pars2 = np.linspace(0.0, 0.3, nstep)  # thresh_geo

    for index in xrange(nstep):  # loop over alpha_max
        par1 = pars1[index]
        # Submit a job which will calculate partial sum 
        mpi.submit_call("loop", (par1, pars2), id=index)

    for i in xrange(nstep):
        print i
        a = mpi.get_result(id=i)
        if i == 0:
            a_alpha_max_thresh_geo = a
        else:
            a_alpha_max_thresh_geo = np.vstack([a_alpha_max_thresh_geo, a])

    print a_alpha_max_thresh_geo
    np.save('/save/a_alpha_max_thresh_geo_c_max=0.2', a_alpha_max_thresh_geo)

nstep = mpi.size  # steps of parameter variation
print nstep
tmax = 600  # length of time series for diff. eq. solution
n_tm = 100  # length of initial condition vector

# boundaries
b_a = 0.21
b_m = 0.31
b_l = 0.59

# model parameters 
main.c_max = 0.2 

mpi.run()



