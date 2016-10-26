# -*- coding: utf-8 -*-
# Author: Vera Heck <heck@pik-potsdam.de>
# Script generates Fig. 7 of Heck et al. 2016 (ESD)

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import brewer2mpl
from scipy import interpolate

array = np.array

# import model code
import main_model as main


nstep = 150  # number of variation points (points form line)
s_geo = 200.

main.c_max= 0.4 
main.alpha_max = 0.02
areas_thresh = np.load('/save/var_thresh_alpha_max_'+ str(main.alpha_max)+ '_c_max_'+str(main.c_max) + '.dat.npy')

thresh = np.linspace(0, 0.3, nstep) 

xnew = np.linspace(thresh.min(),thresh.max(),300)

f_s = interp1d(thresh, areas_thresh[:, 0], kind='linear')
f_l = interp1d(thresh, areas_thresh[:, 3], kind='linear')
f_o = interp1d(thresh, areas_thresh[:, 2], kind='linear')
f_a = interp1d(thresh, areas_thresh[:, 1], kind='linear')
f_n = interp1d(thresh, areas_thresh[:, 4], kind='linear')
set2 = brewer2mpl.get_map('Set3', 'qualitative', 12).mpl_colors


f, (ax1, ax2) = plt.subplots(1, 2, frameon=False, sharey=True, figsize=(12/2.54, 5/2.54))
plt.subplots_adjust(left=None, bottom=0.2, right=None, top=0.9, wspace=0.05, hspace=None)
ax1.spines["top"].set_linewidth(0.5)     
ax1.spines["bottom"].set_linewidth(0.5)     
ax1.spines["right"].set_linewidth(0.5)    
ax1.spines["left"].set_linewidth(0.5)  
ax1.get_xaxis().tick_bottom()    
ax1.get_yaxis().tick_left()     

ax1.plot(xnew, f_s(xnew), '-',color=set2[6], linewidth=2, label="MCSOS")
ax1.plot(xnew, f_l(xnew), '-', color='Gold', label="Land MD")
ax1.plot(xnew, f_o(xnew), '-', color=set2[4], linewidth=0.6, label="Ocean MD")
ax1.plot(xnew, f_a(xnew), '-', color='slategrey', linewidth=0.6, label="Atmosphere MD")

ax1.set_ylabel('relative domain size', fontsize=8)
ax1.set_xlabel('tCDR threshold', fontsize=8)

ax1.set_xlim(0,xnew.max())
ax1.set_ylim([0, 0.97])
  
ax1.set_xticklabels(["","",0.1,"",0.2,"",0.3,""], fontsize=8)
ax1.set_yticklabels([0,0.2,0.4,0.6,0.8], fontsize=8)
ax1.text(0.03, 0.98, "a)", transform=ax1.transAxes, fontsize=8, verticalalignment='top')
    
#############################
main.c_max= 0.4 
main.thresh_geo = 0.2

areas_alpha_max = np.load('/save/var_alpha_max_thresh_'+ str(main.thresh_geo)+ '_c_max_'+ str(main.c_max) + '.dat.npy')

alpha_m = np.linspace(0, 0.03, nstep)

xnew2 = np.linspace(alpha_m.min(),alpha_m.max(),300)

f_s = interp1d(alpha_m, areas_alpha_max[:, 0], kind='linear')
f_l = interp1d(alpha_m, areas_alpha_max[:, 3], kind='linear')
f_o = interp1d(alpha_m, areas_alpha_max[:, 2], kind='linear')
f_a = interp1d(alpha_m, areas_alpha_max[:, 1], kind='linear')
f_n = interp1d(alpha_m, areas_alpha_max[:, 4], kind='linear')
set2 = brewer2mpl.get_map('Set3', 'qualitative', 12).mpl_colors


ax2.spines["top"].set_linewidth(0.5)     
ax2.spines["bottom"].set_linewidth(0.5)     
ax2.spines["right"].set_linewidth(0.5)    
ax2.spines["left"].set_linewidth(0.5)  
ax2.get_xaxis().tick_bottom()    
ax2.get_yaxis().tick_right()     


ax2.plot(xnew2, f_s(xnew2), '-',color=set2[6], linewidth=2, label="MCSOS")
ax2.plot(xnew2, f_l(xnew2), '-', color='Gold', label="Land MD")
ax2.plot(xnew2, f_o(xnew2), '-', color=set2[4], linewidth=0.6, label="Ocean MD")
ax2.plot(xnew2, f_a(xnew2), '-', color='slategrey', linewidth=0.6, label="Atmosphere MD")

ax2.set_xlabel('tCDR rate', fontsize=8)
ax2.set_xlim(0,xnew2.max())
ax2.set_ylim([0, 0.97])

ax2.set_xticklabels(["",0.005,"",0.015,"",0.025,""], fontsize=8)
ax2.set_yticklabels([0,0.2,0.4,0.6,0.8], fontsize=8)
ax2.text(0.03, 0.98, "b)", transform=ax2.transAxes, fontsize=8, verticalalignment='top')



proxy1 = matplotlib.lines.Line2D([1], [0], linestyle="-", linewidth=2,c=set2[6])
proxy2 = matplotlib.lines.Line2D([0], [0], linestyle="-",linewidth=2, c='Gold')
proxy3 = matplotlib.lines.Line2D([0], [0], linestyle="-", linewidth=2,c=set2[4])
proxy4 = matplotlib.lines.Line2D([0], [0], linestyle="-", linewidth=2,c='slategrey')
proxy5 = matplotlib.lines.Line2D([0], [0], linestyle="-", linewidth=2,c=set2[3])
f.legend((proxy1, proxy2,proxy3, proxy4),
      ('MCSOS','land MD','ocean MD','atmosphere MD'),
        'upper center',fontsize=6,ncol=4,frameon=False)

plt.show()
         
saveToFile = 'FIG7.pdf'
if saveToFile: plt.savefig(saveToFile)
