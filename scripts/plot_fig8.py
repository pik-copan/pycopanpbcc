# -*- coding: utf-8 -*-
# Author: Vera Heck <heck@pik-potsdam.de>
# Script generates Fig. 8 of Heck et al. 2016 (ESD)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pylab import *

array = np.array

nstep = 128   # steps of parameter variation 
par1='alpha_max'
par2= 'thresh_geo'
a2 = np.load('/save/a_alpha_max_thresh_geo_c_max=0.2.npy')
a4 = np.load('/save/a_alpha_max_thresh_geo_c_max=0.31.npy')
a6 = np.load('/save/a_alpha_max_thresh_geo_c_max=0.36.npy')
a8 = np.load('/save/a_alpha_max_thresh_geo_c_max=0.51.npy')
    
saveToFile = 'FIG8.pdf'
nrows, ncols = nstep, nstep
x = a2[:, 0]
y = a2[:, 1]

my_cmap = cm.get_cmap('Greens')
my_cmap.set_under('w')

fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, frameon=False, sharey=True, figsize=(18/2.54, 5/2.54))

plt.subplots_adjust(left=0.07, bottom=0.15, right=0.89, top=0.9, wspace=0.05, hspace=None)
area2 = array(a2[:, 2])
area4 = array(a4[:, 2])
area6 = array(a6[:, 2])
area8 = array(a8[:, 2])

titles = ['a) RCP2.6', 'b) RCP4.5', 'c) RCP6.0', 'd) RCP8.5']
            
grid1 = area2.reshape(nrows, ncols).T[::-1]
ax1.set_ylabel('tCDR threshold', fontsize=8)
ax1.set_title(titles[0], fontsize=8)

ax1.spines["top"].set_linewidth(0.5)     
ax1.spines["bottom"].set_linewidth(0.5)     
ax1.spines["right"].set_linewidth(0.5)    
ax1.spines["left"].set_linewidth(0.5)  

ax1.set_xticklabels(["",0.005,"",0.015,"",0.025,""], fontsize=7)
ax1.set_yticklabels([0,"",0.1,"",0.2,"",0.3], fontsize=7)
ax1.get_xaxis().tick_bottom()    
cax = ax1.imshow(grid1, extent=(x.min(), x.max(), y.min(), y.max()),interpolation='none', cmap=my_cmap, alpha=1,aspect ="auto")
cax.set_clim(vmin=(-0.0001), vmax=0.5)
plt.show()


grid2 = area4.reshape(nrows, ncols).T[::-1]
ax2.tick_params(labelsize=7)  # axis='x',
ax2.set_title(titles[1], fontsize=8)
ax2.set_xticklabels(["",0.005,"",0.015,"",0.025,""], fontsize=7)

ax2.spines["top"].set_linewidth(0.5)     
ax2.spines["bottom"].set_linewidth(0.5)     
ax2.spines["right"].set_linewidth(0.5)    
ax2.spines["left"].set_linewidth(0.5)  
ax2.text(0.025, -0.05, 'tCDR rate',color= "black",fontsize=8)
ax2.get_xaxis().tick_bottom()    

cax = ax2.imshow(grid2, extent=(x.min(), x.max(), y.min(), y.max()),interpolation='none', cmap=my_cmap, alpha=1, aspect ="auto")
cax.set_clim(vmin=(-0.0001), vmax=0.5)

grid3 = area6.reshape(nrows, ncols).T[::-1]
ax3.tick_params(labelsize=7)  
ax3.set_title(titles[2], fontsize=8)
ax3.set_xticklabels(["",0.005,"",0.015,"",0.025,""], fontsize=7)

ax3.spines["top"].set_linewidth(0.5)     
ax3.spines["bottom"].set_linewidth(0.5)     
ax3.spines["right"].set_linewidth(0.5)    
ax3.spines["left"].set_linewidth(0.5)  

ax3.get_xaxis().tick_bottom()    
cax = ax3.imshow(grid3, extent=(x.min(), x.max(), y.min(), y.max()), interpolation='none', cmap=my_cmap, alpha=1, aspect ="auto")
cax.set_clim(vmin=(-0.0001), vmax=0.5)

grid4 = area8.reshape(nrows, ncols).T[::-1]
ax4.tick_params(labelsize=7) 
ax4.set_title(titles[3], fontsize=8)
ax4.set_xticklabels(["",0.005,"",0.015,"",0.025,""], fontsize=7)
ax4.set_yticklabels([0,"",0.1,"",0.2,"",0.3], fontsize=7)

ax4.spines["top"].set_linewidth(0.5)     
ax4.spines["bottom"].set_linewidth(0.5)     
ax4.spines["right"].set_linewidth(0.5)    
ax4.spines["left"].set_linewidth(0.5)  

ax4.get_xaxis().tick_bottom() 
ax4.get_yaxis().tick_right() 
ax4.set_xlim([0, 0.03]) 
cax = ax4.imshow(grid4, extent=(x.min(), x.max(), y.min(), y.max()),vmax=0.32,interpolation='nearest', cmap=my_cmap, alpha=1, aspect ="auto")

cbar_ax = fig.add_axes([0.95, 0.155, 0.015, 0.75]) #location and width of cbar
cbar = fig.colorbar(cax, cax=cbar_ax, orientation='vertical',ticks=[0, 0.1, 0.2, 0.3]  )
cbar.set_clim(vmin=(-0.0001), vmax=0.5)
cbar.ax.tick_params(labelsize=8)

cbar.ax.set_title('MCSOS size', fontsize=7)
cbar.ax.set_yticklabels(['0', '0.1', '0.2', '0.3'],fontsize=7)

if saveToFile:
    plt.savefig(saveToFile)

