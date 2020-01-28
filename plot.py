#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 11:45:38 2020

@author: miriam
"""

import matplotlib.pyplot as plt
import pickle
import numpy as np

from mpl_toolkits.axes_grid1 import AxesGrid

####################

path_file = 'files/3D/'
name_file = 'loop_over_theta_and_detuning_with_Finesse_73000_run_2'
filename = path_file + name_file

print_x = True
print_y = True
print_z = True

x_title = '$\Delta$ [kHz/2pi]'
y_title = '$\\theta$ [$\pi$]'
fig_title = 'Finesse_' #+Finesse

###################

### LOAD DATA ###
f = open(filename + '.pckl', 'rb')
param, omega, theta, detuning, n_x, n_y, n_z = pickle.load(f)

n = 0
N = []
title = []
if print_x == True:
    n = n+1
    N.append(n_x)
    title.append('n_x')
if print_y == True:
    n = n+1
    N.append(n_y)
    title.append('n_y')
if print_z == True:
    n = n+1
    N.append(n_z)
    title.append('n_z')

fig = plt.figure(figsize=(6, 4))
fig.suptitle(fig_title + str(param.Finesse))

grid = AxesGrid(fig, 111,
                nrows_ncols=(1, n),
                axes_pad=0.05,
                cbar_mode='single',
                cbar_location='right',
                cbar_pad=0.1
                )
i = 0
for ax in grid:
    #ax.set_axis_off()
    ax.set_xlabel(x_title)
    ax.set_ylabel(y_title)
    ax.set_title(title[i])
    #ax.xaxis.set_xticks(detuning)
    
    ax.set_xticks(np.arange(0, len(detuning) +1, round(len(detuning))/10))
    ax.set_yticks(np.arange(0, len(theta) + 1, round(len(theta))/10))
    ax.set_xticklabels(detuning*1e-3, rotation = 90)
    ax.set_yticklabels(theta)
    im = ax.imshow(N[i], aspect = 'auto', origin = 'lower') #, extent = [detuning[0]*1e-3-0.5, detuning[-1]*1e-3+0.5, theta[0]-0.5, theta[-1]+0.5])
    i = i+1
    

# when cbar_mode is 'single', for ax in grid, ax.cax = grid.cbar_axes[0]

cbar = ax.cax.colorbar(im)
cbar = grid.cbar_axes[0].colorbar(im)

#cbar.ax.set_yticks(np.arange(0, 1.1, 0.5))
#cbar.ax.set_yticklabels(['low', 'medium', 'high'])
plt.show()

plt.savefig('pic/3D/' + name_file)

