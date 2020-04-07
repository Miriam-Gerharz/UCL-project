#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 11:45:38 2020

@author: miriam
"""

import matplotlib.pyplot as plt
import pickle
import numpy as np

import matplotlib.colors

from mpl_toolkits.axes_grid1 import AxesGrid

####################
pic_number = '1'
run = '6'
Finesse = '73000'
path_file = 'files/3D/'
name_file = 'loop_over_x0_and_detuning_with_Finesse_' + Finesse + '_run_' + run
filename = path_file + name_file

# cannot have more than 3
print_x = True
print_y = True
print_xy = False
print_z = True


x_title = '$\Delta$ [2$\pi$ kHz]'
y_title = '$x_0$ [$\lambda$]'
fig_title = 'Finesse_' #+Finesse
print_title = False


number_labels_x = 5
number_labels_y = 5
cbar_ticks = [0.01, 1, 10, 100, 1000, 1e5, 1e6]

save = True
figsize = (4,2.5)
###################

### LOAD DATA ###
f = open(filename + '.pckl', 'rb')
param, omega, x0, detuning, n_x, n_y, n_z = pickle.load(f)
n_xy = np.array(n_x)+np.array(n_y)
   
# Append all data that is to be plotted 
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
if print_xy == True:
    n = n+1
    N.append(n_xy)
    title.append('n_x + n_y')
if print_z == True:
    n = n+1
    N.append(n_z)
    title.append('n_z')
    
print('pressure: ', param.Press)
print('theta: ', param.theta0)
print('Wx: ', param.WX)
print('Wx: ', param.WY)

### PLOT ###
fig = plt.figure(figsize = figsize)
if print_title == True:
    fig.suptitle(fig_title + str(param.Finesse))


grid = AxesGrid(fig, 111,
                nrows_ncols=(1, n),
                axes_pad=0.05,
                cbar_mode='single',
                cbar_location='right',
                cbar_pad=0.1
                )
i = 0

step_x = round(len(detuning)/number_labels_x)
step_y = round(len(x0)/number_labels_y)
for ax in grid:
    ax.set_xlabel(x_title)
    ax.set_ylabel(y_title)
    ax.set_title(title[i])
    
    ax.set_xticks(np.arange(0, len(detuning) +1, step_x))
    ax.set_yticks(np.arange(0, len(x0) + 1, step_y))
    ax.set_xticklabels(detuning[::step_x]*1e-3, rotation = 90)
    ax.set_yticklabels(x0[::step_y])
    im = ax.imshow(N[i], aspect = len(detuning)/len(x0), origin = 'lower', norm=matplotlib.colors.LogNorm()) #, extent = [detuning[0]*1e-3-0.5, detuning[-1]*1e-3+0.5, theta[0]-0.5, theta[-1]+0.5])
    i = i+1
    

cbar = grid.cbar_axes[0].colorbar(im, ticks=cbar_ticks)

if save==False:
    plt.show()
else:
    plt.savefig('pic/3D/loops/' + name_file + '_pic_' + pic_number, bbox_inches='tight')
    
