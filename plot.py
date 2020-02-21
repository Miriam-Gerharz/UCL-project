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
run = '4'
Finesse = '73000'
#Finesse = '10000'
#Finesse = '150000'
path_file = 'files/3D/'
name_file = 'loop_over_theta_and_detuning_with_Finesse_' + Finesse + '_run_' + run
filename = path_file + name_file

# cannot have more than 3
print_x = True
print_y = True
print_xy = False
print_z = True


x_title = '$\Delta$ [kHz/2pi]'
y_title = '$\\theta$ [$\pi$]'
fig_title = 'Finesse_' #+Finesse

number_labels_x = 5
number_labels_y = 8

save = True
###################

### LOAD DATA ###
f = open(filename + '.pckl', 'rb')
param, omega, theta, detuning, n_x, n_y, n_z = pickle.load(f)
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
print('X0: ', param.X0/param.lambda_tw)

### PLOT ###

fig = plt.figure()
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
step_y = round(len(theta)/number_labels_y)
for ax in grid:
    #ax.set_axis_off()
    ax.set_xlabel(x_title)
    ax.set_ylabel(y_title)
    ax.set_title(title[i])
    #ax.xaxis.set_xticks(detuning)
    
    ax.set_xticks(np.arange(0, len(detuning) +1, step_x))
    ax.set_yticks(np.arange(0, len(theta) + 1, step_y))
    ax.set_xticklabels(detuning[::step_x]*1e-3, rotation = 90)
    ax.set_yticklabels(theta[::step_y])
    im = ax.imshow(N[i], aspect = 'auto', origin = 'lower', norm=matplotlib.colors.LogNorm()) #, extent = [detuning[0]*1e-3-0.5, detuning[-1]*1e-3+0.5, theta[0]-0.5, theta[-1]+0.5])
    i = i+1
    

cbar = grid.cbar_axes[0].colorbar(im)
cax = grid.cbar_axes[0]
axis = cax.axis[cax.orientation]
axis.toggle(ticks=True, ticklabels=True, label=True)

if save==False:
    plt.show()
else:
    plt.savefig('pic/3D/loops/' + name_file + '_pic_' + pic_number, bbox_inches='tight')
    
    
# print the minimal values    
def find_min(array):
    minimum = array[0][0]
    theta_min = 0
    detuning_min = 0
    for i in range(len(theta)):
        for j in range(len(detuning)):
            if array[i][j] < minimum:
                minimum = array[i][j]
                theta_min = theta[i]
                detuning_min = detuning[j]
    return minimum, theta_min, detuning_min

print('x', find_min(n_x))
print('y', find_min(n_y))
print('x+y', find_min(n_xy))
print('z', find_min(n_z))



