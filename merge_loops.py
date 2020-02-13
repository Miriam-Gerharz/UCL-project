#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 11:00:59 2020

@author: miriam
"""

import numpy as np
import pickle
import matplotlib.pyplot as plt

def find_min(array, theta, detuning):
    #print(len(detuning))
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




Finesse = np.array([10, 73, 150])*1e3
runs = np.array([1, 2, 3, 4])

n_x = [[0 for i in range(len(Finesse))] for j in range(len(runs))]
n_y = [[0 for i in range(len(Finesse))] for j in range(len(runs))]
n_z = [[0 for i in range(len(Finesse))] for j in range(len(runs))]
n_xy = [[0 for i in range(len(Finesse))] for j in range(len(runs))]


print('run, Finesse, X0, press, n_x, n_y, n_z, n_x+n_y')

for i in range(len(Finesse)):
    for j in range(len(runs)):
        run = str(runs[j])
        _Finesse = str(Finesse[i])[:-2]
        path_file = 'files/3D/'
        name_file = 'loop_over_theta_and_detuning_with_Finesse_' + _Finesse + '_run_' + run
        filename = path_file + name_file
        
        f = open(filename + '.pckl', 'rb')
        param, omega, theta, detuning, _n_x, _n_y, _n_z = pickle.load(f)
        _n_xy = np.array(_n_x)+np.array(_n_y)

        _n_x = np.array(_n_x)
        _n_y = np.array(_n_y)
        _n_z = np.array(_n_z)
        _n_xy = np.array(_n_xy)
        
        #find minima
        n_x[j][i] = find_min(_n_x, theta, detuning)[0]
        n_y[j][i] = find_min(_n_y, theta, detuning)[0]
        n_z[j][i] = find_min(_n_z, theta, detuning)[0]
        n_xy[j][i] = find_min(_n_xy, theta, detuning)[0]
        
        print(run, _Finesse, ' & ', round(param.X0/1064e-9,3),' & ', round(param.Press*1e6,2),' & ', round(n_x[j][i],2),' & ', round(n_y[j][i],2),' & ', round(n_z[j][i],2),' & ', round(n_xy[j][i],2))
    
    
# plot

# X_0 = 0.05: runs 2,3

#for i in range(len(Finesse)):
#    _n_x_1 = []
#    for j in range(len(runs)):
#        _n_x = 

plt.plot(Finesse, n_x[1], linestyle = 'None', marker = 'x', color = 'blue', label = '1e-7, n_x')
plt.plot(Finesse, n_y[1], linestyle = 'None', marker = 'x', color = 'green', label = '1e-7, n_y')
plt.plot(Finesse, n_z[1], linestyle = 'None', marker = 'x', color = 'red', label = '1e-7, n_z')
plt.plot(Finesse, n_xy[1], linestyle = 'None', marker = 'x', color = 'black', label = '1e-7, n_x+n_y')
plt.plot(Finesse, n_x[2], linestyle = 'None', marker = 'x', color = 'cyan', label = '1e-6, n_x')
plt.plot(Finesse, n_y[2], linestyle = 'None', marker = 'x', color = 'lawngreen', label = '1e-6, n_y')
plt.plot(Finesse, n_z[2], linestyle = 'None', marker = 'x', color = 'orange', label = '1e-6, n_z')
plt.plot(Finesse, n_xy[2], linestyle = 'None', marker = 'x', color = 'grey', label = '1e-6, n_x+n_y')

plt.title('X_0 = 0.05 lambda')
plt.xlabel('Finesse')
plt.ylabel('n_i')
plt.yscale('log')
#plt.show()
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('pic/3D/loops/' +'merged_loops_05', bbox_inches='tight')
plt.close()

# X_0 = 0.23: runs 1,4
plt.plot(Finesse, n_x[0], linestyle = 'None', marker = 'x', color = 'blue', label = '1e-7, n_x')
plt.plot(Finesse, n_y[0], linestyle = 'None', marker = 'x', color = 'green', label = '1e-7, n_y')
plt.plot(Finesse, n_z[0], linestyle = 'None', marker = 'x', color = 'red', label = '1e-7, n_z')
plt.plot(Finesse, n_xy[0], linestyle = 'None', marker = 'x', color = 'black', label = '1e-7, n_x+n_y')
plt.plot(Finesse, n_x[3], linestyle = 'None', marker = 'x', color = 'cyan', label = '1e-6, n_x')
plt.plot(Finesse, n_y[3], linestyle = 'None', marker = 'x', color = 'lawngreen', label = '1e-6, n_y')
plt.plot(Finesse, n_z[3], linestyle = 'None', marker = 'x', color = 'orange', label = '1e-6, n_z')
plt.plot(Finesse, n_xy[3], linestyle = 'None', marker = 'x', color = 'grey', label = '1e-6, n_x+n_y')

plt.title('X_0 = 0.23 lambda')
plt.xlabel('Finesse')
plt.ylabel('n_i')
plt.yscale('log')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#plt.show()
plt.savefig('pic/3D/loops/' +'merged_loops_23', bbox_inches='tight')