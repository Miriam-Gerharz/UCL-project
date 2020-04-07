#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 10:52:53 2020

@author: miriam
"""

import tools
import numpy as np
import pickle
import time

###########################

param = tools.parameters()
omega = np.arange(0, 700, 10e-2)*1e3*2*np.pi # freq for spectrum


##################################################
### PARAMETER ###

# slightly adjust w_x and w_y (change default values)
param.WX = 0.6e-6 #m
param.WY = 0.705e-6 #m


# LOOP PARAMETER
theta = np.arange(0.2, 0.4001, 0.01) #[pi]
detuning = np.arange(-400, -200 +0.00001, 10)*1e3 #[Hz/2pi]
x0 = np.arange(0.15, 0.35001, 0.01) # [lambda]

# run
run = '3'
filename = 'files/3D/loop_over_theta_and_x0_and_detuning_with_Finesse_' + str(round(param.Finesse)) + '_run_' + run

#############################################################################

print('Finesse: ', str(param.Finesse), ', run: ', run)

# initialize arrays for minima
M = 1e10
name = ['n_x', 'n_y', 'n_z', 'n_x+n_y', 'n_x*n_y', 'n_x+n_z', 'n_x*n_z']
n_min = np.array([M, M, M, M, M, M, M]) #(n_x, n_y, n_z, n_xplusy, n_xy, n_xplusz, n_xz)
n_current = np.zeros(7)
theta_min = np.zeros(7)
x0_min = np.zeros(7)
det_min = np.zeros(7)

# step width omega
Delta_omega = omega[1] - omega[0]

### CALCULATIONS ####
start = time.time()
for k in range(len(theta)):
    for i in range(len(x0)):
        # print progress
        tools.loop_progress(len(x0), len(theta), i, k, start)
        for j in range(len(detuning)):
            # Set parameters
            param.X0 = x0[i]*param.lambda_tw
            param.detuning = detuning[j]
            param.theta0 = theta[k]
            
            # Calculate resulting parameters
            param.prepare_calc()  
            
            # actually the detuning is given as an angular freq
            param.detuning = param.detuning *2*np.pi
    
            #param.print_param()
            param.Gamma = param.Gamma*2 # the calculated gamma is actually gamma/2
    
            # Calculate spectra
            SXX_plus_3D = tools.spectrum_output(omega, 0, param, True)
            SXX_minus_3D = tools.spectrum_output(-omega, 0, param, True)
            SYY_plus_3D = tools.spectrum_output(omega, 1, param, True)
            SYY_minus_3D = tools.spectrum_output(-omega, 1, param, True)
            SZZ_plus_3D = tools.spectrum_output(omega, 2, param, True)
            SZZ_minus_3D = tools.spectrum_output(-omega, 2, param, True)
            
            # Calculate n_i
            n_current[0] = tools.n_from_area(SXX_plus_3D, SXX_minus_3D, Delta_omega, _N=0, _name='', printing = False)[2]
            n_current[1] = tools.n_from_area(SYY_plus_3D, SYY_minus_3D, Delta_omega, _N=0, _name='', printing = False)[2]
            n_current[2] = tools.n_from_area(SZZ_plus_3D, SZZ_minus_3D, Delta_omega, _N=0, _name='', printing = False)[2]
            
            n_current[3] = n_current[0] + n_current[1] # n_x + n_y
            n_current[4] = n_current[0] * n_current[1] # n_x * n_y
            n_current[5] = n_current[0] + n_current[2] # n_x + n_z
            n_current[6] = n_current[0] * n_current[2] # n_x * n_z
            
            for a in range(7):
                if n_current[a] < n_min[a]:
                    n_min[a] = n_current[a]
                    theta_min[a] = theta[k]
                    det_min[a] = detuning[j]
                    x0_min[a] = x0[i]
            
### PRINT RESULTS ###
print(' & n & $\\theta$ & $\Delta$ & $x_0$')
for i in range(7):
    print('$', name[i],'$ & ', round(n_min[i], 2), ' & ', theta_min [i], ' & ', det_min[i]/1e3, ' & ', x0_min[i])


### SAVE INTO FILE ###
to_save = [param, omega, name, n_min, theta_min, det_min, x0_min]
f = open(filename + '.pckl', 'wb')
pickle.dump(to_save, f)
f.close()
