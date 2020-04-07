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
param.theta0= 0.47 #pi, angle between tweezer polarization and cavity axis

# LOOP PARAMETER
x0 = np.arange(0.00, 1.0001, 0.05) # [pi]
detuning = np.arange(-500, -100, 5)*1e3 #[Hz/2pi]

run = '6'
filename = 'files/3D/loop_over_x0_and_detuning_with_Finesse_' + str(round(param.Finesse)) + '_run_' + run

#############################################################################

print('Finesse: ', str(param.Finesse), ', run: ', run)

# initialize arrays for n_i
n_x = [[0 for j in range(len(detuning))] for i in range(len(x0))]
n_y = [[0 for j in range(len(detuning))] for i in range(len(x0))]
n_z = [[0 for j in range(len(detuning))] for i in range(len(x0))]

# step with omega
Delta_omega = omega[1] - omega[0]

### CALCULATIONS ####
start = time.time()

for i in range(len(x0)):
    for j in range(len(detuning)):
        # Set parameters
        param.X0 = x0[i]
        param.detuning = detuning[j]
        
        # actually the detuning is given as an angular freq
        param.detuning = param.detuning *2*np.pi

        # Calculate resulting parameters
        param.prepare_calc()  
        param.Gamma = param.Gamma*2 # the calculated gamma is actually gamma/2

        # Calculate spectra
        SXX_plus_3D = tools.spectrum_output(omega, 0, param, True)
        SXX_minus_3D = tools.spectrum_output(-omega, 0, param, True)
        SYY_plus_3D = tools.spectrum_output(omega, 1, param, True)
        SYY_minus_3D = tools.spectrum_output(-omega, 1, param, True)
        SZZ_plus_3D = tools.spectrum_output(omega, 2, param, True)
        SZZ_minus_3D = tools.spectrum_output(-omega, 2, param, True)
        
        # Calculate n_i
        n_x[i][j] = tools.n_from_area(SXX_plus_3D, SXX_minus_3D, Delta_omega, _N=0, _name='', printing = False)[2]
        n_y[i][j] = tools.n_from_area(SYY_plus_3D, SYY_minus_3D, Delta_omega, _N=0, _name='', printing = False)[2]
        n_z[i][j] = tools.n_from_area(SZZ_plus_3D, SZZ_minus_3D, Delta_omega, _N=0, _name='', printing = False)[2]
        
        # print progress
        tools.loop_progress(len(detuning), len(x0), j, i, start)
        
### SAVE INTO FILE ###
to_save = [param, omega, x0, detuning, n_x, n_y, n_z]
f = open(filename + '.pckl', 'wb')
pickle.dump(to_save, f)
f.close()

# save settings of run
f = open('files/3D/overview_loop_over_x0_and_detuning', 'a')
f.write('\n********************************************************************')
f.write('\nFinesse: '+ str(param.Finesse) + ', run: ' + str(run))
f.write('\npressure: ' + str(param.Press))
f.write('\ntheta: ' + str(param.theta0))
f.write('\nX0: ' + str(x0[0])+ ' ' + str(x0[-1])+ ' ' + str(x0[1]-x0[0]) + ' (start, end, step)')
f.write('\ndetuning: '+ str(detuning[0]*1e-3)+ ' ' + str(detuning[-1]*1e-3)+ ' '  + str(detuning[1]*1e-3-detuning[0]*1e-3) + ' (start, end, step)')
f.write('\n********************************************************************')
f.write('\n')

f.close()
