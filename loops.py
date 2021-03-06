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

param.WX = 0.6e-6 #m, adjusted parameters
param.WY = 0.705e-6 #m, adjusted parameters
param.theta0 = 0.47 #pi, adjusted parameters
param.X0 = 0.03 * param.lambda_tw

# LOOP PARAMETER
theta = np.arange(0.00, 1.00001, 0.005) # [pi]
detuning = np.arange(-500, -100, 5)*1e3 #[Hz/2pi]

# run
run = '6'
filename = 'files/3D/loop_over_theta_and_detuning_with_Finesse_' + str(round(param.Finesse)) + '_run_' + run

#############################################################################

print('Finesse: ', str(param.Finesse), ', run: ', run)

# initialize arrays for n_i
n_x = [[0 for j in range(len(detuning))] for i in range(len(theta))]
n_y = [[0 for j in range(len(detuning))] for i in range(len(theta))]
n_z = [[0 for j in range(len(detuning))] for i in range(len(theta))]

# step with omega
Delta_omega = omega[1] - omega[0]

### CALCULATIONS ####
start = time.time()

for i in range(len(theta)):
    for j in range(len(detuning)):
        # Set parameters
        param.theta0 = theta[i]
        param.detuning = detuning[j]
        
        # Calculate resulting parameters
        param.prepare_calc()  
        param.Gamma = param.Gamma*2 # the calculated gamma is actually gamma/2

        # actually the detuning is given as an angular freq
        param.detuning = param.detuning *2*np.pi

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
        tools.loop_progress(len(detuning), len(theta), j, i, start)
        
### SAVE INTO FILE ###
to_save = [param, omega, theta, detuning, n_x, n_y, n_z]
f = open(filename + '.pckl', 'wb')

pickle.dump(to_save, f)
f.close()

# save settings of run
f = open('files/3D/overview_loop_over_theta_and_detuning', 'a')
f.write('\n********************************************************************')
f.write('\nFinesse: '+ str(param.Finesse) + ', run: ' + str(run))
f.write('\npressure: ' + str(param.Press))
f.write('\nX0: ' + str(param.X0/param.lambda_tw) + ' [lambda]')
f.write('\ntheta: ' + str(theta[0])+ ' ' + str(theta[-1])+ ' ' + str(theta[1]-theta[0]) + ' (start, end, step)')
f.write('\ndetuning: '+ str(detuning[0]*1e-3)+ ' ' + str(detuning[-1]*1e-3)+ ' '  + str(detuning[1]*1e-3-detuning[0]*1e-3) + ' (start, end, step)')
f.write('\n********************************************************************')
f.write('\n')

f.close()
