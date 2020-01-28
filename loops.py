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

#################
### CONSTANTS ###
k = 1.380649e-23 #J/K
hbar = 1.054571817e-34 #Js
hbar = 1.05e-34 #Js
c = 3e8 #m/s
grav = 9.8 #m/s^2
Epsi0=8.854e-12 # vacuum permitivity [F m^-1]
###########################

param = tools.parameters()
omega = np.arange(0, 300, 0.5e-2)*1e3*2*np.pi # freq for spectrum


##################################################
### PARAMETER ###

# EXPERIMENTAL PARAMETER
param.T = 300 #temperatur [K]
param.R0= 0.5*143e-9#sphere radius [m]
param.RHO= 2198 #sphere density [kg m^-3]
param.EPSR= 2.1  #parameter used to obtain the refractive index
#    EPSR= 1.5  #parameter used to obtain the refractive index
#              rho=2198.d0,EPSR=1.45d0**2,Epsi0=8.854d-12, &
param.lambda_tw = 1064e-9 # wavelength [m]
param.waist=41.1e-6  #waist radius of cavity mode
param.WX= 0.67e-6# waist size of tweezer [m]
param.WY= 0.77e-6 # waist size of tweezer [m]
param.XL=  1.07e-2 #cavity length [m]
param.Finesse = 73e3
#Finesse = 15e4
param.Press=1e-6 #air pressure in millibars
param.Pin1= 400e-3 #input power in Watts tweezer beam
#Pin1= 0.17 #input power in Watts tweezer beam
#param.detuning=-580e3 #detuning in KHz trap beam (omega_cav - omega_tw)
#    detuning=-380e3 #detuning in KHz trap beam (omega_cav - omega_tw)
#    detuning=-315e3 #detuning in KHz trap beam (omega_cav - omega_tw)
param.DelFSR= 14e9 # 1FSR= 14 GHz, free spectral range =separation between cavity modes, you don't use these if you input g_x
#param.theta0= 0.25#angle between tweezer polarization and cavity axis. Given as as FRACTION of pi so pi/4  is theta0=0.25

# Equilibrium positions 
# X0=0.125*lambda ,Y0=waist/sqrt(2)
param.Y0=0.0e-6
#X0=0.125*1.064e-6
#    X0=0.25*1064e-9 + 3e-9 #x position [m]
param.X0 = 0.23 * 1064e-9
#    X0 = 0.05*1064e-9
param.Z0=0e-6


# LOOP PARAMETER
theta = np.arange(0.05, 0.95, 0.025) # [pi]
detuning = np.arange(-600, -100, 25)*1e3 #[Hz/2pi]

run = '3'
filename = 'files/3D/loop_over_theta_and_detuning_with_Finesse_' + str(round(param.Finesse)) + '_run_' + run

#############################################################################

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

        # Calculate spectra
        SXX_plus_3D = tools.spectrum_output(omega, 0, param, True)
        SXX_minus_3D = tools.spectrum_output(-omega, 0, param, True)
        SYY_plus_3D = tools.spectrum_output(omega, 1, param, True)
        SYY_minus_3D = tools.spectrum_output(-omega, 1, param, True)
        SZZ_plus_3D = tools.spectrum_output(omega, 2, param, True)
        SZZ_minus_3D = tools.spectrum_output(-omega, 2, param, True)
        
        # Calculate n_i
        n_x[i][j] = tools.n_from_area(SXX_plus_3D, SXX_minus_3D, Delta_omega, printing = False)[2]
        n_y[i][j] = tools.n_from_area(SYY_plus_3D, SYY_minus_3D, Delta_omega, printing = False)[2]
        n_z[i][j] = tools.n_from_area(SZZ_plus_3D, SZZ_minus_3D, Delta_omega, printing = False)[2]
        
        # print progress
        tools.loop_progress(len(detuning), len(theta), j, i, start)
        
### SAVE INTO FILE ###
to_save = [param, omega, theta, detuning, n_x, n_y, n_z]
f = open(filename + '.pckl', 'wb')

pickle.dump(to_save, f)
f.close()

