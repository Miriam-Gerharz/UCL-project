#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 17:55:53 2020

@author: miriam
"""


import numpy as np
import matplotlib.pyplot as plt

import tools

import time

#################
### CONSTANTS ###
#################

k = 1.380649e-23 #J/K
hbar = 1.054571817e-34 #Js
hbar = 1.05e-34 #Js
c = 3e8 #m/s
grav = 9.8 #m/s^2
Epsi0=8.854e-12 # vacuum permitivity [F m^-1]
    

###########################################################################
##################
### PARAMETERS ###
##################


omega = np.arange(0, 400, 0.5e-2)*1e3*2*np.pi # freq for spectrum

#settings = 'Tania'
#settings = 'Delic'

#consider_3D = False
#log_plot = False
   
        
param = tools.parameters()
#finess = np.array([1e3, 5e3, 1e4, 5e4, 1e5, 5e5])
detuning = np.linspace(-7e5, -1e5, 20)


### setiings Tania ###
#param.T = 300 #temperatur [K]
#param.R0=71.5e-9 #sphere radius
#param.Finesse=73e3
#param.Press=1e-6 #air pressure in millibars
#param.Pin1=0.4e0 #input power in Watts tweezer beam
param.theta0=0.47 #angle between tweezer polarization and cavity axis. Given as as FRACTION of pi so pi/4  is theta0=0.25
#param.X0=0.25*1.064e-6
    
# slightly adjust w_x and w_y (change default values)
param.WX = 0.6e-6 #m
param.WY = 0.705e-6 #m
    

filename = 'pic/3D/phonon_numbers_for_diff_detunings_3'
figsize = (6,3.5)

################################################################

    #filename = 'pic/3D_setup_Tania_1'

attrs = vars(param)
print(', \n'.join("%s: %s" % item for item in attrs.items()))

n_x_equation = np.zeros(len(detuning))
n_x_area = np.zeros(len(detuning))
n_y_equation = np.zeros(len(detuning))
n_y_area = np.zeros(len(detuning))
n_z_equation = np.zeros(len(detuning))
n_z_area = np.zeros(len(detuning))
n_x_area_3D = np.zeros(len(detuning))
n_y_area_3D = np.zeros(len(detuning))
n_z_area_3D = np.zeros(len(detuning))


start = time.time()

for i in range(len(detuning)):
    # nice progress control
    tools.loop_progress(len(detuning), 1, i, 0, start)
    
    # define pressure
    param.detuning = detuning[i]
    
    param.prepare_calc()
    
    #param.print_param()
        
    # optical damping rates
    Gamma_opt = param.opt_damp_rate()
    
    # photon numbers at equiv
    N = tools.photon_number(param.n_mech, Gamma_opt, param.Gamma, printing = False)
    
    param.Gamma = param.Gamma*2 # the calculated gamma is actually gamma/2
    
    #param = (param.omega_mech, param.detuning, param.g, param.Gamma, param.kappa, param.n_opt, param.n_mech)
    #print_parameter(param)
    
    # actually the detuning is given as an angular freq
    param.detuning = param.detuning *2*np.pi

    #param2 = (param.omega_mech, param.detuning, param.g, param.Gamma, param.kappa, param.n_opt, param.n_mech)
    #print_parameter(param)
    
    #param.Gamma = param.Gamma*2 # the calculated gamma is actually gamma/2
    
    ### 1D calculations ###
    SXX_plus = tools.spectrum_output(omega, 0, param, False)
    SXX_minus = tools.spectrum_output(-omega, 0, param, False)
    SYY_plus = tools.spectrum_output(omega, 1, param, False)
    SYY_minus = tools.spectrum_output(-omega, 1, param, False)
    SZZ_plus = tools.spectrum_output(omega, 2, param, False)
    SZZ_minus = tools.spectrum_output(-omega, 2, param, False)
    
    # calculate photon numbers from area
    Delta = np.abs(omega[1])-np.abs(omega[0])
    NX_from_area = tools.n_from_area(SXX_plus, SXX_minus, Delta, N[0], 'X', printing = False)[2]
    NY_from_area = tools.n_from_area(SYY_plus, SYY_minus, Delta, N[1], 'Y', printing = False)[2]
    NZ_from_area = tools.n_from_area(SZZ_plus, SZZ_minus, Delta, N[2], 'Z', printing = False)[2]
    
    # save photon numbers from area 1D
    n_x_area[i] = NX_from_area
    n_y_area[i] = NY_from_area
    n_z_area[i] = NZ_from_area
    
    ### 3D calculations ###
    SXX_plus = tools.spectrum_output(omega, 0, param, True)
    SXX_minus = tools.spectrum_output(-omega, 0, param, True)
    SYY_plus = tools.spectrum_output(omega, 1, param, True)
    SYY_minus = tools.spectrum_output(-omega, 1, param, True)
    SZZ_plus = tools.spectrum_output(omega, 2, param, True)
    SZZ_minus = tools.spectrum_output(-omega, 2, param, True)
    
    # calculate photon numbers from area
    Delta = np.abs(omega[1])-np.abs(omega[0])
    NX_from_area = tools.n_from_area(SXX_plus, SXX_minus, Delta, N[0], 'X', printing = False)[2]
    NY_from_area = tools.n_from_area(SYY_plus, SYY_minus, Delta, N[1], 'Y', printing = False)[2]
    NZ_from_area = tools.n_from_area(SZZ_plus, SZZ_minus, Delta, N[2], 'Z', printing = False)[2]
    
    # save photon numbers from area 3D
    n_x_area_3D[i] = NX_from_area
    n_y_area_3D[i] = NY_from_area
    n_z_area_3D[i] = NZ_from_area
    
    
    # save phonon numbers from formula
    n_x_equation[i] = N[0]
    n_y_equation[i] = N[1]
    n_z_equation[i] = N[2]
    
    # reset the detuning for the next loop
    param.detuning = param.detuning / (2*np.pi) 
    
    
    
    
    
#print(n_x_area, n_y_area, n_z_area)

fig = plt.figure(figsize = figsize)

plt.plot(detuning/1e3, n_x_equation, color = 'orange', label = 'x eq')
plt.plot(detuning/1e3, n_y_equation, color = 'cyan', label = 'y eq')
#plt.plot(detuning/1e3, n_z_equation, color = 'lawngreen', label = 'z equ', marker = 'x')
plt.plot(detuning/1e3, n_x_area, color = 'red', label = 'x area 1D', linestyle = '--')
plt.plot(detuning/1e3, n_y_area, color = 'blue', label = 'y area 1D',  linestyle = '--')
#plt.plot(detuning/1e3, n_x_area, color = 'green', label = 'z area', marker = 'x')
plt.plot(detuning/1e3, n_x_area_3D, color = 'red', label = 'x area 3D')
plt.plot(detuning/1e3, n_y_area_3D, color = 'blue', label = 'y area 3D')
#plt.plot(detuning/1e3, n_x_area_3D, color = 'green', label = 'z area', marker = 'x')
#plt.plot(detuning/1e3, n_y_area_3D/n_y_area, color = 'black', label = '3D / 1D (y)')

plt.hlines(y=1, xmin = detuning[0]/1e3, xmax = detuning[-1]/1e3, color = 'gray')
plt.xlabel('$\Delta$ [2$\pi$ kHz]')
plt.ylabel('Phonon numbers')
#plt.xscale('log')
plt.yscale('log')
plt.legend(loc = 'best')

plt.savefig(filename, bbox_inches='tight')


# save settings
f = open(filename+'_parameters', 'w')
f.write('T ' + str(param.T) + ' #temperatur [K] \n')
f.write('R0 ' + str(param.R0) + ' #sphere radius [m] \n')
f.write('RHO ' + str(param.RHO) + ' #sphere density [kg m^-3] \n')
f.write('EPSR ' + str(param.EPSR) + '# #parameter used to obtain the refractive index \n' )
f.write('lambda_tw ' + str(param.lambda_tw) + ' # wavelength [m] \n')
#f.write('WK ' + str(param.WK) + ' #=2*pi/lambda=k \n')
f.write('waist ' + str(param.waist) + ' #waist radius of cavity mode [m] \n')
f.write('WX ' + str(param.WX) + ' #waist size of tweezer [m] \n')
f.write('WY ' + str(param.WY) + '# waist size of tweezer [m] \n')
f.write('XL ' + str(param.XL) + ' #cavity length [m] \n')
#f.write('Finesse ' + str(param.Finesse) +' \n')
f.write('Press ' + str(param.Press) + ' #air pressure [mbar] \n')
f.write('Pin1 ' + str(param.Pin1) + ' #input power tweezer beam [W] \n')
f.write('detuning ' + str(param.detuning) + ' #detuning of trap beam (omega_cav - omega_tw) [kHz] \n')
f.write('DelFSR ' + str(param.DelFSR) + ' #1 FSR= 14 GHz, free spectral range =separation between cavity modes, you dont use these if you input g_x \n')
f.write('theta0 ' + str(param.theta0) + ' #angle between tweezer polarization and cavity axis. Given as as FRACTION of pi so pi/4  is theta0=0.25 \n')
f.write('Y0 ' + str(param.Y0) + ' #y position [m] \n')
f.write('X0 ' + str(param.X0) + ' #x position [m] \n')
f.write('Z0 ' + str(param.Z0) + ' #z position [m] \n')
f.close()

