#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 17:55:53 2020

@author: miriam
"""


import numpy as np
import matplotlib.pyplot as plt

import tools

#################
### CONSTANTS ###
#################

k = 1.380649e-23 #J/K
hbar = 1.054571817e-34 #Js
#hbar = 1.05e-34 #Js
c = 3e8 #m/s
grav = 9.8 #m/s^2
Epsi0=8.854e-12 # vacuum permitivity [F m^-1]
    

###########################################################################
##################
### PARAMETERS ###
##################


omega = np.arange(0, 700, 10e-2)*1e3*2*np.pi # freq for spectrum

 ########################################################################  
        
param = tools.parameters()
pressures = np.logspace(-8, -3, 30, endpoint = True)

filename = 'pic/1D/phonon_numbers_for_diff_pressures'
figsize = (6,3.5)
    #filename = 'pic/3D_setup_Tania_1'


##########################################################################
attrs = vars(param)
print(', \n'.join("%s: %s" % item for item in attrs.items()))

n_x_equation = np.zeros(len(pressures))
n_x_area = np.zeros(len(pressures))
n_y_equation = np.zeros(len(pressures))
n_y_area = np.zeros(len(pressures))
n_z_equation = np.zeros(len(pressures))
n_z_area = np.zeros(len(pressures))

for i in range(len(pressures)):
    # define pressure
    param.Press = pressures[i]
    
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
    
    # 1D calculations
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
    
    # save phonon numbers
    n_x_equation[i] = N[0]
    n_y_equation[i] = N[1]
    n_z_equation[i] = N[2]
    
    n_x_area[i] = NX_from_area
    n_y_area[i] = NY_from_area
    n_z_area[i] = NZ_from_area
    
    # reset the detuning for the next loop
    param.detuning = param.detuning / (2*np.pi) 
    
    
#print(n_x_area, n_y_area, n_z_area)
fig = plt.figure(figsize = figsize)
plt.plot(pressures, n_x_area, color = 'red', label = 'x area')
plt.plot(pressures, n_y_area, color = 'blue', label = 'y area')
plt.plot(pressures, n_z_area, color = 'green', label = 'z area')
plt.plot(pressures, n_x_equation, color = 'orange', label = 'x eq', linestyle = '--')
plt.plot(pressures, n_y_equation, color = 'cyan', label = 'y eq', linestyle = '--')
plt.plot(pressures, n_z_equation, color = 'lawngreen', label = 'z eq', linestyle = '--')
plt.xlabel('Pressure [mbar]')
plt.ylabel('Phonon numbers')
plt.xscale('log')
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
f.write('Finesse ' + str(param.Finesse) +' \n')
#f.write('Press ' + str(param.Press) + ' #air pressure [mbar] \n')
f.write('Pin1 ' + str(param.Pin1) + ' #input power tweezer beam [W] \n')
f.write('detuning ' + str(param.detuning) + ' #detuning of trap beam (omega_cav - omega_tw) [kHz] \n')
f.write('DelFSR ' + str(param.DelFSR) + ' #1 FSR= 14 GHz, free spectral range =separation between cavity modes, you dont use these if you input g_x \n')
f.write('theta0 ' + str(param.theta0) + ' #angle between tweezer polarization and cavity axis. Given as as FRACTION of pi so pi/4  is theta0=0.25 \n')
f.write('Y0 ' + str(param.Y0) + ' #y position [m] \n')
f.write('X0 ' + str(param.X0) + ' #x position [m] \n')
f.write('Z0 ' + str(param.Z0) + ' #z position [m] \n')
f.close()


