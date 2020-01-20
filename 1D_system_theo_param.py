#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 20:05:27 2019

@author: miriam
"""

import numpy as np
import matplotlib.pyplot as plt

import tools

#################
### CONSTANTS ###
#################

k = 1.380649e-23 #J/K
k = 1.4e-23
hbar = 1.054571817e-34 #Js
hbar = 1.05e-34 #Js
c = 3e8 #m/s
grav = 9.8 #m/s^2

###########################################################################
##################
### PARAMETERS ###
##################

param = tools.parameters()

omega = np.arange(0, 300, 1e-2)*1e3*2*np.pi # freq for spectrum

# temperature 
param.T = 300 #K

# frequencies
#omega = np.arange(0, 200, 1e-2)*1e3*2*np.pi # freq for spectrum
#omega_laser = 10 # freq of laser
param.omega_mech = np.array([174, 157, 78])*1e3*2*np.pi # freq of mechanical modes

#detuning = omega_j[0] - omega_laser

param.detuning = -300e3 * 2*np.pi #Hz
#omega_laser = 3e8/(1064e-9)*2*np.pi + detuning
#print("laser freq: ", round(omega_laser*1e-9), 'GHz')

# damping
param.Gamma = 1e-8 # mechanical damping
param.kappa = 2*np.pi*93e3 #Hz
print('kapp2/2pi', param.kappa/4/np.pi)
print(param.Gamma)

# coupling
param.g = np.array([-25,-39,57, 0, 0, 0])*1e3 *2*np.pi # Hz ????
#param.g = np.array([-25,-39,57, 0, 0, 0])*1e1 *2*np.pi # Hz ????

# phases
#phi = np.array([0,0,np.pi])
#############################################################################

### resulting variables ###
#n_opt = 0
param.n_mech = k*param.T/(hbar * param.omega_mech)

omega_minus = -1*omega # freq for spectrum

filename = 'pic/1D/theo_param_spectrum_1D_calc_low'
    
############    
### MAIN ###
############
   
# optical damping rates
Gamma_opt = tools.opt_damp_rate(param.kappa, param.detuning, param.g, param.omega_mech)

# photon numbers at equiv
N = tools.photon_number(param.n_mech, Gamma_opt, param.Gamma)


#param = (omega_j, detuning, g, Gamma, kappa, n_opt, n_mech)
#print_parameter(param)
    
SXX_plus = tools.spectrum_output(omega, 0, param, False)
SXX_minus = tools.spectrum_output(-omega, 0, param, False)
SYY_plus = tools.spectrum_output(omega, 1, param, False)
SYY_minus = tools.spectrum_output(-omega, 1, param, False)
SZZ_plus = tools.spectrum_output(omega, 2, param, False)
SZZ_minus = tools.spectrum_output(-omega, 2, param, False)

plt.plot(omega/2/np.pi*1e-3, SXX_plus, color = 'red', label = 'x')
plt.plot(-omega/2/np.pi*1e-3, SXX_minus, color = 'red')
plt.plot(omega/2/np.pi*1e-3, SYY_plus, color = 'blue', label = 'y')
plt.plot(-omega/2/np.pi*1e-3, SYY_minus, color = 'blue')
plt.plot(omega/2/np.pi*1e-3, SZZ_plus, color = 'lawngreen', label = 'z')
plt.plot(-omega/2/np.pi*1e-3, SZZ_minus, color = 'lawngreen')

plt.xlabel('$\omega/(2\pi)$ [Hz]')
plt.ylabel('$S_{qq}$ [a.u.]')
plt.yscale('log')
plt.legend(loc = 'best')
plt.savefig(filename, bbox_inches='tight')
plt.show()


# calculate photon numbers from area
Delta = np.abs(omega[1])-np.abs(omega[0])
NX_from_area = tools.n_from_area(SXX_plus, SXX_minus, Delta, N[0], 'X')
NY_from_area = tools.n_from_area(SYY_plus, SYY_minus, Delta, N[1], 'Y')
NZ_from_area = tools.n_from_area(SZZ_plus, SZZ_minus, Delta, N[2], 'Z')

