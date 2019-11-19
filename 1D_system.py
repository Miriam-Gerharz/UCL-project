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
hbar = 1.054571817e-34 #Js

###########################################################################
##################
### PARAMETERS ###
##################

# temperature 
T = 300 #K

# frequencies
omega = np.arange(0, 200, 1e-2)*1e3*2*np.pi # freq for spectrum
omega_j = np.array([174260, 157e3, 78e3])*2*np.pi # freq of mechanical modes
detuning = -300e3 * 2*np.pi #Hz

# damping
Gamma = 0.2317e-2*np.array([1, 1, 1])*2 # mechanical damping
kappa = 2*np.pi*93458 #Hz

# coupling
g = np.array([-25412,-39e3,57e3]) *2*np.pi # Hz ????

# phases
phi = np.array([0,0,np.pi])
#############################################################################

### resulting variables ###
n_opt = 0
n_mech = k*T/(hbar * omega_j)

#############################################################################

#################
### FUNCTIONS ###
#################



def print_parameter(param):
    print()
    print('Parameter:')
    print('mechanical frequencies/2pi: ', param[0]/2/np.pi)
    print('detuning/2pi: ', param[1]/2/np.pi)
    print('Gamma/2pi: ', param[3]/2/np.pi)
    print('kappa/2pi: ', param[4]/2/np.pi)
    print('couplings:')
    print('g_x, g_y, g_z: ', param[2]/2/np.pi, '*2pi')

########################################################################
    
############    
### MAIN ###
############


    
# optical damping rates
Gamma_opt = tools.opt_damp_rate(kappa, detuning, g, omega_j)

# photon numbers at equiv
N = tools.photon_number(n_mech, Gamma_opt, Gamma)


param = (omega_j, detuning, g, Gamma, kappa, n_opt, n_mech)
print_parameter(param)
    
SXX_plus = tools.spectrum_output(omega, 0, param)
SXX_minus = tools.spectrum_output(-omega, 0, param)
SYY_plus = tools.spectrum_output(omega, 1, param)
SYY_minus = tools.spectrum_output(-omega, 1, param)
SZZ_plus = tools.spectrum_output(omega, 2, param)
SZZ_minus = tools.spectrum_output(-omega, 2, param)

plt.plot(omega/2/np.pi*1e-3, SXX_plus, color = 'red', label = 'x')
plt.plot(-omega/2/np.pi*1e-3, SXX_minus, color = 'orange')
plt.plot(omega/2/np.pi*1e-3, SYY_plus, color = 'blue', label = 'y')
plt.plot(-omega/2/np.pi*1e-3, SYY_minus, color = 'cyan')
plt.plot(omega/2/np.pi*1e-3, SZZ_plus, color = 'green', label = 'z')
plt.plot(-omega/2/np.pi*1e-3, SZZ_minus, color = 'lawngreen')

plt.xlabel('$\omega/(2\pi)$ [Hz]')
plt.ylabel('$S_{ii}$ [a.u.]')
plt.legend(loc = 'best')
plt.savefig('pic/spectrum_1D_calc')
plt.show()


# calculate photon numbers from area
Delta = np.abs(omega[1])-np.abs(omega[0])
NX_from_area = tools.n_from_area(SXX_plus, SXX_minus, Delta, N[0], 'X')
NY_from_area = tools.n_from_area(SYY_plus, SYY_minus, Delta, N[1], 'Y')
NZ_from_area = tools.n_from_area(SZZ_plus, SZZ_minus, Delta, N[2], 'Z')

