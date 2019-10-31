
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 17:25:19 2019

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
omega = np.arange(20, 200, 1)*1e3*2*np.pi # freq for spectrum
omega_minus = -1*omega # freq for spectrum
#omega_laser = 10 # freq of laser
omega_j = np.array([190, 170, 38])*1e3*2*np.pi # freq of mechanical modes

#detuning = omega_j[0] - omega_laser

detuning = 300e3 * 2*np.pi #Hz
omega_laser = 3e8/(1064e-9)*2*np.pi + detuning
print("laser freq: ", round(omega_laser*1e-9), 'GHz')

# damping
Gamma = 5e-7
kappa = 2*np.pi*193e3 #Hz

# coupling
g = np.array([20,30,7])*1e3 *2*np.pi # Hz ????

# phases
phi = np.array([0,0,np.pi])
#############################################################################

### resulting variables ###
n_opt = 0

n_mech = k*T/(hbar * omega_j)


    
####################
### PLOT EXAMPLE ###
####################
    
### omega plus

# calculate spectrum of q
param = (omega, omega_j, detuning, g, Gamma, kappa, n_opt, n_mech)
S_x = tools.spectrum_output(0, param)
S_y = tools.spectrum_output(1, param)
S_z = tools.spectrum_output(2, param)


# plot
plt.plot(omega/(2*np.pi)*1e-3, S_x, label = 'x')
plt.plot(omega/(2*np.pi)*1e-3, S_y, label = 'y')
plt.plot(omega/(2*np.pi)*1e-3, S_z, label = 'z')

plt.axvline(omega_j[0]/(2*np.pi)*1e-3, color='gray')
plt.axvline(omega_j[1]/(2*np.pi)*1e-3, color='gray')
plt.axvline(omega_j[2]/(2*np.pi)*1e-3, color='gray')


plt.yscale('log')
#plt.ylim([1e-25,1e3])

plt.xlabel('omega/2pi [kHz]')
plt.ylabel('spectrum')
plt.title('first example')
plt.legend(loc='best')

#plt.show()
plt.savefig('pic/omega_pos')
plt.close()
        

### omega minus
# calculate spectrum of q
param = (omega_minus, omega_j, detuning, g, Gamma, kappa, n_opt, n_mech)
S_x = tools.spectrum_output(0, param)
S_y = tools.spectrum_output(1, param)
S_z = tools.spectrum_output(2, param)


# plot
plt.plot(omega_minus/(2*np.pi)*1e-3, S_x, label = 'x')
plt.plot(omega_minus/(2*np.pi)*1e-3, S_y, label = 'y')
plt.plot(omega_minus/(2*np.pi)*1e-3, S_z, label = 'z')

plt.axvline(-omega_j[0]/(2*np.pi)*1e-3, color='gray')
plt.axvline(-omega_j[1]/(2*np.pi)*1e-3, color='gray')
plt.axvline(-omega_j[2]/(2*np.pi)*1e-3, color='gray')


plt.yscale('log')
#plt.ylim([1e-25,1e3])

plt.xlabel('omega/2pi [kHz]')
plt.ylabel('spectrum')
plt.title('first example')
plt.legend(loc='best')

#plt.show()
plt.savefig('pic/omega_neg')
        