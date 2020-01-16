
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
omega = np.arange(0, 200, 1e-2)*1e3*2*np.pi # freq for spectrum
#omega_laser = 10 # freq of laser
omega_j = np.array([174, 157, 78])*1e3*2*np.pi # freq of mechanical modes

#detuning = omega_j[0] - omega_laser

detuning = -300e3 * 2*np.pi #Hz
#omega_laser = 3e8/(1064e-9)*2*np.pi + detuning
#print("laser freq: ", round(omega_laser*1e-9), 'GHz')

# damping
Gamma = 0.5e-2*np.array([1, 1, 1]) # mechanical damping
kappa = 2*np.pi*93e3 #Hz
print('kapp2/2pi', kappa/4/np.pi)
print(Gamma)

# coupling
g = np.array([-25,-39,57])*1e3 *2*np.pi # Hz ????

# phases
phi = np.array([0,0,np.pi])
#############################################################################

### resulting variables ###
n_opt = 0
n_mech = k*T/(hbar * omega_j)

omega_minus = -1*omega # freq for spectrum


#################
### FUNCTIONS ###
#################

def opt_damp_rate(_kappa, _detuning, _g, _omega_j):
    return -_g**2 * _kappa * (1/(_kappa**2/4 + (_detuning+_omega_j)**2) - 1/(_kappa**2/4 + (_detuning-_omega_j)**2) )  

def photon_number(_n_j, _Gamma_opt, _Gamma_j):
    return _n_j * _Gamma_j / (abs(_Gamma_opt) + 2*_Gamma_j)

 
####################
### PLOT EXAMPLE ###
####################
    
### omega plus

# calculate spectrum of q
param = (omega, omega_j, detuning, g, Gamma, kappa, n_opt, n_mech)
#param2 = (omega, omega_j, detuning, g, Gamma*1e-6, kappa, n_opt, n_mech)
S_x = tools.spectrum_output(0, param)
#S_y = tools.spectrum_output(1, param)
#S_z = tools.spectrum_output(2, param)

#S_x2 = tools.spectrum_output(0, param2)
#print(S_x-S_x2)


# plot
plt.plot(omega/(2*np.pi)*1e-3, S_x, label = 'x')
#plt.plot(omega/(2*np.pi)*1e-3, S_y, label = 'y')
#plt.plot(omega/(2*np.pi)*1e-3, S_z, label = 'z')
#plt.plot(omega/(2*np.pi)*1e-3, S_x2, label = 'x2')

plt.axvline(omega_j[0]/(2*np.pi)*1e-3, color='gray')
#plt.axvline(omega_j[1]/(2*np.pi)*1e-3, color='gray')
#plt.axvline(omega_j[2]/(2*np.pi)*1e-3, color='gray')


plt.yscale('log')
#plt.ylim([1e-25,1e3])

plt.xlabel('omega/2pi [kHz]')
plt.ylabel('spectrum')
plt.title('first example')
plt.legend(loc='best')

#plt.show()
#plt.savefig('pic/omega_pos')
#plt.close()
        

# calculate the area under the curve
delta_omega = omega[1]-omega[0]
I_x = np.sum(S_x)*delta_omega
#I_y = np.sum(S_y)*delta_omega
#I_z = np.sum(S_z)*delta_omega

print('n_x (area +): ', round(I_x))
#print('n_y (area +): ', round(I_y))
#print('n_z (area +): ', round(I_z))


### omega minus
# calculate spectrum of q
param = (-omega, omega_j, detuning, g, Gamma, kappa, n_opt, n_mech)
S_x_minus = tools.spectrum_output(0, param)
#S_y = tools.spectrum_output(1, param)
#S_z = tools.spectrum_output(2, param)


# plot
plt.plot(omega/(2*np.pi)*1e-3, S_x_minus, label = 'x-')
#plt.plot(-omega_minus/(2*np.pi)*1e-3, S_y, label = 'y-')
#plt.plot(-omega_minus/(2*np.pi)*1e-3, S_z, label = 'z-')

#plt.axvline(-omega_j[0]/(2*np.pi)*1e-3, color='gray')
#plt.axvline(-omega_j[1]/(2*np.pi)*1e-3, color='gray')
#plt.axvline(-omega_j[2]/(2*np.pi)*1e-3, color='gray')


plt.yscale('log')
#plt.ylim([1e-25,1e3])

plt.xlabel('omega/2pi [kHz]')
plt.ylabel('spectrum')
plt.title('first example')
plt.legend(loc='best')

#plt.show()
plt.savefig('pic/omega_neg')
        

# calculate the area under the curve
delta_omega = omega[1]-omega[0]
I_x = np.sum(S_x_minus)*delta_omega
#I_y = np.sum(S_y)*delta_omega
#I_z = np.sum(S_z)*delta_omega

print('n_x (area -): ', round(I_x))
#print('n_y (area -): ', round(I_y))
#print('n_z (area -): ', round(I_z))




# optical damping rates
Gamma_opt_x = opt_damp_rate(kappa, detuning, g[0], omega_j[0])
Gamma_opt_y = opt_damp_rate(kappa, detuning, g[1], omega_j[1])
Gamma_opt_z = opt_damp_rate(kappa, detuning, g[2], omega_j[2])

print()
print('optical damping rates')
print('$\Gamma_{opt, x}$:', Gamma_opt_x/1e5, '1e5')
print('$\Gamma_{opt, y}$:', Gamma_opt_y)
print('$\Gamma_{opt, z}$:', Gamma_opt_z)

print('mechanical damping')
print('X', Gamma[0])

# photon numbers
n_theo_x = photon_number(n_mech[0], Gamma_opt_x, Gamma[0])
n_theo_y = photon_number(n_mech[1], Gamma_opt_y, Gamma[1])
n_theo_z = photon_number(n_mech[2], Gamma_opt_z, Gamma[2])

print('theoretical photon numbers at equiv')
print('n_x theo: ', round(n_theo_x, 1))
print('n_y theo: ', round(n_theo_y, 1))
print('n_z theo: ', round(n_theo_z, 1))

print('theoretical photon numbers at room')
print('n_x theo: ', round(n_mech[0]/1e8, 1), '1e8')
print('n_y theo: ', round(n_mech[1], 1))
print('n_z theo: ', round(n_mech[2], 1))