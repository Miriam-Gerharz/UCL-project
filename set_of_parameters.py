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
#hbar = 1.05e-34 #Js
c = 3e8 #m/s
grav = 9.8 #m/s^2
Epsi0=8.854e-12 # vacuum permitivity [F m^-1]
    
param = tools.parameters()

###########################################################################
##################
### PARAMETERS ###
##################


omega = np.arange(0, 400, 0.5e-2)*1e3*2*np.pi # freq for spectrum

consider_3D = True
log_plot = True
plot_z = False
       

### setiings ###
param.T = 300 #temperatur [K]
param.R0= 0.5*143e-9#sphere radius [m]
param.Finesse = 73e3
param.Press=1e-7 #air pressure in millibars
param.Pin1= 400e-3 #input power in Watts tweezer beam
param.detuning=-300e3 #detuning in KHz trap beam (omega_cav - omega_tw)
param.theta0= 0.75#angle between tweezer polarization and cavity axis. Given as as FRACTION of pi so pi/4  is theta0=0.25
param.X0 = 0.23*1064e-9  # 5.32e-8 #0.05 * 1064e-9

pic = '1'
#filename = 'pic/3D/3D_setup_Delic_det_'+str(round(param.detuning*1e-3))+'kHz_' + str(round(param.X0/1064e-9*1e2)) + '_' + pic 
#filename = 'pic/low_pressure'
filename = 'pic/test'
###############################################################


param.prepare_calc()

param.print_param()
    
# optical damping rates
Gamma_opt = param.opt_damp_rate()

# photon numbers at equiv
N = tools.photon_number(param.n_mech, Gamma_opt, param.Gamma)

# the calculated gamma is actually gamma/2
param.Gamma = param.Gamma*2 

# actually the detuning is given as an angular freq
param.detuning = param.detuning *2*np.pi

    
# 1D calculations
SXX_plus = tools.spectrum_output(omega, 0, param, False)
SXX_minus = tools.spectrum_output(-omega, 0, param, False)
SYY_plus = tools.spectrum_output(omega, 1, param, False)
SYY_minus = tools.spectrum_output(-omega, 1, param, False)
SZZ_plus = tools.spectrum_output(omega, 2, param, False)
SZZ_minus = tools.spectrum_output(-omega, 2, param, False)

# 3D calculations
if consider_3D == True:
    SXX_plus_3D = tools.spectrum_output(omega, 0, param, True)
    SXX_minus_3D = tools.spectrum_output(-omega, 0, param, True)
    SYY_plus_3D = tools.spectrum_output(omega, 1, param, True)
    SYY_minus_3D = tools.spectrum_output(-omega, 1, param, True)
    SZZ_plus_3D = tools.spectrum_output(omega, 2, param, True)
    SZZ_minus_3D = tools.spectrum_output(-omega, 2, param, True)

    # plotting
    plt.plot(omega/2/np.pi*1e-3, SXX_plus_3D, color = 'red', label = 'x 3D', linestyle = '-')
    plt.plot(-omega/2/np.pi*1e-3, SXX_minus_3D, color = 'red', linestyle = '-')
    plt.plot(omega/2/np.pi*1e-3, SYY_plus_3D, color = 'blue', label = 'x 3D', linestyle = '-')
    plt.plot(-omega/2/np.pi*1e-3, SYY_minus_3D, color = 'blue', linestyle = '-')

#plt.plot(-omega/2/np.pi*1e-3, SZZ_minus, color = 'lawngreen', linestyle = '--')
plt.plot(omega/2/np.pi*1e-3, SXX_plus, color = 'orange', label = 'x 1D', linestyle = '--')
plt.plot(-omega/2/np.pi*1e-3, SXX_minus, color = 'orange', linestyle = '--')
plt.plot(omega/2/np.pi*1e-3, SYY_plus, color = 'cyan', label = 'y 1D', linestyle = '--')
plt.plot(-omega/2/np.pi*1e-3, SYY_minus, color = 'cyan', linestyle = '--')
if plot_z ==True:
    plt.plot(omega/2/np.pi*1e-3, SZZ_plus, color = 'lawngreen', label = 'z 1D', linestyle = '--')
    plt.plot(-omega/2/np.pi*1e-3, SZZ_minus, color = 'lawngreen', label = 'z 1D', linestyle = '--')
#plt.plot(omega/2/np.pi*1e-3, SYY_plus_3D, color = 'blue', label = 'y', linestyle = '-')
#plt.plot(-omega/2/np.pi*1e-3, SYY_minus_3D, color = 'cyan', linestyle = '-')
#plt.plot(omega/2/np.pi*1e-3, SZZ_plus_3D, color = 'green', label = 'z', linestyle = '-')
#plt.plot(-omega/2/np.pi*1e-3, SZZ_minus_3D, color = 'lawngreen', linestyle = '-')

plt.xlabel('$\omega/(2\pi)$ [Hz]')
plt.ylabel('$S_{qq}$ [a.u.]')
if log_plot == True:
    plt.yscale('log')
#plt.title('3D calculation for $\phi=$'+str(round(Wkx0/2/np.pi,2))+'$\pi$')
plt.legend(loc = 'best')
plt.savefig(filename, bbox_inches='tight')
plt.show()


# calculate photon numbers from area
Delta = np.abs(omega[1])-np.abs(omega[0])
print('***1D spectra***')
NX_from_area = tools.n_from_area(SXX_plus, SXX_minus, Delta, N[0], 'X')
NY_from_area = tools.n_from_area(SYY_plus, SYY_minus, Delta, N[1], 'Y')
NZ_from_area = tools.n_from_area(SZZ_plus, SZZ_minus, Delta, N[2], 'Z')

if consider_3D == True:
    print('\n ***3D spectra (ignore difference)***')
    NX_from_area_3D = tools.n_from_area(SXX_plus_3D, SXX_minus_3D, Delta, N[0], 'X')
    NY_from_area_3D = tools.n_from_area(SYY_plus_3D, SYY_minus_3D, Delta, N[1], 'Y')
    NZ_from_area_3D = tools.n_from_area(SZZ_plus_3D, SZZ_minus_3D, Delta, N[2], 'Z')


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
f.write('Press ' + str(param.Press) + ' #air pressure [mbar] \n')
f.write('Pin1 ' + str(param.Pin1) + ' #input power tweezer beam [W] \n')
f.write('detuning ' + str(param.detuning) + ' #detuning of trap beam (omega_cav - omega_tw) [kHz] \n')
f.write('DelFSR ' + str(param.DelFSR) + ' #1 FSR= 14 GHz, free spectral range =separation between cavity modes, you dont use these if you input g_x \n')
f.write('theta0 ' + str(param.theta0) + ' #angle between tweezer polarization and cavity axis. Given as as FRACTION of pi so pi/4  is theta0=0.25 \n')
f.write('Y0 ' + str(param.Y0) + ' #y position [m] \n')
f.write('X0 ' + str(param.X0) + ' #x position [m] \n')
f.write('Z0 ' + str(param.Z0) + ' #z position [m] \n')
f.close()


