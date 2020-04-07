#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 20:05:27 2019

@author: miriam
"""

import numpy as np
import matplotlib.pyplot as plt
import tools
   
param = tools.parameters()
omega = np.arange(0, 700, 10e-2)*1e3*2*np.pi # freq for spectrum


###########################################################################

### PARAMETERS ###
'''
# for high an low pressure 1D
param.Press = 1e-3
'''

'''
# for 3D spectra
param.X0 = 0.37 * param.lambda_tw
'''

'''
# for recovering experimental results
param.detuning = -315e3
'''


# for optimization, single phonon numbers
param.WX = 0.6e-6 #m, adjusted parameters
param.WY = 0.705e-6 #m, adjusted parameters

param.theta0 = 0.5
param.detuning = -125e3
param.X0 = 0.03 * param.lambda_tw

#######################

# calculation
consider_3D = True
calc_1D = False

# plotting
plot_1D_and_3D = True
log_plot = True
plot_z = False
legend = True
plot_range = [-400, 400] #kHz/2pi
figsize = (4,2.5)

# saving
pic = 'Fig1'
#filename = 'pic/3D/3D_setup_Delic_det_'+str(round(param.detuning*1e-3))+'kHz_' + str(round(param.X0/1064e-9*1e2)) + '_' + pic 
filename = 'pic/test'
###############################################################
if plot_1D_and_3D == True:
    consider_3D = True
    calc_1D = True

# Calculate the needed parameters and print all
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

if calc_1D == True:    
    # 1D calculations
    SXX_plus = tools.spectrum_output(omega, 0, param, False)
    SXX_minus = tools.spectrum_output(-omega, 0, param, False)
    SYY_plus = tools.spectrum_output(omega, 1, param, False)
    SYY_minus = tools.spectrum_output(-omega, 1, param, False)
    SZZ_plus = tools.spectrum_output(omega, 2, param, False)
    SZZ_minus = tools.spectrum_output(-omega, 2, param, False)
    
    plot_1D = True
    plot_3D = False
    style = '-'

# 3D calculations
if consider_3D == True:
    SXX_plus_3D = tools.spectrum_output(omega, 0, param, True)
    SXX_minus_3D = tools.spectrum_output(-omega, 0, param, True)
    SYY_plus_3D = tools.spectrum_output(omega, 1, param, True)
    SYY_minus_3D = tools.spectrum_output(-omega, 1, param, True)
    SZZ_plus_3D = tools.spectrum_output(omega, 2, param, True)
    SZZ_minus_3D = tools.spectrum_output(-omega, 2, param, True)
    
    plot_3D = True
    plot_1D = False


### PLOTTING ###
if plot_1D_and_3D == True:
    plot_1D = True
    plot_3D = True
    style = '--' #linestyle of 1D spectra

fig = plt.figure(figsize=figsize)    
  
# plot 3D parts  
if plot_3D == True:
    plt.plot(omega/2/np.pi*1e-3, SXX_plus_3D, color = 'red', label = 'x 3D', linestyle = '-')
    plt.plot(-omega/2/np.pi*1e-3, SXX_minus_3D, color = 'red', linestyle = '-')
    plt.plot(omega/2/np.pi*1e-3, SYY_plus_3D, color = 'blue', label = 'y 3D', linestyle = '-')
    plt.plot(-omega/2/np.pi*1e-3, SYY_minus_3D, color = 'blue', linestyle = '-')
    if plot_z ==True:
        plt.plot(omega/2/np.pi*1e-3, SZZ_plus_3D, color = 'green', label = 'z 3D', linestyle = '-')
        plt.plot(-omega/2/np.pi*1e-3, SZZ_minus_3D, color = 'green', label = 'z 3D', linestyle = '-')

# plot 1D parts
if plot_1D == True:    
    plt.plot(omega/2/np.pi*1e-3, SXX_plus, color = 'orange', label = 'x 1D', linestyle = style)
    plt.plot(-omega/2/np.pi*1e-3, SXX_minus, color = 'orange', linestyle = style)
    plt.plot(omega/2/np.pi*1e-3, SYY_plus, color = 'cyan', label = 'y 1D', linestyle = style)
    plt.plot(-omega/2/np.pi*1e-3, SYY_minus, color = 'cyan', linestyle = style)
    if plot_z ==True:
        plt.plot(omega/2/np.pi*1e-3, SZZ_plus, color = 'lawngreen', label = 'z 1D', linestyle = style)
        plt.plot(-omega/2/np.pi*1e-3, SZZ_minus, color = 'lawngreen', label = 'z 1D', linestyle = style)

# plotting settings
plt.xlabel('$\omega/(2\pi)$ [kHz]')
plt.ylabel('$S_{qq}$ [a.u.]')
if log_plot == True:
    plt.yscale('log')
plt.xlim(plot_range)
if legend == True:
    plt.legend(loc = 'best')
plt.savefig(filename, bbox_inches='tight')
plt.show()


# calculate photon numbers from area
Delta = np.abs(omega[1])-np.abs(omega[0])
if calc_1D == True:
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


