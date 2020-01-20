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
hbar = 1.05e-34 #Js
c = 3e8 #m/s
grav = 9.8 #m/s^2
Epsi0=8.854e-12 # vacuum permitivity [F m^-1]
    

###########################################################################
##################
### PARAMETERS ###
##################


omega = np.arange(0, 300, 0.5e-2)*1e3*2*np.pi # freq for spectrum

settings = 'Tania'
#settings = 'Delic'

   
        
param = tools.parameters()

### setiings Tania ###
if settings == 'Tania':
    param.T = 300 #temperatur [K]
    param.R0=100e-9 #sphere radius
    param.RHO=2198 #sphere density
    param.EPSR=2.1 #parameter used to obtain the refractive index
    #              rho=2198.d0,EPSR=1.45d0**2,Epsi0=8.854d-12, &
    param.lambda_tw = 1064e-9 # wavelength [m]   
    #WK=5.9e6 #=2*pi/lambda=k
    param.waist=41.1e-6 #waist radius
    param.WX=0.67e-6
    param.WY=0.77e-6
    param.XL=1.07e-2 #cavity length 
    param.Finesse=15e4
    param.Press=1e-6 #air pressure in millibars
    param.Pin1=0.17e0 #input power in Watts tweezer beam
    param.detuning=-300e3 #detuning in KHz trap beam
    param.DelFSR=14e9 #1 FSR= 14 GHz, free spectral range =separation between cavity modes, you don't use these if you input g_x
    param.theta0=0.2 #angle between tweezer polarization and cavity axis. Given as as FRACTION of pi so pi/4  is theta0=0.25
    
    # Equilibrium positions 
    # X0=0.125*lambda ,Y0=waist/sqrt(2)
    param.Y0=0.0e-6
    #X0=0.125*1.064e-6
    param.X0=0.73*1.064e-6
    param.Z0=0e-6
    
    
    filename = 'pic/3D/3D_setup_Tania_det_'+str(round(param.detuning*1e-3))+'kHz' + str(round(param.X0/1064e-9*1e2)) + 'pi'

    #filename = 'pic/3D_setup_Tania_1'




### setiings Delic ###
if settings == 'Delic':
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
    param.detuning=-580e3 #detuning in KHz trap beam (omega_cav - omega_tw)
#    detuning=-380e3 #detuning in KHz trap beam (omega_cav - omega_tw)
#    detuning=-315e3 #detuning in KHz trap beam (omega_cav - omega_tw)
    param.DelFSR= 14e9 # 1FSR= 14 GHz, free spectral range =separation between cavity modes, you don't use these if you input g_x
    param.theta0= 0.25#angle between tweezer polarization and cavity axis. Given as as FRACTION of pi so pi/4  is theta0=0.25
    
    # Equilibrium positions 
    # X0=0.125*lambda ,Y0=waist/sqrt(2)
    param.Y0=0.0e-6
    #X0=0.125*1.064e-6
#    X0=0.25*1064e-9 + 3e-9 #x position [m]
    param.X0 = 0.23 * 1064e-9
#    X0 = 0.05*1064e-9
    param.Z0=0e-6
    
    filename = 'pic/3D/3D_setup_Delic_det_'+str(round(param.detuning*1e-3))+'kHz_' + str(round(param.X0/1064e-9*1e2)) 

###############################################################



param.prepare_calc()
param.print_param()




     

########################################################################
    
############    
### MAIN ###
############


    
# optical damping rates
Gamma_opt = param.opt_damp_rate()

# photon numbers at equiv
N = tools.photon_number(param.n_mech, Gamma_opt, param.Gamma)


#param2 = (param.omega_mech, param.detuning, param.g, param.Gamma, param.kappa, param.n_opt, param.n_mech)
#print_parameter(param)
    
# 1D calculations
SXX_plus = tools.spectrum_output(omega, 0, param, False)
SXX_minus = tools.spectrum_output(-omega, 0, param, False)
SYY_plus = tools.spectrum_output(omega, 1, param, False)
SYY_minus = tools.spectrum_output(-omega, 1, param, False)
SZZ_plus = tools.spectrum_output(omega, 2, param, False)
SZZ_minus = tools.spectrum_output(-omega, 2, param, False)

# 3D calculations
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
#plt.plot(omega/2/np.pi*1e-3, SZZ_plus, color = 'lawngreen', label = 'z 1D', linestyle = '--')
#plt.plot(omega/2/np.pi*1e-3, SYY_plus_3D, color = 'blue', label = 'y', linestyle = '-')
#plt.plot(-omega/2/np.pi*1e-3, SYY_minus_3D, color = 'cyan', linestyle = '-')
#plt.plot(omega/2/np.pi*1e-3, SZZ_plus_3D, color = 'green', label = 'z', linestyle = '-')
#plt.plot(-omega/2/np.pi*1e-3, SZZ_minus_3D, color = 'lawngreen', linestyle = '-')

plt.xlabel('$\omega/(2\pi)$ [Hz]')
plt.ylabel('$S_{qq}$ [a.u.]')
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


