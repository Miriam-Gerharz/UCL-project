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

#####################################

omega = np.arange(0, 400, 1e-2)*1e3*2*np.pi # freq for spectrum

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
param.Press=1e-5 #air pressure in millibars
param.Pin1=0.17e0 #input power in Watts tweezer beam
param.detuning=-300e3 #detuning in KHz/2pi trap beam
param.DelFSR=14e9 #1 FSR= 14 GHz, free spectral range =separation between cavity modes, you don't use these if you input g_x
param.theta0=0.2 #angle between tweezer polarization and cavity axis. Given as as FRACTION of pi so pi/4  is theta0=0.25

# Equilibrium positions 
# X0=0.125*lambda ,Y0=waist/sqrt(2)
param.Y0=0.0e-6
#X0=0.125*1.064e-6
param.X0=0.125*1.064e-6
param.Z0=0e-6


###############################################################





#################
### FUNCTIONS ###
#################


'''
def print_parameter(param):
    print()
    print('Parameter:')
    print('mechanical frequencies/2pi: ', param[0]/2/np.pi)
    print('detuning/2pi: ', param[1]/2/np.pi)
    print('Gamma/2pi: ', param[3]/2/np.pi)
    print('kappa/2pi: ', param[4]/2/np.pi)
    print('couplings:')
    print('g_x, g_y, g_z: ', param[2]/2/np.pi, '*2pi')
    
'''
    
############################
### PREPARE CALCULATIONS ###
############################                
      
'''   
# zero eq. initial values
XM=RHO*4.*np.pi/3.*R0**3
print('Mass of bead= (Kg)', XM)

Polaris=4.*np.pi*Epsi0*(EPSR-1)/(EPSR+2)*R0**3
print('Polarisability of bead=', Polaris)

OMOPT=c*WK
# note waist is waist radius
W2=waist**2

# add a factor of 4 here. Not in GALAX1-5 routines!!!
VOL=XL*np.pi*W2/4

# &&&&&&&&&&&&&
# Depth of cavity field. Weak and unimportant for CS case
A=OMOPT*Polaris/2./VOL/Epsi0
print('Cavity trap A/(2pi)=  (Hz); zero-th shift',  A/np.pi/2, A*np.cos(WK*X0)**2)
#&&&&&&&&&&&&&&&&&&
KAPPin=np.pi*c/Finesse/XL
print('kappaIn=  (Hz)', KAPPin/2/np.pi)

epsTW=4*Pin1/(WX*WY*np.pi*c*Epsi0)
epsTW = np.sqrt(epsTW)
epsCAV=hbar*OMOPT/(2.*VOL*Epsi0)
epsCAV = np.sqrt(epsCAV)
print('epsTW,epsCAV', epsTW,epsCAV)

coeff=WK*Polaris/Epsi0/OMOPT**2/np.pi
kappnano=4*coeff**2*DelFSR*np.cos(WK*X0)*np.cos(WK*X0)
print('kappnano',kappnano)
kappa=kappnano+KAPPin
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
ZR=WX*WY*WK/2
print('ZR=',ZR)
 
# Pressure 1.d-4 mBar => ~ 0.125Hz in old expt
# now take usual expression eg Levitated review by Li Geraci etc
# 1 bar= 10^ 5 pascal; Press is in mbar = 10^ 2 pascal
#gamma=16 * P/(np.pi*v*RHO*R)
# v=speed of air=500 /s
GAMMAM=1600*Press/np.pi
GAMMAM=GAMMAM/500/RHO/R0
#Fix of Feb.2016 our GAMMAM => GAMMAM/2!!
# gamma/2=8 P/(pi*v*rho*R)
Gamma=GAMMAM/2
print('mechanical damping* 2pi', Gamma)

print('theta [pi]', theta0)
#*************************************************************************
#*************************************************************************
# subroutine below obtains the optomechanical parameters
# SUBROUTINE EQUIL_PAR(thet,Wkx0,Polaris,epsTW,epsCAV,XM,ZR,kappnano,kappin,GammaM,OMX,OMY,OMZ,GMAT,PHON)
            
#************************************************************************
#************************************************************************
iwrite=6
Det2pi=detuning*2*np.pi

kappa=kappnano+KAPPin
kapp2=0.5*kappa
print('kappa/2/pi (kHz)=',kappa/2/np.pi)
#     note that detunings include zeroth order correction for linearised versions
#   as a first pass work out frequencies with Wk0=0
Wkx0=WK*X0 #was commented out
OmX=Polaris*epsTW**2/XM/WX**2
OmY=Polaris*epsTW**2/XM/WY**2
OmZ=0.5*Polaris*epsTW**2/XM/ZR**2

print('Wkxeq',Wkx0)
# OPtomechanical drive frequency
print('epsTW,epsCAV', epsTW,epsCAV)

# theta vs thet
thet = theta0 * np.pi

# Sept 5 we will use negative Edip
Edip=-0.5*Polaris*epsTW*epsCAV*np.sin(thet)
Ediph=Edip/hbar
print('Edip/2/pi/hbar=', Ediph/2/np.pi)
# photon number in cavity
# real part of photon field
ALPRe=Det2pi*Ediph*np.cos(Wkx0)/(kapp2**2+Det2pi**2)
ALPim=-kapp2*Ediph*np.cos(Wkx0)/(kapp2**2+Det2pi**2)
Nphoton=Ediph*Ediph*np.cos(Wkx0)*np.cos(Wkx0)
Nphoton=Nphoton/(kapp2**2+Det2pi**2)
print('delta,kappa/2,number of photons in cavity', Det2pi/2/np.pi,kapp2/2/np.pi,Nphoton)

### ADD CAVITY CORRECTION to frequency squared ###
C1=-Edip/XM*2.*ALPRe*WK**2*np.cos(Wkx0)
OmX=OmX+C1*np.sin(thet)*np.sin(thet)
OmY=OmY+C1*np.cos(thet)*np.cos(thet)
OmZ=OmZ-2.*Edip/XM*ALPRe*(WK-1/ZR)**2*np.cos(Wkx0)

omega_j = np.array([np.sqrt(OmX), np.sqrt(OmY), np.sqrt(OmZ)])
print('mech freq Omx/2pi=',omega_j[0]/2/np.pi)
print('mech freq Omy/2pi=',omega_j[1]/2/np.pi)
print('mech freq Omz/2pi=',omega_j[2]/2/np.pi)

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    
### resulting variables ###
n_opt = 0
n_mech = k*T/(hbar * omega_j)
    
    
    
def calculate_couplings():
# Optomechanical couplings
    XZPF = np.sqrt(hbar/(2*XM*omega_j[0]))
    YZPF = np.sqrt(hbar/(2*XM*omega_j[1]))
    ZZPF = np.sqrt(hbar/(2*XM*omega_j[2]))

    GX = Ediph*WK*XZPF*np.sin(thet)*np.sin(Wkx0)
    GY = Ediph*WK*YZPF*np.cos(thet)*np.sin(Wkx0)
    GZ = -Ediph*(WK-1/ZR)*ZZPF*np.cos(Wkx0)
    # write(iwrite,*)'GX , GY, GZ in Hz',GX/2/pi,GY/2/pi,Gz/2/pi
# corrected sign on 29/8/2019
    GXY = Ediph*WK*XZPF*WK*YZPF*ALPRe*np.sin(2*thet)*np.cos(Wkx0)
    GZX = 2*Ediph*(WK-1/ZR)*ZZPF*WK*XZPF*ALPim*np.sin(Wkx0)*np.sin(thet)
    GYZ = 2*Ediph*(WK-1/ZR)*ZZPF*WK*YZPF*ALPim*np.sin(Wkx0)*np.cos(thet)
  #write(iwrite,*)'GXY , GYZ, GZX in Hz',GXY/2/pi,GYZ/2/pi,GZX/2/pi
  
    print('couplings')
    print('GX, GY, GZ', GX/2/np.pi, GY/2/np.pi, GZ/2/np.pi)
    print('GXY, GYZ, GZX', GXY/2/np.pi, GYZ/2/np.pi, GZX/2/np.pi)
    return [GX, GY, GZ, GXY, GYZ, GZX]

########################################################################
    
############    
### MAIN ###
############

# calculate coupling 
g = calculate_couplings()
'''

#param.Gamma = 2*param.Gamma 

param.prepare_calc()

param.print_param()
    
# optical damping rates
Gamma_opt = param.opt_damp_rate()

# photon numbers at equiv
N = tools.photon_number(param.n_mech, Gamma_opt, param.Gamma)


#param = (param.omega_mech, param.detuning, param.g, param.Gamma, param.kappa, param.n_opt, param.n_mech)
#print_parameter(param)

# actually the detuning is given as an angular freq
param.detuning = param.detuning *2*np.pi
    
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

plt.xlabel('$\omega/(2\pi)$ [kHz]')
plt.ylabel('$S_{ii}$ [a.u.]')
#plt.yscale('log')
plt.legend(loc = 'best')
plt.savefig('pic/spectrum_1D_calc')
plt.show()


# calculate photon numbers from area
Delta = np.abs(omega[1])-np.abs(omega[0])
NX_from_area = tools.n_from_area(SXX_plus, SXX_minus, Delta, N[0], 'X')
NY_from_area = tools.n_from_area(SYY_plus, SYY_minus, Delta, N[1], 'Y')
NZ_from_area = tools.n_from_area(SZZ_plus, SZZ_minus, Delta, N[2], 'Z')

