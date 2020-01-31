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
Epsi0=8.854e-12

###########################################################################
##################
### PARAMETERS ###
##################


omega = np.arange(0, 400, 1e-2)*1e3*2*np.pi # freq for spectrum


param = tools.parameters()

param.T = 300 #temperatur [K]
#param.NPERIOD=160000 #NPERIOD=NOT RELEVANT TO ANALYTICAL CODE: ignore
#NTOT=8 #NTOT=number of equations so 8=>1 optical mode +3D
param.R0=100e-9 #sphere radius
param.RHO=2198 #sphere density
param.EPSR=2.1
#              rho=2198.d0,EPSR=1.45d0**2,Epsi0=8.854d-12, &
param.lambda_tw = 1064e-9
param.WK=5.9e6 #=2*pi/lambda=k
param.waist=41.1e-6 #waist radius
param.WX=0.67e-6
param.WY=0.77e-6
param.XL=1.07e-2 #cavity length 
param.Finesse=15e4
param.Press=1e-5 #air pressure in millibars
param.Pin1=0.17e0 #input power in Watts tweezer beam
param.detuning=-300e3 #detuning in KHz trap beam
param.DelFSR=14e9 #1 FSR= 14 GHz
param.theta0=0.2 #angle between tweezer polarization and cavity axis. Given as as FRACTION of pi so pi/4  is theta0=0.25

# Equilibrium positions 
# X0=0.125*lambda ,Y0=waist/sqrt(2)
param.Y0=0.0e-6
param.X0=0.125*1.064e-6
param.Z0=0e-6



###############################################################



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
    
 
    
############################
### PREPARE CALCULATIONS ###
############################                

param.prepare_calc()

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
KAPP2 = kapp2
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

omega_mech = np.array([np.sqrt(OmX), np.sqrt(OmY), np.sqrt(OmZ)])
print('mech freq Omx/2pi=',omega_mech[0]/2/np.pi)
print('mech freq Omy/2pi=',omega_mech[1]/2/np.pi)
print('mech freq Omz/2pi=',omega_mech[2]/2/np.pi)

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    
### resulting variables ###
n_opt = 0
n_mech = k*T/(hbar * omega_mech)
  
'''
 
NPTS = 2000/2
omega = np.linspace(0, 2*param.omega_mech[0]+1, NPTS)
 
    
'''
def calculate_couplings():
# Optomechanical couplings
    XZPF = np.sqrt(hbar/(2*XM*omega_mech[0]))
    YZPF = np.sqrt(hbar/(2*XM*omega_mech[1]))
    ZZPF = np.sqrt(hbar/(2*XM*omega_mech[2]))

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
'''

########################################################################
    
############    
### MAIN ###
############

#################
### FUNCTIONS ###
#################
    
#g = np.array(calculate_couplings()[0:3] )
GMAT = param.g
#print(GMAT/2/np.pi)
Gamm2 = param.Gamma/2

#print(g)

#################
### FUNCTIONS ###
#################

def opt_damp_rate(_kappa, _detuning, _g, _omega_mech):
    det2 = _detuning * 2*np.pi
    Gamma_opt = -_g**2 * _kappa * (1/(_kappa**2/4 + (det2+_omega_mech)**2) - 1/(_kappa**2/4 + (det2-_omega_mech)**2) )  

    print()
    print('optical damping rates')
    print('$\Gamma_{opt, x}$:', round(Gamma_opt[0]/1e5, 5), '1e5')
    print('$\Gamma_{opt, y}$:', round(Gamma_opt[1]/1e5, 5), '1e5')
    print('$\Gamma_{opt, z}$:', round(Gamma_opt[2]/1e5, 5), '1e5')
    
#    print('mechanical damping')
#    print('X', Gamma[0])
    
    return Gamma_opt


def photon_number(_n_j, _Gamma_opt, _Gamma_j):
    N = 2*_n_j * _Gamma_j / (abs(_Gamma_opt) + 2*_Gamma_j)
    
    print()
    print('theoretical photon numbers at equiv')
    print('n_x theo: ', round(N[0], 3))
    print('n_y theo: ', round(N[1], 3))
    print('n_z theo: ', round(N[2], 3))
    
    print('theoretical photon numbers at room temperature')
    print('n_x theo: ', round(param.n_mech[0]/1e8, 5), '1e8')
    print('n_y theo: ', round(param.n_mech[1]/1e8, 5), '1e8')
    print('n_z theo: ', round(param.n_mech[2]/1e8, 5), '1e8')
    
    return N

#CHIR1 = 0
#CHISTMOM1 = 0
#CHIMX = 0
#CHIMMOMX = 0
#CHIMY = 0
#CHIMMOMY = 0
#CHIMZ = 0
#CHIMMOMZ = 0

#********************************************************************
#  Generic routine for noise spectra of trap and probe beams
#********************************************************************
#def Suscept(OMEGA,DET1x,Kapp2,gamm,OMX,OMY,OMZ,CHIR1,CHISTMOM1,CHIMX,CHIMMOMX,CHIMY,CHIMMOMY,CHIMZ,CHIMMOMZ):
    #Suscept(OMEGA,DET,Kapp2,GAMMAM,OMX,OMY,OMZ,CHIR1,CHISTMOM1,CHIMX,CHIMMOMX,CHIMY,CHIMMOMY,CHIMZ,CHIMMOMZ)
    
    #detuning = det

#def Avect(GMAT,Gamm,kapp2,CHIR1,CHISTMOM1,CHIMX,CHIMMOMX,CHIMY,CHIMMOMY,CHIMZ,CHIMMOMZ,A1,A1dagg,BAX,BAY,BAZ):
#    Avect(GMAT,GAMMAM,Kapp2,CHIR1,CHISTMOM1,CHIMX,CHIMMOMX,CHIMY,CHIMMOMY,CHIMZ,CHIMMOMZ,A1,A1dagg,BAX,BAY,BAZ)
    



#def Homodyne(NTOT,AVNX,AVNY,AVNZ,AVPHOT,THETA,A1,A1dagg,SHOM1):
#    Homodyne(NT,AVNX,AVNY,AVNZ,AVPHOT,THETA,A1,A1dagg,SHOM1)
    
XXF = 0
    
def ANALYT(GMAT,AVNX,AVNY,AVNZ,AVPHOT,THETA,DET,Kapp2,GAMMAM,OMX,OMY,OMZ,omega):
    detuning = DET

#************************************************************************
#************************************************************************
# First do susceptibilities
##### ROUTINE SUCEPT #########
    kappa = 2*Kapp2
    #Gamm2 = GAMMAM # GAMMA is actually Gamma/2
   
    #****************************************************
########################################


#########################################
### ROUTINE AVECT ####
#    GMAT[1] = 0
#    GMAT[2] = 0
#    GMAT[3] = 0
#    GMAT[4] = 0
#    GMAT[5] = 0
    
    
    #Gamm = GAMMAM

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
# actually Gamm is gamma/2
    #Gamm2=Gamm
    
    Gamma = GAMMAM

    phi = np.array([0,0,+np.pi/2])
    omega_a = np.array([omega])
    BAX, BAY, BAZ = tools.q_1D(omega_a, param.omega_mech, detuning, GMAT, Gamma, kappa, phi)
    #BAX[0] = tools.q_1D(omega_a, omega_mech, detuning, GMAT, Gamma, kappa, phi)[0][0]
    #BAX[1] = tools.q_1D(omega_a, omega_mech, detuning, GMAT, Gamma, kappa, phi)[0][1]
    #BAX[3] = tools.q_1D(omega_a, omega_mech, detuning, GMAT, Gamma, kappa, phi)[0][3]
    
 

###########################################################################    
   

    XXF = tools.spectrum(BAX, param.n_opt, param.n_mech)
    YYF = tools.spectrum(BAY, param.n_opt, param.n_mech)
    ZZF = tools.spectrum(BAZ, param.n_opt, param.n_mech)
    return XXF, YYF, ZZF
 

def area(S, omega):     

#  integrate the position spectrum of bead
# quick hack - use trapezoidal rule- improve later
    summe=0
    Delta=np.abs(omega[1]-omega[0])
    for i in range(len(S)-1):
        Tem = 0.5*(S[i]+S[i+1])
        summe = summe + Tem
    return summe*Delta / (2*np.pi)
   



#***************************************************************
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# open loop over frequency for noise spectra eg S_xx(\omega)=FT(autocorrelation <X(t)X^T(t+tau)>)


# optical damping rates
#Gamma_opt = opt_damp_rate(param.kappa, param.detuning, param.g, param.omega_mech)
Gamma_opt = param.opt_damp_rate()

param.print_param()

# photon numbers at equiv
#N = photon_number(param.n_mech, Gamma_opt, param.Gamma)
N = tools.photon_number(param.n_mech, Gamma_opt, param.Gamma)

#AVNX = n_mech[0]
#AVNY = n_mech[1]
#AVNZ = n_mech[2]
AVPHOT = 0
THETAHOM = 0
DET2pi = param.detuning*2*np.pi
Kapp2 = param.kappa/2
#GAMMAM = Gamma
#OMX = omega_mech[0]
#OMY = omega_mech[1]
#OMZ = omega_mech[2]
#XXQM = 0
#YYQM = 0 
#ZZQM = 0
#SHOM1 = 0

SXX_plus = np.zeros(len(omega))
SXX_minus = np.zeros(len(omega))
for i in range(len(omega)):
    OMsweep = omega[i]
    XXF, YYF, ZFF = ANALYT(param.g,param.n_mech[0],param.n_mech[1],param.n_mech[2],AVPHOT,THETAHOM,DET2pi,Kapp2,param.Gamma,param.omega_mech[0],param.omega_mech[1],param.omega_mech[2],-OMsweep)
    #XXF = tools.spectrum_output(-omega[i], 0, param, False)
    SXX_minus[i] = XXF 
    XXF, YYF, ZZF = ANALYT(GMAT,param.n_mech[0],param.n_mech[1],param.n_mech[2],AVPHOT,THETAHOM,DET2pi,Kapp2,param.Gamma,param.omega_mech[0],param.omega_mech[1],param.omega_mech[2],OMsweep)
    SXX_plus[i] = XXF

#XXF, YYF, ZFF = ANALYT(GMAT,AVNX,AVNY,AVNZ,AVPHOT,THETAHOM,DET2pi,Kapp2,GAMMAM,OMX,OMY,OMZ,-omega,XXQM,YYQM,ZZQM,SHOM1)
#SXX_minus = XXF 
#XXF, YYF, ZZF = ANALYT(GMAT,AVNX,AVNY,AVNZ,AVPHOT,THETAHOM,DET2pi,Kapp2,GAMMAM,OMX,OMY,OMZ,omega,XXQM,YYQM,ZZQM,SHOM1)
#SXX_plus = XXF



plt.plot(omega/2/np.pi*1e-3, SXX_plus)
plt.plot(omega/2/np.pi*1e-3, SXX_minus)
plt.show()

#delta_omega = omega[1]-omega[0]

# calculate photon numbers from area
N_X_plus = area(SXX_plus, omega)
N_X_minus = area(SXX_minus, omega)
N_X = (N_X_plus + N_X_minus -1)/2

print()
print('photon number from area, difference')
print('X +:', round(N_X_plus, 4), '(', round((N_X_plus-N[0])/(N[0]+1)*100, 2), '% )')
print('X -:', round(N_X_minus, 4), '(', round((N_X_minus-N[0])/(N[0])*100, 2), '% )')
print('X summed:', round(N_X, 4), '(', round((N_X-N[0])/(N[0])*100, 2), '% )')
#print('area -', area(SXX_minus, omega))


#for i in range(len(omega)):
#    print(omega[i]/2/np.pi*1e-6, SXX_plus[i])