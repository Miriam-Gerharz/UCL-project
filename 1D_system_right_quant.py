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


omega = np.arange(0, 400, 1e-2)*1e3*2*np.pi # freq for spectrum

T = 300 #temperatur [K]
NPERIOD=160000 #NPERIOD=NOT RELEVANT TO ANALYTICAL CODE: ignore
NTOT=8 #NTOT=number of equations so 8=>1 optical mode +3D
R0=100e-9 #sphere radius
RHO=2198 #sphere density
EPSR=2.1
Epsi0=8.854e-12
#              rho=2198.d0,EPSR=1.45d0**2,Epsi0=8.854d-12, &
WK=5.9e6 #=2*pi/lambda=k
waist=41.1e-6 #waist radius
WX=0.67e-6
WY=0.77e-6
XL=1.07e-2 #cavity length 
Finesse=15e4
Press=1e-5 #air pressure in millibars
Pin1=0.17e0 #input power in Watts tweezer beam
detuning=-300e3 #detuning in KHz trap beam
DelFSR=14e9 #1 FSR= 14 GHz
theta0=0.2 #angle between tweezer polarization and cavity axis. Given as as FRACTION of pi so pi/4  is theta0=0.25

# Equilibrium positions 
# X0=0.125*lambda ,Y0=waist/sqrt(2)
Y0=0.0e-6
X0=0.125*1.064e-6
Z0=0e-6



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

omega_j = np.array([np.sqrt(OmX), np.sqrt(OmY), np.sqrt(OmZ)])
print('mech freq Omx/2pi=',omega_j[0]/2/np.pi)
print('mech freq Omy/2pi=',omega_j[1]/2/np.pi)
print('mech freq Omz/2pi=',omega_j[2]/2/np.pi)

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    
### resulting variables ###
n_opt = 0
n_mech = k*T/(hbar * omega_j)
   
NPTS = 2000/2
omega = np.linspace(0, 2*omega_j[0]+1, NPTS)
 
    
    
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

#################
### FUNCTIONS ###
#################
    
g = np.array(calculate_couplings()[0:3] )
Gamm2 = Gamma/2

#print(g)

#################
### FUNCTIONS ###
#################

def opt_damp_rate(_kappa, _detuning, _g, _omega_j):
    det2 = _detuning * 2*np.pi
    Gamma_opt = -_g**2 * _kappa * (1/(_kappa**2/4 + (det2+_omega_j)**2) - 1/(_kappa**2/4 + (det2-_omega_j)**2) )  

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
    print('n_x theo: ', round(n_mech[0]/1e8, 5), '1e8')
    print('n_y theo: ', round(n_mech[1]/1e8, 5), '1e8')
    print('n_z theo: ', round(n_mech[2]/1e8, 5), '1e8')
    
    return N

CHIR1 = 0
CHISTMOM1 = 0
CHIMX = 0
CHIMMOMX = 0
CHIMY = 0
CHIMMOMY = 0
CHIMZ = 0
CHIMMOMZ = 0

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
    
def ANALYT(GMAT,AVNX,AVNY,AVNZ,AVPHOT,THETA,DET,Kapp2,GAMMAM,OMX,OMY,OMZ,omega,XXF,YYF,ZZF,SHOM1):
    detuning = DET

#************************************************************************
#************************************************************************
# First do susceptibilities
##### ROUTINE SUCEPT #########
    Gamm2 = GAMMAM
    ANORM = KAPP2**2+ (omega+detuning)**2
    t1=KAPP2/ANORM
    t2= (omega+detuning)/ANORM
    CHIR1 = t1 + t2*1j
    
    #print(CHIR1)
# chi_r^*(-omega)
    ANORM = KAPP2**2+ (-omega+detuning)**2
    t1=KAPP2/ANORM
    t2= (-omega+detuning)/ANORM
    CHISTMOM1=t1-t2*1j
    #******************************************
# X MECHANICAL susceptibilities
# chi_M X
    BNORM=(Gamm2)**2 + (omega-OMX)**2
    t1= (Gamm2)/BNORM
    t2=(omega-OMX)/BNORM
    CHIMX=t1+t2*1j
#CHI_M*(-om) X
    T1=(Gamm2)**2 + (-omega-OMX)**2
    CHIMMOMX= Gamm2/T1 + 1j*(omega+OMX)/T1
    
    
# Y MECHANICAL susceptibilities
# chi_M Y
    BNORM=(Gamm2)**2 + (omega-OMY)**2
    t1= (Gamm2)/BNORM
    t2=(omega-OMY)/BNORM
    CHIMY=t1 + 1j*t2
    # CHI_M*(-om) Y
    T1=(Gamm2)**2 + (-omega-OMY)**2
    CHIMMOMY=Gamm2/T1 + 1j* (omega+OMY)/T1
    #****************************************************
    #******************************************
    # Z MECHANICAL susceptibilities
    # chi_M Z
    BNORM=(Gamm2)**2 + (omega-OMZ)**2
    t1= (Gamm2)/BNORM
    t2=(omega-OMZ)/BNORM
    CHIMZ=t1 + 1j*t2
    # CHI_M*(-om) Z
    T1=(Gamm2)**2 + (-omega-OMZ)**2
    CHIMMOMZ=Gamm2/T1 + 1j* (omega+OMZ)/T1
    #****************************************************
########################################


#########################################
### ROUTINE AVECT ####
    Gamm = GAMMAM
    GX=GMAT[0]
    GY=GMAT[1]
    GZ=GMAT[2]
    GXY=0
    GYZ=0
    GZX=0
#&&&&&&&&& 2D  &&&&&&&&&&&&&&&&&&
# as we only want x spectrum, switch off other couplings
    GY=0
    GZ=0
    GXY=0
    GYZ=0
    GZX=0

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
# actually Gamm is gamma/2
    Gamm2=Gamm
# zero arrays
    BAX=np.array([0+0j,0,0,0,0,0,0,0])
    BAY=np.array([0+0j,0,0,0,0,0,0,0])
    BAZ=np.array([0+0j,0,0,0,0,0,0,0])

    NTOT = 8
    N0X=np.array([0+1j,0,0,0,0,0,0,0])
    N0Y=np.array([0+1j,0,0,0,0,0,0,0])
    N0Z=np.array([0+1j,0,0,0,0,0,0,0])
    #A1=0
    #A1dagg=0
# i
    XI=1j
    ONE=1

#     SUSCEPTIBILITIES
    eta0c=CHIR1-CHISTMOM1
    etaMpi2c=XI*(CHIR1+CHISTMOM1)
    etaPpi2c=-XI*(CHIR1+CHISTMOM1)
    etax=CHIMX-CHIMMOMX
    etay=CHIMY-CHIMMOMY
    etaz=CHIMZ-CHIMMOMZ

# coeff of X-Y coupling- Combines direct and indirect paths
    GcoefXY=(GXY+1j*eta0c*GX*GY)
    GcoefYX=(GXY+1j*eta0c*GX*GY)
    GcoefXZ=(GZX+1j*etaMpi2c*GZ*GX)
    GcoefZX=(GZX+1j*etaPpi2c*GZ*GX)
    GcoefYZ=(GYZ+1j*etaMpi2c*GY*GZ)
    GcoefZY=(GYZ+1j*etaPpi2c*GY*GZ)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#        GcoefXY=CMPLX(0.d0,0.d0)
 #      GcoefZX=CMPLX(0.d0,0.d0)
 #    GcoefYZ=CMPLX(0.d0,0.d0)
# NORMALIZATIONS
    CMX=1+GX**2*etax*eta0c
    CMY=1+GY**2*etay*eta0c
    CMZ=1+GZ**2*etaz*eta0c
    

    Sqrtkapp=np.sqrt(2*KAPP2)
    Sqrtgamm=np.sqrt(2*Gamm2)
    BETX=1j*Sqrtkapp*etax*GX
    BETY=1j*Sqrtkapp*etay*GY
    BETZ=1j*Sqrtkapp*etaz*GZ
#&&&&&&&&&&&&&&&&&&&&&&&&
# zero-th order X noise vector; weights of a1,a1*,bx,bx*,by,by*,bz,bz*
   # print(CHIR1.size)
    N0X[0]=BETX*CHIR1/CMX
    N0X[1]=BETX*CHISTMOM1/CMX
    N0X[2]=Sqrtgamm*CHIMX/CMX
    N0X[3]=Sqrtgamm*CHIMMOMX/CMX

# Zero-th order Y noise vector;weights of a1,a1*,bx,bx*,by,by*,bz,bz*
    N0Y[0]=BETY*CHIR1/CMY
    N0Y[1]=BETY*CHISTMOM1/CMY
    N0Y[4]=Sqrtgamm*CHIMY/CMY
    N0Y[5]=Sqrtgamm*CHIMMOMY/CMY
# Zero-th Z noise vector;weights of a1,a1*,bx,bx*,by,by*,bz,bz*
    N0Z[0]=-1j*BETZ*CHIR1/CMZ
    N0Z[1]=1j*BETZ*CHISTMOM1/CMZ
    N0Z[6]=Sqrtgamm*CHIMZ/CMZ
    N0Z[7]=Sqrtgamm*CHIMMOMZ/CMZ

    
#    for i in range(NTOT):
#       BAX[i]=N0X[i]
#       BAY[i]=N0Y[i]
#       BAZ[i]=N0Z[i]
    
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# Higher order
    RXY=1j*etax*GcoefXY/CMX
    RYX=1j*etay*GcoefYX/CMY
    RXZ=1j*etax*GcoefXZ/CMX
    RZX=1j*etaz*GcoefZX/CMZ
    RYZ=1j*etay*GcoefYZ/CMY
    RZY=1j*etaz*GcoefZY/CMZ

    #print(RXY, RYX, RXZ, RZX, RYZ, RZY)
#      RXY=CMPLX(0.d0,0.d0)
 #     RYX=CMPLX(0.d0,0.d0)
#     RXZ=CMPLX(0.d0,0.d0)
#     RZX=CMPLX(0.d0,0.d0)
#     RYZ=CMPLX(0.d0,0.d0)
#     RZY=CMPLX(0.d0,0.d0)

    CNORM=1-RZX*RXZ-RZY*RYZ-RYX*RXY-RZX*RXY*RYZ-RYX*RXZ*RZY
    
# ADD 3D BACK-ACTION TERMS
    for i in range(NTOT):
        CSUM=(1-RZY*RYZ)*N0X[i]+(RXY+RXZ*RZY)*N0Y[i]+(RXZ+RXY*RYZ)*N0Z[i]
        BAX[i]=BAX[i]+CSUM/CNORM
        CSUM=(1-RZX*RXZ)*N0Y[i]+(RYX+RYZ*RZX)*N0X[i]+(RYZ+RYX*RXZ)*N0Z[i]
        BAY[i]=BAY[i]+CSUM/CNORM
        CSUM=(1-RYX*RXY)*N0Z[i]+(RZX+RZY*RYX)*N0X[i]+(RZY+RZX*RXY)*N0Y[i]
        BAZ[i]=BAZ[i]+CSUM/CNORM
    
    #print(omega, BAX[0])
    
# now work out the optical trap field =a1
    CA1=XI*CHIR1
# now work out the photon field a1dagger
    CA1dagg=-XI*CHISTMOM1
    
    A1 = np.array([0+1j,0,0,0,0,0,0,0])
    A1dagg = np.array([0+1j,0,0,0,0,0,0,0])
#
    for i in range(NTOT):
        A1[i]=CA1*(GX*BAX[i]+GY*BAY[i]+GZ*BAZ[i])
        A1dagg[i]=CA1dagg*(GX*BAX[i]+GY*BAY[i]+GZ*BAZ[i])
         
# add shot or incoming noise
# trap beam: add cavity-filtered contribution
    A1[1]=A1[1]+Sqrtkapp*CHIR1
    A1dagg[2]=A1dagg[2]+Sqrtkapp*CHISTMOM1

# cavity output : add incoming imprecision
# work out a_out=a_in-Sqrtkapp(a)
    for i in range(NTOT):
        A1[i]=-A1[i]*Sqrtkapp
        A1dagg[i]=-A1dagg[i]*Sqrtkapp
      
    A1[1]=ONE+A1[1]
    A1dagg[2]=ONE+A1dagg[2]
#####################################################################
    
#####################################################################
    ### ROUTINE HOMODYNE ###
    NTOT = 8
    
    theta = 0
    XAM1 = np.array([0+1j,0,0,0,0,0,0,0])
    XPM1 = np.array([0+1j,0,0,0,0,0,0,0])
    XTHET1 = np.array([0+1j,0,0,0,0,0,0,0]) 
    for i in range(8):
        XAM1[i]=A1[i]+A1dagg[i]
        XPM1[i]=1j*(A1[i]-A1dagg[i])
        
        XTHET1[i]=XPM1[i]*np.sin(theta)+XAM1[i]*np.cos(theta)
           
    SHOM1=0
    
    SHOM1=SHOM1+abs(XTHET1[0])**2*AVPHOT+abs(XTHET1[1])**2*(AVPHOT+1)
    SHOM1=SHOM1+abs(XTHET1[2])**2*AVNX+abs(XTHET1[3])**2*(AVNX+1)
    SHOM1=SHOM1+abs(XTHET1[4])**2*AVNY+abs(XTHET1[5])**2*(AVNY+1)
    SHOM1=SHOM1+abs(XTHET1[6])**2*AVNZ+abs(XTHET1[7])**2*(AVNZ+1)
    #SHOM1 = tools.spectrum(XTHET1, n_opt, n_mech)

###########################################################################    
    #      G2NORM=G2*G2*(ABS(CHIR2-cos(2*theta)*CHISTMOM2))**2
    
    # work out noise vector for x, a1 and a2
    #Avect(GMAT,GAMMAM,Kapp2,CHIR1,CHISTMOM1,CHIMX,CHIMMOMX,CHIMY,CHIMMOMY,CHIMZ,CHIMMOMZ,A1,A1dagg,BAX,BAY,BAZ)
    # work out homodyne spectra
    
    XXF=(abs(BAX[0]))**2
    XXF=XXF+ (AVNX+1)*(abs(BAX[2]))**2+AVNX*(abs(BAX[3]))**2+(AVNY+1)*abs(BAX[4])**2+AVNY*abs(BAX[5])**2
    XXF=XXF+ (AVNZ+1)*(abs(BAX[6]))**2+AVNZ*(abs(BAX[7]))**2

    YYF=(abs(BAY[0]))**2
    YYF=YYF+ (AVNX+1)*(abs(BAY[2]))**2+AVNX*(abs(BAY[3]))**2+(AVNY+1)*abs(BAY[4])**2+AVNY*abs(BAY[5])**2
    YYF=YYF+ (AVNZ+1)*(abs(BAY[6]))**2+AVNZ*(abs(BAY[7]))**2

    ZZF=(abs(BAZ[0]))**2
    ZZF=ZZF+ (AVNX+1)*(abs(BAZ[2]))**2+AVNX*(abs(BAZ[3]))**2+(AVNY+1)*abs(BAZ[4])**2+AVNY*abs(BAZ[5])**2
    ZZF=ZZF+ (AVNZ+1)*(abs(BAZ[6]))**2+AVNZ*(abs(BAZ[7]))**2


    #XXF = tools.spectrum(BAX, n_opt, n_mech)
    #XXF = tools.spectrum(BAX, n_opt, n_mech)
    #XXF = tools.spectrum(BAX, n_opt, n_mech)
    return XXF, YYF, ZZF
    

'''
#      SUSCEPTIBILITIES
    eta0c=CHIR1-CHISTMOM1
    etaMpi2c=1j*(CHIR1+CHISTMOM1)
    etaPpi2c=-1j*(CHIR1+CHISTMOM1)
    etaX=CHIMX-CHIMMOMX
    #etaY=CHIMY-CHIMMOMY
    #etaZ=CHIMZ-CHIMMOMZ

# NORMALIZATIONS
    CMX=1+g[0]**2*etaX*eta0c
    

    Sqrtkapp=np.sqrt(2*KAPP2)
    Sqrtgamm=np.sqrt(2*Gamm2)
    BETX=1j*Sqrtkapp*etaX*g[0]
    
#&&&&&&&&&&&&&&&&&&&&&&&&
# zero-th order X noise vector; weights of a1,a1*,bx,bx*,by,by*,bz,bz*
    N0X1=BETX*CHIR1/CMX
    N0X2=BETX*CHISTMOM1/CMX
    N0X3=Sqrtgamm*CHIMX/CMX
    N0X4=Sqrtgamm*CHIMMOMX/CMX


    XPM1=0
    XAM1=0

    XTHET1=0

      
#XX= sqrt(0.5) (b+b^dagg) so halve XX,
    XXF=np.abs(N0X1)**2
    XXF=XXF+ (n_mech[0]+1)*np.abs(N0X3)**2+n_mech[0]*np.abs(N0X4)**2
    return XXF
'''


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
Gamma_opt = opt_damp_rate(kappa, detuning, g, omega_j)

# photon numbers at equiv
N = photon_number(n_mech, Gamma_opt, Gamma)

AVNX = n_mech[0]
AVNY = n_mech[1]
AVNZ = n_mech[2]
AVPHOT = 0
THETAHOM = 0
DET2pi = detuning*2*np.pi
Kapp2 = kappa/2
GAMMAM = Gamma
OMX = omega_j[0]
OMY = omega_j[1]
OMZ = omega_j[2]
XXQM = 0
YYQM = 0 
ZZQM = 0
SHOM1 = 0
GMAT = g

SXX_plus = np.zeros(len(omega))
SXX_minus = np.zeros(len(omega))
for i in range(len(omega)):
    OMsweep = omega[i]
    XXF, YYF, ZFF = ANALYT(GMAT,AVNX,AVNY,AVNZ,AVPHOT,THETAHOM,DET2pi,Kapp2,GAMMAM,OMX,OMY,OMZ,-OMsweep,XXQM,YYQM,ZZQM,SHOM1)
    SXX_minus[i] = XXF 
    XXF, YYF, ZZF = ANALYT(GMAT,AVNX,AVNY,AVNZ,AVPHOT,THETAHOM,DET2pi,Kapp2,GAMMAM,OMX,OMY,OMZ,OMsweep,XXQM,YYQM,ZZQM,SHOM1)
    SXX_plus[i] = XXF

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