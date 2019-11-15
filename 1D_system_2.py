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

KAPP2 = kappa/2
Gamm2 = Gamma[0]/2
#################

#################
### FUNCTIONS ###
#################

def opt_damp_rate(_kappa, _detuning, _g, _omega_j):
    Gamma_opt = -_g**2 * _kappa * (1/(_kappa**2/4 + (_detuning+_omega_j)**2) - 1/(_kappa**2/4 + (_detuning-_omega_j)**2) )  

    print()
    print('optical damping rates')
    print('$\Gamma_{opt, x}$:', round(Gamma_opt[0]/1e5, 3), '1e5')
    print('$\Gamma_{opt, y}$:', round(Gamma_opt[1]/1e5, 3), '1e5')
    print('$\Gamma_{opt, z}$:', round(Gamma_opt[2]/1e5, 3), '1e5')
    
#    print('mechanical damping')
#    print('X', Gamma[0])
    
    return Gamma_opt


def photon_number(_n_j, _Gamma_opt, _Gamma_j):
    N = _n_j * _Gamma_j / (abs(_Gamma_opt) + 2*_Gamma_j)
    
    print()
    print('theoretical photon numbers at equiv')
    print('n_x theo: ', round(N[0], 1))
    print('n_y theo: ', round(N[1], 1))
    print('n_z theo: ', round(N[2], 1))
    
    print('theoretical photon numbers at room temperature')
    print('n_x theo: ', round(n_mech[0]/1e8, 3), '1e8')
    print('n_y theo: ', round(n_mech[1]/1e8, 3), '1e8')
    print('n_z theo: ', round(n_mech[2]/1e8, 3), '1e8')
    
    return N



#********************************************************************
#  Generic routine for noise spectra of trap and probe beams
#********************************************************************

def ANALYT(omega):

#************************************************************************
#************************************************************************
    
# *******WORK OUT NOISE SPECTRA
#or i in omega:
# First do susceptibilities
    #*****************************************
#Chi_R
    ANORM = KAPP2**2+ (omega+detuning)**2
    t1=KAPP2/ANORM
    t2= (omega+detuning)/ANORM
    CHIR1 = t1 + t2*1j
# chi_r^*(-omega)
    ANORM = KAPP2**2+ (-omega+detuning)**2
    t1=KAPP2/ANORM
    t2= (-omega+detuning)/ANORM
    CHISTMOM1=t1-t2*1j
    #******************************************
# X MECHANICAL susceptibilities
# chi_M X
    BNORM=(Gamm2)**2 + (omega-omega_j[0])**2
    t1= (Gamm2)/BNORM
    t2=(omega-omega_j[0])/BNORM
    CHIMX=t1+t2*1j
#CHI_M*(-om) X
    T1=(Gamm2)**2 + (-omega-omega_j[0])**2
    CHIMMOMX= Gamm2/T1 + 1j*(omega+omega_j[0])/T1


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



SXX_plus = np.zeros(len(omega))
SXX_minus = np.zeros(len(omega))

for i in range(len(omega)):
    SXX_plus[i] = ANALYT(omega[i])
    SXX_minus[i] = ANALYT(-omega[i])

plt.plot(omega/2/np.pi*1e-3, SXX_plus)
plt.plot(omega/2/np.pi*1e-3, SXX_minus)
plt.show()

#delta_omega = omega[1]-omega[0]

# calculate photon numbers from area
N_X_plus = area(SXX_plus, omega)
N_X_minus = area(SXX_minus, omega)

print()
print('photon number from area, difference')
print('X +:', round(N_X_plus, 2), '(', round((N_X_plus-N[0])/(N[0]+1)*100, 2), '% )')
print('X -:', round(N_X_minus, 2), '(', round((N_X_minus-N[0])/(N[0])*100, 2), '% )')
#print('area -', area(SXX_minus, omega))

'''
AVRE=(omx+Det2pi)**2+kapp2**2
AVRE=AVRE/4/omx/(-det2pi)
    
write(6,*)'X: back action limited phonons',AVRE
AVRE=(omy+Det2pi)**2+kapp2**2
AVRE=AVRE/4/omY/(-det2pi)
write(6,*)'Y: back action limited phonons',AVRE
AVRE=(omz+Det2pi)**2+kapp2**2
AVRE=AVRE/4/omz/(-det2pi)
write(6,*)'Z: back action limited phonons',AVRE
! integrate and normalise the quantum  noise spectra. Get temperature and quanta

      CALL NORM1(NPTS,OMX,TBATH,GAMMAM,TEMPX,SXXQM,OMSTOR,AVRE)
! Area AVRE corresponds to 2n+1 so convert to get n
                PHONONS=PHON(1)
     AVRE=0.5d0*(AVRE-1.d0)
     AV(1)=AVRE
      write(6,*)'X phonons from formula,  from SXX FT'
      write(6,200)PHONONS,AVRE
       CALL NORM1(NPTS,OMY,TBATH,GAMMAM,TEMPY,SYYQM,OMSTOR,AVRE)
! Area AVRE corresponds to 2n+1 so convert to get n
        PHONONS=PHON(2)
        AVRE=0.5d0*(AVRE-1.d0)
        AV(2)=AVRE
!write(6,*)'Y phonons from formula,  from SXX FT'
!        write(6,200)PHONONS,AVRE
write(6,*)'Y phonons from formula'
write(6,200)PHONONS

      CALL NORM1(NPTS,OMZ,TBATH,GAMMAM,TEMPZ,SZZQM,OMSTOR,AVRE)
! Area AVRE corresponds to 2n+1 so convert to get n
         PHONONS=PHON(3)
       AVRE=0.5d0*(AVRE-1.d0)
         AV(3)=AVRE
 !  write(6,*)'Z phonons from formula,  from SXX FT'
 !    write(6,200)PHONONS,AVRE
 write(6,*)'Z phonons from formula'
    write(6,200)PHONONS
!     write(14,120)thet,PHON(1),AV(1),PHON(2),AV(2),PHON(3),AV(3)
! loop over theta
10    enddo

!***********************************************************************
100   FORMAT(I3,3E14.6,1x,2(E14.6))
200  FORMAT(7E14.6)

    STOP
    END
!********************************************************************
!  Function to evaluate optomechanical cooling formula
!********************************************************************

     FUNCTION optocool(G1,Det1X,KAPP2,OMEGAM,GAMMAM)
       Implicit None
      double precision::G1,Det1X,KAPP2,OMEGAM

       double precision::C1,C2,C3,C4
       double precision::OPTOCOOL,COOL1,COOL2,GAMMAM
       double precision:: hbar,BOLTZ,TBATH
       PARAMETER(BOLTZ=1.4d-23,hbar=1.05d-34,TBATH=300.0)
        COOL1=0.d0


! now work out opto cooling expression 
! Trap beam
       C1=2.*KAPP2*G1*G1
       C2=DET1X


      C3=(C2+omegam)**2+kapp2**2
       C3=1.d0/C3
       C4=(C2-omegam)**2+kapp2**2
       C4=1.d0/C4
       COOL1=-C1*(C3-C4)

        optocool=COOL1
100   format(4D16.6)
   return
    end








'''

