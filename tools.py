import numpy as np
import inspect
import time
import datetime


'''
comments:
index i: function of omega
index j: different operators

'''


k = 1.380649e-23 #J/K
hbar = 1.054571817e-34 #Js
hbar = 1.05e-34 #Js
c = 3e8 #m/s
grav = 9.8 #m/s^2
Epsi0=8.854e-12 # vacuum permitivity [F m^-1]
 


### operators (only 1D by now)
a_in = np.zeros(8)
a_in_d = np.zeros(8)
b1_in = np.zeros(8)
b1_in_d = np.zeros(8)
b2_in = np.zeros(8)
b2_in_d = np.zeros(8)
b3_in = np.zeros(8)
b3_in_d = np.zeros(8)

a_in[0] = 1
a_in_d[1] = 1
b1_in[2] = 1
b1_in_d[3] = 1
b2_in[4] = 1
b2_in_d[5] = 1
b3_in[6] = 1
b3_in_d[7] = 1

### SUSCEPTIBILITIES ###
# chi
def chi(_omega, _omega_j, _Gamma):
    return 1/(-1j*(_omega-_omega_j) + _Gamma/2.0)

# optical
def eta(_omega, _detuning, _phi, _kappa):
    return np.exp(-1j*_phi)*chi(_omega, -_detuning, _kappa) - np.exp(1j*_phi)*np.conj(chi(-_omega, -_detuning, _kappa))
                  
# mechanical 
def mu(_omega, _omega_j, _Gamma):
    mu1 = chi(_omega, _omega_j[0], _Gamma) - np.conj(chi(-_omega, _omega_j[0], _Gamma))
    mu2 = chi(_omega, _omega_j[1], _Gamma) - np.conj(chi(-_omega, _omega_j[1], _Gamma))
    mu3 = chi(_omega, _omega_j[2], _Gamma) - np.conj(chi(-_omega, _omega_j[2], _Gamma))
    return [mu1, mu2, mu3]
    
### NOISES ###
# optical
def Q_opt(_omega, _detuning, _kappa, _phi):
    return np.exp(-1j*_phi)*np.einsum('i, j->ji', chi(_omega, -_detuning, _kappa), a_in) + np.exp(1j*_phi)*np.einsum('i, j->ji', np.conj(chi(-_omega, -_detuning, _kappa)), a_in_d)

# mechanical
def Q_mech(_omega, _omega_j, _Gamma):
    Q1 = np.einsum('i, j-> ji', chi(_omega, _omega_j[0], _Gamma),b1_in) + np.einsum('i,j -> ji', np.conj(chi(-_omega, _omega_j[0], _Gamma)),b1_in_d)
    Q2 = np.einsum('i, j-> ji', chi(_omega, _omega_j[1], _Gamma),b2_in) + np.einsum('i,j -> ji', np.conj(chi(-_omega, _omega_j[1], _Gamma)),b2_in_d)
    Q3 = np.einsum('i, j-> ji', chi(_omega, _omega_j[2], _Gamma),b3_in) + np.einsum('i,j -> ji', np.conj(chi(-_omega, _omega_j[2], _Gamma)),b3_in_d)
    return [Q1, Q2, Q3]
    #return [Q1]

### normalization factor
def M(_omega, _omega_j, _detuning, _phi, _Gamma, _kappa, _g):
    M1 = 1+ _g[0]**2 *mu(_omega, _omega_j, _Gamma)[0]*eta(_omega, _detuning, _phi[0], _kappa)
    M2 = 1+ _g[1]**2 *mu(_omega, _omega_j, _Gamma)[1]*eta(_omega, _detuning, _phi[1], _kappa)
    M3 = 1+ _g[2]**2 *mu(_omega, _omega_j, _Gamma)[2]*eta(_omega, _detuning, _phi[2], _kappa)
    return [M1, M2, M3]

### displacement operator
def q_1D(_omega, _omega_j, _detuning, _g, _Gamma, _kappa, _phi):
    _M = M(_omega, _omega_j, _detuning, _phi, _Gamma, _kappa, _g)
    _Q_mech = Q_mech(_omega, _omega_j, _Gamma)
    _mu = mu(_omega, _omega_j, _Gamma)
    
    q1 = np.sqrt(_Gamma)*np.einsum('i,ji -> ji',1/_M[0], _Q_mech[0]) + 1j*np.sqrt(_kappa)*_g[0]*np.einsum('i, i, ji->ji',1/_M[0], _mu[0], Q_opt(_omega, _detuning, _kappa, _phi[0]))
    q2 = np.sqrt(_Gamma)*np.einsum('i,ji -> ji',1/_M[1], _Q_mech[1]) + 1j*np.sqrt(_kappa)*_g[1]*np.einsum('i, i, ji->ji',1/_M[1], _mu[1], Q_opt(_omega, _detuning, _kappa, _phi[1]))
    q3 = np.sqrt(_Gamma)*np.einsum('i,ji -> ji',1/_M[2], _Q_mech[2]) + 1j*np.sqrt(_kappa)*_g[2]*np.einsum('i, i, ji->ji',1/_M[2], _mu[2], Q_opt(_omega, _detuning, _kappa, _phi[2]))

    return [q1, q2, q3]

def q_3D(_omega, _omega_j, _detuning, _g, _Gamma, _kappa, _phi):
    _M = M(_omega, _omega_j, _detuning, _phi, _Gamma, _kappa, _g)
    #_Q_mech = Q_mech(_omega, _omega_j, _Gamma)
    _mu = mu(_omega, _omega_j, _Gamma)
    q = q_1D(_omega, _omega_j, _detuning, _g, _Gamma, _kappa, _phi)
    
    # coupling terms
    GXY = 1j*eta(_omega, _detuning, 0, _kappa)*_g[0]*_g[1] + _g[3]
    GYX = GXY
    GYZ = -1j*eta(_omega, _detuning, np.pi/2, _kappa)*_g[1]*_g[2] + _g[4]
    #GYZ = 1j*eta(_omega, _detuning, -np.pi/2, _kappa)*_g[1]*_g[2] + _g[4]
    GZY = 1j*eta(_omega, _detuning, np.pi/2, _kappa)*_g[2]*_g[1] + _g[4]
    GXZ = -1j*eta(_omega, _detuning, np.pi/2, _kappa)*_g[0]*_g[2] + _g[5]
    #GXZ = 1j*eta(_omega, _detuning, -np.pi/2, _kappa)*_g[0]*_g[2] + _g[5]
    GZX = 1j*eta(_omega, _detuning, np.pi/2, _kappa)*_g[2]*_g[0] + _g[5]
    
    
    RXY = 1j * _mu[0] * GXY / _M[0]
    RYX = 1j * _mu[1] * GYX / _M[1]
    RXZ = 1j * _mu[0] * GXZ / _M[0]
    RZX = 1j * _mu[2] * GZX / _M[2]
    RYZ = 1j * _mu[1] * GYZ / _M[1]
    RZY = 1j * _mu[2] * GZY / _M[2]
    
    # set coupling to z mode to 0
    #RXZ = 0
    #RZX = 0
    #RYZ = 0
    #RZY = 0
    
    
    q1 = q[0] + RXY * q[1] + RXZ * q[2]
    q2 = q[1] + RYX * q[0] + RYZ * q[2]
    q3 = q[2] + RZX * q[0] + RZY * q[1]
    
    return [q1, q2, q3]


    
### helper
def expectation_value(_operator, _n, _pair):
    return (_n+1)* np.abs(_operator[2*_pair])**2 + _n * np.abs(_operator[2*_pair+1])**2

def spectrum(_operator, _n_opt, _n_mech):
    s_a = expectation_value(_operator, _n_opt, 0)
    s_b1 = expectation_value(_operator, _n_mech[0], 1)
    s_b2 = expectation_value(_operator, _n_mech[1], 2)
    s_b3 = expectation_value(_operator, _n_mech[2], 3)
    return s_a + s_b1 + s_b2 + s_b3


def spectrum_output(omega, _i, param, ThreeD):
    # define phases
    _phi = np.array([0,0,np.pi/2])
    n_opt = 0
    
    # calculate q operator
    if ThreeD == False:
        operator = q_1D(omega, param.omega_mech, param.detuning, param.g, param.Gamma, param.kappa, _phi)[_i]
    if ThreeD == True:
        operator = q_3D(omega, param.omega_mech, param.detuning, param.g, param.Gamma, param.kappa, _phi)[_i]
    
    # calculate spectrum
    _spectrum = spectrum(operator, n_opt, param.n_mech)
    
    return _spectrum




def area(_S, _Delta):     

#  integrate the position spectrum of bead
# quick hack - use trapezoidal rule- improve later
    summe=0
    #Delta=np.abs(omega[1]-omega[0])
    for i in range(len(_S)-1):
        Tem = 0.5*(_S[i]+_S[i+1])
        summe = summe + Tem
    return summe*_Delta / (2*np.pi) #/2pi for normalization
   
    

def n_from_area(_S_plus, _S_minus, _Delta_omega, _N = 0, _name = '', printing = True):
    N_X_plus = area(_S_plus, _Delta_omega)
    N_X_minus = area(_S_minus, _Delta_omega)
    N_X = (N_X_plus + N_X_minus -1)/2
    
    if printing == True:
        print()
        print(_name, 'photon number from area, difference')
        print('+:', round(N_X_plus, 2), '(', round((N_X_plus-_N)/(_N+1)*100, 2), '% )')
        print('-:', round(N_X_minus, 2), '(', round((N_X_minus-_N)/(_N)*100, 2), '% )')
        print('summed:', round(N_X, 2), '(', round((N_X-_N)/(_N)*100, 2), '% )')
    #print('area -', area(SXX_minus, omega))
    return [N_X_plus, N_X_minus, N_X]



def photon_number(_n_j, _Gamma_opt, _Gamma):
    N = 2*_n_j * _Gamma / (abs(_Gamma_opt) + 2*_Gamma)
    
    print()
    print('theoretical photon numbers at equiv')
    print('n_x theo: ', round(N[0], 4))
    print('n_y theo: ', round(N[1], 4))
    print('n_z theo: ', round(N[2], 4))
    
    print('theoretical photon numbers at room temperature')
    print('n_x theo: ', round(_n_j[0]/1e8, 4), '1e8')
    print('n_y theo: ', round(_n_j[1]/1e8, 4), '1e8')
    print('n_z theo: ', round(_n_j[2]/1e8, 4), '1e8')
    
    return N













class parameters:
    "This class contains all relevant parameters"
    
    # class variables
    #T = 300 #temperatur [K]
    #R0=100e-9 #sphere radius
    #RHO=2198 #sphere density
    #EPSR=2.1 #parameter used to obtain the refractive index
    #              rho=2198.d0,EPSR=1.45d0**2,Epsi0=8.854d-12, &
    #lambda_tw = 1064e-9 # wavelength [m]   
    #WK=5.9e6 #=2*pi/lambda=k
    #waist=41.1e-6 #waist radius
    #WX=0.67e-6
    #WY=0.77e-6
    #XL=1.07e-2 #cavity length 
    #Finesse=15e4
    #Press=1e-6 #air pressure in millibars
    #Pin1=0.17e0 #input power in Watts tweezer beam
    #detuning=-300e3 #detuning in KHz trap beam
    #DelFSR=14e9 #1 FSR= 14 GHz, free spectral range =separation between cavity modes, you don't use these if you input g_x
    #theta0=0.2 #angle between tweezer polarization and cavity axis. Given as as FRACTION of pi so pi/4  is theta0=0.25
    
    # Equilibrium positions 
    # X0=0.125*lambda ,Y0=waist/sqrt(2)
    #Y0=0.0e-6
    #X0=0.125*1.064e-6
    #X0=0.73*1.064e-6
    #Z0=0e-6
    
    n_opt = 0
        
    
    def print_param(self):
        print('\n *** PARAMETERS ***')
        print('T [K]: ', self.T)
        print('R0 [m]: ', self.R0)
        print('RHO [kg/m^3]: ', self.RHO)
        print('EPSR: ', self.EPSR)
        print('lambda_tw [m]: ', self.lambda_tw)
        print('waist [m]: ', self.waist)
        print('W_x [m]: ', self.WX)
        print('W_y [m]: ', self.WY)
        print('L_cav [m]: ', self.XL)
        print('Finesse: ', self.Finesse)
        print('P [mbar]: ', self.Press)
        print('P_in [W]: ', self.Pin1)
        print('detuning [kHz]: ', self.detuning*1e-3)
        print('Delta_FSR [Hz]: ', self.DelFSR)
        print('theta [pi]: ', self.theta0)
        print('X_0 [lambda_tw]: ', self.X0/self.lambda_tw)
        print('Y_0 [lambda_tw]: ', self.Y0/self.lambda_tw)
        print('Z_0 [lambda_tw]: ', self.Z0/self.lambda_tw)
        print('mass [kg]: ', self.XM)
        print('Polarisibility :', self.Polaris)
        print('epsilon_tw: ', self.epsTW)
        print('epsilon_c: ', self.epsCAV)
        print('kappa/2pi [Hz]: ', self.kappa /(2*np.pi))
        print('Gamma: ', self.Gamma)
        print('omega_mech/2pi [kHz]: ', self.omega_mech/(2*np.pi)*1e-3)
        print('photons in cavity: ', self.n_photon)
        print('GX, GY, GZ', self.g[0]/2/np.pi, self.g[1]/2/np.pi, self.g[2]/2/np.pi)
        print('GXY, GYZ, GZX', self.g[3]/2/np.pi, self.g[4]/2/np.pi, self.g[5]/2/np.pi)
        print('***************')
        
    def prepare_calc(self):
        # mass
        self.XM=self.RHO*4.*np.pi/3.*self.R0**3
        
        # polarisibility
        self.Polaris=4.*np.pi*Epsi0*(self.EPSR-1)/(self.EPSR+2)*self.R0**3 # from Tania
        
        # tweezer properties
        WK= 2*np.pi / self.lambda_tw#=2*pi/lambda=k
        OMOPT=c*WK
        W2=self.waist**2
        _epsTW=4*self.Pin1/(self.WX*self.WY*np.pi*c*Epsi0)
        self.epsTW = np.sqrt(_epsTW)
        
        # cavity properties
        VOL=self.XL*np.pi*W2/4 # add a factor of 4 here.
        KAPPin=np.pi*c/self.Finesse/self.XL
        _epsCAV=hbar*OMOPT/(2.*VOL*Epsi0)
        self.epsCAV = np.sqrt(_epsCAV)
        ZR=self.WX*self.WY*WK/2
        
        # linewiddth
        coeff=WK*self.Polaris/Epsi0/OMOPT**2/np.pi
        kappnano=4*coeff**2*self.DelFSR*np.cos(WK*self.X0)*np.cos(WK*self.X0)
        self.kappa=kappnano+KAPPin
         
        # Pressure 1.d-4 mBar => ~ 0.125Hz in old expt
        # now take usual expression eg Levitated review by Li Geraci etc
        # 1 bar= 10^ 5 pascal; Press is in mbar = 10^ 2 pascal
        #gamma=16 * P/(np.pi*v*RHO*R)
        # v=speed of air=500 /s
        GAMMAM=1600*self.Press/np.pi
        GAMMAM=GAMMAM/500/self.RHO/self.R0
        #Fix of Feb.2016 our GAMMAM => GAMMAM/2!!
        self.Gamma=GAMMAM/2
        
        # mechanical frequencies
        Det2pi=self.detuning*2*np.pi
        kapp2=0.5*self.kappa
        Wkx0=WK*self.X0 #was commented out
        OmX=self.Polaris*self.epsTW**2/self.XM/self.WX**2
        OmY=self.Polaris*self.epsTW**2/self.XM/self.WY**2
        OmZ=0.5*self.Polaris*self.epsTW**2/self.XM/ZR**2
        
        # theta vs thet
        thet = self.theta0 * np.pi
        
        # photon field
        Edip=-0.5*self.Polaris*self.epsTW*self.epsCAV*np.sin(thet)
        Ediph=Edip/hbar
        ALPRe=Det2pi*Ediph*np.cos(Wkx0)/(kapp2**2+Det2pi**2)
        ALPim=-kapp2*Ediph*np.cos(Wkx0)/(kapp2**2+Det2pi**2)
        Nphoton=Ediph*Ediph*np.cos(Wkx0)*np.cos(Wkx0)
        self.n_photon=Nphoton/(kapp2**2+Det2pi**2)
        
        # corrections to frequencies due to cavity
        C1=-Edip/self.XM*2.*ALPRe*WK**2*np.cos(Wkx0)
        OmX=OmX+C1*np.sin(thet)*np.sin(thet)
        OmY=OmY+C1*np.cos(thet)*np.cos(thet)
        OmZ=OmZ-2.*Edip/self.XM*ALPRe*(WK-1/ZR)**2*np.cos(Wkx0)
        self.omega_mech = np.array([np.sqrt(OmX), np.sqrt(OmY), np.sqrt(OmZ)])
        
        # phonon number at equilibrium
        self.n_mech = k*self.T/(hbar * self.omega_mech)
            
        ### COUPLINGS
        # Optomechanical couplings
        XZPF = np.sqrt(hbar/(2*self.XM*self.omega_mech[0]))
        YZPF = np.sqrt(hbar/(2*self.XM*self.omega_mech[1]))
        ZZPF = np.sqrt(hbar/(2*self.XM*self.omega_mech[2]))
    
        # light-matter couplings
        GX = Ediph*WK*XZPF*np.sin(thet)*np.sin(Wkx0)
        GY = Ediph*WK*YZPF*np.cos(thet)*np.sin(Wkx0)
        GZ = -Ediph*(WK-1/ZR)*ZZPF*np.cos(Wkx0)
        
        # matter-matter couplings
        GXY = Ediph*WK*XZPF*WK*YZPF*ALPRe*np.sin(2*thet)*np.cos(Wkx0)
        GZX = 2*Ediph*(WK-1/ZR)*ZZPF*WK*XZPF*ALPim*np.sin(Wkx0)*np.sin(thet)
        GYZ = 2*Ediph*(WK-1/ZR)*ZZPF*WK*YZPF*ALPim*np.sin(Wkx0)*np.cos(thet)
      
        self.g = np.array([GX, GY, GZ, GXY, GYZ, GZX])
     
        
    def opt_damp_rate(self):
        Det2pi = self.detuning * 2*np.pi
        g_1D = np.array([self.g[0], self.g[1], self.g[2]])
        Gamma_opt = -g_1D**2 * self.kappa * (1/(self.kappa**2/4 + (Det2pi+self.omega_mech)**2) - 1/(self.kappa**2/4 + (Det2pi-self.omega_mech)**2) )  
    
        print()
        print('optical cooling rates')
        print('$\Gamma_{opt, x}$:', round(Gamma_opt[0]/1e5, 3), '1e5')
        print('$\Gamma_{opt, y}$:', round(Gamma_opt[1]/1e5, 3), '1e5')
        print('$\Gamma_{opt, z}$:', round(Gamma_opt[2]/1e5, 3), '1e5')
            
        return Gamma_opt

def loop_progress(L_inner, L_outer, inner, outer, start_time):
    length = L_inner * L_outer
    progress = (inner + L_inner*outer)/length
    current = time.time()
    diff = current - start_time
    if progress != 0:
        rest_time = round((diff)/progress * (1-progress))
    else: 
        rest_time = 0
        
    current_time_form = str(datetime.timedelta(seconds=round(diff)))
    remaining_time_form = str(datetime.timedelta(seconds=rest_time))
    #print('\n completed: ',  round(progress*100, 2), '%, remaining time: ', str(datetime.timedelta(seconds=rest_time)) )
    print('completed: {0:.3f}%, running: {1:6}s, remaining: {2:6}s \r'.format(progress, current_time_form, remaining_time_form), end = '\r')
    
    
    