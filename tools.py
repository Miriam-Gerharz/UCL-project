import numpy as np
import time
import datetime


'''
comments:
index i: function of omega
index j: different operators

'''
#############################################
### CONSTANTS ###
k = 1.380649e-23 #J/K
k = 1.4e-23
hbar = 1.054571817e-34 #Js
hbar = 1.05e-34 #Js
c = 3e8 #m/s
grav = 9.8 #m/s^2
Epsi0=8.854e-12 # vacuum permittivity [F m^-1]
#############################################



### SUSCEPTIBILITIES ###
# chi
def chi(_omega, _omega_j, _Gamma):
    """Calculates what is defined as chi in the paper
    
    Parameters
    ----------
    _omega : 1D numpy array
        The frequency range of which chi shall be calculated
    _omega_j : float
        respective mechanical frequency
    _Gamma : float
        Damping (either $\Gamma$ or $\kappa$)
        
    Returns
    -------
     : np.array 
         chi(omega)
    """
    #print('Gamma, chi', _Gamma)
    return 1/(-1j*(_omega-_omega_j) + _Gamma/2.0)

# optical
def eta(_omega, _detuning, _phi, _kappa):
    """Calculates optical susceptibility
    
    Parameters
    ---------
    _omega : numpy array
        The frequency range of which chi shall be calculated
    _detuning : float
        The detuning
    _phi : float
        Phase ???
    _kappa : float
        cavity linewidth
        
    Returns
    -------
     : np.array
        eta(omega)
    """
    return np.exp(-1j*_phi)*chi(_omega, -_detuning, _kappa) - np.exp(1j*_phi)*np.conj(chi(-_omega, -_detuning, _kappa))
                  
# mechanical 
def mu(_omega, _omega_j, _Gamma):
    """Calculates the mechanical susceptibilities
    
    Parameters
    ----------
    _omega : 1D numpy array
        The frequency range of which chi shall be calculated
    _omega_j : numpy array of length 3
        mechanical frequencies
    _Gamma : float
        Damping (either $\Gamma$ or $\kappa$)
        
    Returns
    -------
     : np.array 
        mu(omega)
    """
    mu1 = chi(_omega, _omega_j[0], _Gamma) - np.conj(chi(-_omega, _omega_j[0], _Gamma))
    mu2 = chi(_omega, _omega_j[1], _Gamma) - np.conj(chi(-_omega, _omega_j[1], _Gamma))
    mu3 = chi(_omega, _omega_j[2], _Gamma) - np.conj(chi(-_omega, _omega_j[2], _Gamma))
    return [mu1, mu2, mu3]
    
### NOISES ###
# optical
def Q_opt(_omega, _detuning, _kappa, _phi):
    """Calculates optical noise
    
    Parameters
    ---------
    _omega : 1D numpy array
        The frequency range of which chi shall be calculated
    _detuning: float
        The detuning
    _kappa : float
        cavity linewidth
    _phi : float
        Phase ???
       
    Returns
    -------
     : 2D np.array
         Q_opt(omega, mode)
    """
    
    # define operators
    a_in = np.zeros(8)
    a_in_d = np.zeros(8)
    
    a_in[0] = 1
    a_in_d[1] = 1
    
    return np.exp(-1j*_phi)*np.einsum('i, j->ji', chi(_omega, -_detuning, _kappa), a_in) + np.exp(1j*_phi)*np.einsum('i, j->ji', np.conj(chi(-_omega, -_detuning, _kappa)), a_in_d)

# mechanical
def Q_mech(_omega, _omega_j, _Gamma):
    """Calculates the mechanical noises
    
    Parameters
    ----------
    _omega : 1D numpy array
        The frequency range of which chi shall be calculated
    _omega_j : numpy array of length 3
        mechanical frequencies
    _Gamma : float
        Damping (either $\Gamma$ or $\kappa$)
        
    Returns
    -------
     : 2D np.array 
         Q_mech(omega, mode) 
    """
    ### operators
    b1_in = np.zeros(8)
    b1_in_d = np.zeros(8)
    b2_in = np.zeros(8)
    b2_in_d = np.zeros(8)
    b3_in = np.zeros(8)
    b3_in_d = np.zeros(8)
    
    b1_in[2] = 1
    b1_in_d[3] = 1
    b2_in[4] = 1
    b2_in_d[5] = 1
    b3_in[6] = 1
    b3_in_d[7] = 1
    
    Q1 = np.einsum('i, j-> ji', chi(_omega, _omega_j[0], _Gamma),b1_in) + np.einsum('i,j -> ji', np.conj(chi(-_omega, _omega_j[0], _Gamma)),b1_in_d)
    Q2 = np.einsum('i, j-> ji', chi(_omega, _omega_j[1], _Gamma),b2_in) + np.einsum('i,j -> ji', np.conj(chi(-_omega, _omega_j[1], _Gamma)),b2_in_d)
    Q3 = np.einsum('i, j-> ji', chi(_omega, _omega_j[2], _Gamma),b3_in) + np.einsum('i,j -> ji', np.conj(chi(-_omega, _omega_j[2], _Gamma)),b3_in_d)
    return [Q1, Q2, Q3]

### normalization factor
def M(_omega, _omega_j, _detuning, _phi, _Gamma, _kappa, _g):
    """Calculates the normalization factor
    
    Parameters
    ----------
    _omega : 1D numpy array
        The frequency range of which chi shall be calculated
    _omega_j : numpy array of length 3
        mechanical frequencies
    _detuning : float
        Detuning
    _phi : np.array
        [0,0,pi/2]
    _Gamma : float
        Damping (either $\Gamma$ or $\kappa$)
    _kappa : float
        linewidth of cavity
    _g : np.array
        Couplings (g_x, g_y, g_z, g_xy, g_yz, g_zx)
        
    Returns
    -------
     : 2D np.array 
         M(omega, mode)
    """
    M1 = 1+ _g[0]**2 *mu(_omega, _omega_j, _Gamma)[0]*eta(_omega, _detuning, _phi[0], _kappa)
    M2 = 1+ _g[1]**2 *mu(_omega, _omega_j, _Gamma)[1]*eta(_omega, _detuning, _phi[1], _kappa)
    M3 = 1+ _g[2]**2 *mu(_omega, _omega_j, _Gamma)[2]*eta(_omega, _detuning, _phi[2], _kappa)
    return [M1, M2, M3]

### displacement operator
def q_1D(_omega, _omega_j, _detuning, _g, _Gamma, _kappa, _phi):
    """Calculates the operator q_j $\propto$(b_j+b_j^\dagger) without taking into account the 3D contributions
    
    Parameters
    ----------
    _omega : 1D numpy array
        The frequency range of which chi shall be calculated
    _omega_j : numpy array of length 3
        mechanical frequencies
    _detuning : float
        Detuning
    _g : np.array
        Couplings (g_x, g_y, g_z, g_xy, g_yz, g_zx)
    _Gamma : float
        Damping (either $\Gamma$ or $\kappa$)
    _kappa : float
        linewidth of cavity
    _phi : np.array
        [0,0,pi/2]
        
    Returns
    -------
     : 2D np.array 
         q(omega, mode) 1D
    """
    _M = M(_omega, _omega_j, _detuning, _phi, _Gamma, _kappa, _g)
    _Q_mech = Q_mech(_omega, _omega_j, _Gamma)
    _mu = mu(_omega, _omega_j, _Gamma)
    
    q1 = np.sqrt(_Gamma)*np.einsum('i,ji -> ji',1/_M[0], _Q_mech[0]) + 1j*np.sqrt(_kappa)*_g[0]*np.einsum('i, i, ji->ji',1/_M[0], _mu[0], Q_opt(_omega, _detuning, _kappa, _phi[0]))
    q2 = np.sqrt(_Gamma)*np.einsum('i,ji -> ji',1/_M[1], _Q_mech[1]) + 1j*np.sqrt(_kappa)*_g[1]*np.einsum('i, i, ji->ji',1/_M[1], _mu[1], Q_opt(_omega, _detuning, _kappa, _phi[1]))
    q3 = np.sqrt(_Gamma)*np.einsum('i,ji -> ji',1/_M[2], _Q_mech[2]) + 1j*np.sqrt(_kappa)*_g[2]*np.einsum('i, i, ji->ji',1/_M[2], _mu[2], Q_opt(_omega, _detuning, _kappa, _phi[2]))

    return [q1, q2, q3]

def q_3D(_omega, _omega_j, _detuning, _g, _Gamma, _kappa, _phi):
    """Calculates the operator q_j $\propto$(b_j+b_j^\dagger) with taking into account the 3D contributions
    
    Parameters
    ----------
    _omega : 1D numpy array
        The frequency range of which chi shall be calculated
    _omega_j : numpy array of length 3
        mechanical frequencies
    _detuning : float
        Detuning
    _g : np.array
        Couplings (g_x, g_y, g_z, g_xy, g_yz, g_zx)
    _Gamma : float
        Damping (either $\Gamma$ or $\kappa$)
    _kappa : float
        linewidth of cavity
    _phi : np.array
        [0,0,pi/2]
        
    Returns
    -------
     : 2D np.array 
         q(omega, mode) 3D
    """

    _M = M(_omega, _omega_j, _detuning, _phi, _Gamma, _kappa, _g)
    _mu = mu(_omega, _omega_j, _Gamma)
    q = q_1D(_omega, _omega_j, _detuning, _g, _Gamma, _kappa, _phi)
    
    # coupling terms
    GXY = 1j*eta(_omega, _detuning, 0, _kappa)*_g[0]*_g[1] + _g[3]
    GYX = GXY
    GYZ = -1j*eta(_omega, _detuning, np.pi/2, _kappa)*_g[1]*_g[2] + _g[4]
    GZY = 1j*eta(_omega, _detuning, np.pi/2, _kappa)*_g[2]*_g[1] + _g[4]
    GXZ = -1j*eta(_omega, _detuning, np.pi/2, _kappa)*_g[0]*_g[2] + _g[5]
    GZX = 1j*eta(_omega, _detuning, np.pi/2, _kappa)*_g[2]*_g[0] + _g[5]
    
    # prefactors of 3D contributions
    RXY = 1j * _mu[0] * GXY / _M[0]
    RYX = 1j * _mu[1] * GYX / _M[1]
    RXZ = 1j * _mu[0] * GXZ / _M[0]
    RZX = 1j * _mu[2] * GZX / _M[2]
    RYZ = 1j * _mu[1] * GYZ / _M[1]
    RZY = 1j * _mu[2] * GZY / _M[2]
    
    # 3D modes
    Denum = RXY*(RYX+RYZ*RZX) + RXZ*(RYX*RZY+RZX) + RYZ*RZY -1
    
    q1 = (q[0]*(RYZ*RZY-1) - q[1]*(RXY+RXZ*RZY) - q[2]*(RXY*RYZ+RXZ)) / Denum
    q2 = (-q[0]*(RYX+RYZ*RZX) + q[1]*(RXZ*RZX-1) - q[2]*(+RXZ*RYX+RYZ)) / Denum
    q3 = (-q[0]*(RYX*RZY+RZX) - q[1]*(RXY*RZX+RZY) + q[2]*(RXY*RYX-1)) / Denum
    
    return [q1, q2, q3]

    
### helpers
def expectation_value(_operator, _n, _pair):
    """Calculates the expectation value of an operator by analyzing the noises
    
    Parameters
    ----------
    _operator : np.array
        operator as function of omega (containing all directions)
    _n : float
        Expectation value of the respective noise
    _pair : integer
        select pair (0=photon, 1=x, 2=y, 3=z)
    
    Returns
    -------
     : np.array
        <operator>(omega)
    """
    return (_n+1)* np.abs(_operator[2*_pair])**2 + _n * np.abs(_operator[2*_pair+1])**2

def spectrum(_operator, _n_opt, _n_mech):
    """Calculates the PSD
    
    Parameters
    ----------
    _operator : np.array
        operator as function of omega (containing all directions)
    _n_opt : float
        optical photon number (n_opt=0)
    _n_mech : np.array
        phonon numbers (n_x, n_y, n_z)
    
    Returns
    -------
     : np.array
        <operator>_total(omega) (sum over all modes)
    """
    s_a = expectation_value(_operator, _n_opt, 0)
    s_b1 = expectation_value(_operator, _n_mech[0], 1)
    s_b2 = expectation_value(_operator, _n_mech[1], 2)
    s_b3 = expectation_value(_operator, _n_mech[2], 3)
    return s_a + s_b1 + s_b2 + s_b3


def spectrum_output(omega, _i, param, ThreeD):
    """Calculates the PSD for a given omega regime and set of parameters
    
    Parameters
    ----------
    omega: np.array
        Frequency range in which the spectrum is to be computed
    _i : integer
        selection of operator (0=photon, 1=x, 2=y, 3=z)
    param : class param
        set of parameters
    ThreedD : boolean
        Consider 3D contribution (True) or not (False)
        
    Returns
    -------
     : np.array
         PSD(omega)
    """
    
    # define phases and n_opt
    _phi = np.array([0, 0, np.pi/2])
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
    """Calculates area under curve by using the trapezoidal rule

    Parameters
    ----------
    _S : np.array
        spectrum
    _Delta : float
        spacing of omega
    
    Returns
    -------
     : float
        Area under the spectrum
    """

    summe=0
    for i in range(len(_S)-1):
        Tem = 0.5*(_S[i]+_S[i+1])
        summe = summe + Tem
    return summe*_Delta / (2*np.pi) #/2pi for normalization
   
    
def n_from_area(_S_plus, _S_minus, _Delta_omega, _N = 0, _name = '', printing = True):
    """Calculates phonon number from area and compares it to the one from the formula
    
    Parameters
    ----------
    _S_plus : np.array
        Spectrum for positive omega
    _S_minus : np.array
        Spectrum for negative omega
    _Delta_omega : float
        Spacing of omega
    _N : float
        Phonon number from formula
    _name : str
        Name of respective operator (x, y or z)
    printing : boolean
        Print the result (True), default is True
    
    Returns
    -------
     : list
        Phonon numbers (N_plus, N_minus, N_total)
    """
    N_X_plus = area(_S_plus, _Delta_omega)
    N_X_minus = area(_S_minus, _Delta_omega)
    N_X = (N_X_plus + N_X_minus -1)/2
    
    if printing == True:
        print()
        print(_name, 'photon number from area, difference')
        print('+:', round(N_X_plus, 2), '(', round((N_X_plus-_N)/(_N+1)*100, 2), '% )')
        print('-:', round(N_X_minus, 2), '(', round((N_X_minus-_N)/(_N)*100, 2), '% )')
        print('summed:', round(N_X, 2), '(', round((N_X-_N)/(_N)*100, 2), '% )')
    
    return [N_X_plus, N_X_minus, N_X]


def photon_number(_n_j, _Gamma_opt, _Gamma, printing = True):
    """Calculates the phonon number from the formula
    
    Parameters
    ----------
    _n_j : np.array
        phonon numbers at room temperature
    _Gamma_opt : np.array
        optical damping rate (x,y,z)
    _Gamma : float
        mechanical damping rate
    printing : boolean
        Print the result (True), default is True
    
    Returns
    -------
     : np.array
        Phonon numbers (N_plus, N_minus, N_total)
    """
    N = 2*_n_j * _Gamma / (abs(_Gamma_opt) + 2*_Gamma)
    
    if printing == True:
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
    """This class contains all relevant parameters
    """
    n_opt = 0 #: Photon number at room temperature
    T = 300 # room temperature [K]
    RHO = 2198 #: sphere density [kg/m^3]
    EPSR = 2.1 #: relative permittivity [F m^-1]
    lambda_tw = 1064e-9 #: wavelength of tweezer [m]   
    waist = 41.1e-6 #:  waist radius [m]
    WX = 0.67e-6 #: focus of tweezer in x-direction [m]
    WY = 0.77e-6 #: focus of tweezer in y-direction [m] 
    XL = 1.07e-2 #: cavity length [m]
    DelFSR = 14.0e9 #:Free spectral range [Hz], not used if couplings are given
    Y0 = 0 #: y_0, equilibrium position in x-direction
    Z0 = 0 #: z_0, equilibrium position in x-direction
    R0= 0.5*143e-9 #: sphere radius [m]
    Finesse = 73e3 #: Finesse
    Press=1e-6 #: air pressure [mbar]
    Pin1= 0.4 #: input power tweezer beam [W]
    detuning=-300e3 #: detuning of trap beam (omega_cav - omega_tw) [2pi kHz]
    theta0= 0.25 #: angle between tweezer polarization and cavity axis [pi]
    X0 = 0.23*1064e-9  #: equilibrium position in x [m]
    
    def print_param(self):
        """Prints all the parameters in a nice fashion"""
        
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
        print('ZR [m]', self.ZR)
        print('Gamma: ', self.Gamma)
        print('omega_mech/2pi [kHz]: ', self.omega_mech/(2*np.pi)*1e-3)
        print('photons in cavity: ', self.n_photon)
        print('GX, GY, GZ [2pi Hz]', self.g[0]/2/np.pi, self.g[1]/2/np.pi, self.g[2]/2/np.pi)
        print('GXY, GYZ, GZX [2pi Hz]', self.g[3]/2/np.pi, self.g[4]/2/np.pi, self.g[5]/2/np.pi)
        print('***************')
        
    def prepare_calc(self):
        """Calculates all the theoretical relevant parameters if only the experimental ones are given
        
        Warning
        ------
        Detuning has to be given in 2pi Hz
        """
        
        # mass
        self.XM=self.RHO*4.*np.pi/3.*self.R0**3
        
        # polarisibility
        self.Polaris=4.*np.pi*Epsi0*(self.EPSR-1)/(self.EPSR+2)*self.R0**3 # from Tania
        
        # tweezer properties
        WK= 2*np.pi / self.lambda_tw #=2*pi/lambda=k
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
        self.ZR = ZR
        
        # linewiddth
        coeff=WK*self.Polaris/Epsi0/OMOPT**2/np.pi
        kappnano=4*coeff**2*self.DelFSR*np.cos(WK*self.X0)*np.cos(WK*self.X0)
        self.kappa=kappnano+KAPPin
         
        # damping rate
        GAMMAM=1600*self.Press/np.pi
        GAMMAM=GAMMAM/500/self.RHO/self.R0
        self.Gamma=GAMMAM/2 # our Gamma is Gammam/2
        
        # mechanical frequencies
        Det2pi=self.detuning*2*np.pi
        kapp2=0.5*self.kappa
        Wkx0=WK*self.X0
        OmX=self.Polaris*self.epsTW**2/self.XM/self.WX**2
        OmY=self.Polaris*self.epsTW**2/self.XM/self.WY**2
        OmZ=0.5*self.Polaris*self.epsTW**2/self.XM/ZR**2
        
        # theta[rad]
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
     
        
    def opt_damp_rate(self, printing = False):
        """Calculates the optical damping rate
        
        Parameters
        ----------
        printing : boolean
            Result is printed (True) or not, default True
            
        Returns
        -------
         : np.array
            Optical damping rate for (x,y and z)
        
        Warning
        -------
        Detuning has to be given in 2pi Hz
        """
        Det2pi = self.detuning * 2*np.pi
        g_1D = np.array([self.g[0], self.g[1], self.g[2]])
        Gamma_opt = -g_1D**2 * self.kappa * (1/(self.kappa**2/4 + (Det2pi+self.omega_mech)**2) - 1/(self.kappa**2/4 + (Det2pi-self.omega_mech)**2) )  
    
        if printing == True:
            print()
            print('optical cooling rates')
            print('$\Gamma_{opt, x}$:', round(Gamma_opt[0]/1e5, 7), '1e5')
            print('$\Gamma_{opt, y}$:', round(Gamma_opt[1]/1e5, 7), '1e5')
            print('$\Gamma_{opt, z}$:', round(Gamma_opt[2]/1e5, 7), '1e5')
            
        return Gamma_opt

def loop_progress(L_inner, L_outer, inner, outer, start_time):
    """Print nice progress control in terminal
    
    Parameters
    ----------
    L_inner : integer
        length of inner loop
    L_outer : integer
        length of outer loop
    inner : integer
        current value of loop parameter of inner loop
    outer : integer
        current value of loop parameter of outer loop
    start_time : float
        time when loops where started
    """
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
    print('\n completed: {0:.3f}%, running: {1:6}s, remaining: {2:6}s \r'.format(progress*100, current_time_form, remaining_time_form))
    
    
    