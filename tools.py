import numpy as np


'''
comments:
index i: function of omega
index j: different operators

'''


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
    #_phi = np.array([0,0,0])
    
    # unpack parameters
    omega_j = param[0]
    detuning = param[1]
    g = param[2]
    Gamma = param[3]
    kappa = param[4]
    n_opt = param[5]
    n_mech = param[6]
    
    # calculate q operator
    if ThreeD == False:
        operator = q_1D(omega, omega_j, detuning, g, Gamma, kappa, _phi)[_i]
    if ThreeD == True:
        operator = q_3D(omega, omega_j, detuning, g, Gamma, kappa, _phi)[_i]
    # calculate spectrum
    _spectrum = spectrum(operator, n_opt, n_mech)
    
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
   
    

def n_from_area(_S_plus, _S_minus, _Delta, _N, _name):
    N_X_plus = area(_S_plus, _Delta)
    N_X_minus = area(_S_minus, _Delta)
    N_X = (N_X_plus + N_X_minus -1)/2
    
    print()
    print(_name, 'photon number from area, difference')
    print('+:', round(N_X_plus, 2), '(', round((N_X_plus-_N)/(_N+1)*100, 2), '% )')
    print('-:', round(N_X_minus, 2), '(', round((N_X_minus-_N)/(_N)*100, 2), '% )')
    print('summed:', round(N_X, 2), '(', round((N_X-_N)/(_N)*100, 2), '% )')
    #print('area -', area(SXX_minus, omega))
    return [N_X_plus, N_X_minus, N_X]

def opt_damp_rate(_kappa, _detuning, _g, _omega_j):
    Det2pi = _detuning * 2*np.pi
    g_1D = np.array([_g[0], _g[1], _g[2]])
    Gamma_opt = -g_1D**2 * _kappa * (1/(_kappa**2/4 + (Det2pi+_omega_j)**2) - 1/(_kappa**2/4 + (Det2pi-_omega_j)**2) )  

    print()
    print('optical cooling rates')
    print('$\Gamma_{opt, x}$:', round(Gamma_opt[0]/1e5, 3), '1e5')
    print('$\Gamma_{opt, y}$:', round(Gamma_opt[1]/1e5, 3), '1e5')
    print('$\Gamma_{opt, z}$:', round(Gamma_opt[2]/1e5, 3), '1e5')
    
#    print('mechanical damping')
#    print('X', Gamma[0])
    
    return Gamma_opt


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
