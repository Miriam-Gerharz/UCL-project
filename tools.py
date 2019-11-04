import numpy as np


### operators (only 1D by now)
def a_in(_omega):
    array = np.array([[0 for j in range(len(_omega))] for i in range(8)])    
    array[0] = 1 * np.ones(len(_omega))
    return array

def a_in_d(_omega):
    array = np.array([[0 for j in range(len(_omega))] for i in range(8)])    
    array[1] = 1 * np.ones(len(_omega))
    return array

def b1_in(_omega):
    array = np.array([[0 for j in range(len(_omega))] for i in range(8)])    
    array[2] = 1 * np.ones(len(_omega))
    return array

def b1_in_d(_omega):
    array = np.array([[0 for j in range(len(_omega))] for i in range(8)])    
    array[3] = 1 * np.ones(len(_omega))
    return array

def b2_in(_omega):
    array = np.array([[0 for j in range(len(_omega))] for i in range(8)])    
    array[4] = 1 * np.ones(len(_omega))
    return array

def b2_in_d(_omega):
    array = np.array([[0 for j in range(len(_omega))] for i in range(8)])    
    array[5] = 1 * np.ones(len(_omega))
    return array

def b3_in(_omega):
    array = np.array([[0 for j in range(len(_omega))] for i in range(8)])    
    array[6] = 1 * np.ones(len(_omega))
    return array

def b3_in_d(_omega):
    array = np.array([[0 for j in range(len(_omega))] for i in range(8)])    
    array[7] = 1 * np.ones(len(_omega))
    return array

### SUSCEPTIBILITIES ###
# chi
def chi(_omega, _omega_j, _Gamma):
    return 1/(-1j*(_omega-_omega_j) + _Gamma/2.0)

# optical
def eta(_omega, detuning, _phi, _kappa):
    return np.exp(-1j*_phi)*chi(_omega, -detuning, _kappa) - np.exp(1j*_phi)*np.conj(chi(-_omega, -detuning, _kappa))
                  
# mechanical 
def mu(_omega, _omega_j, _Gamma):
    mu1 = chi(_omega, _omega_j[0], _Gamma[0]) - np.conj(chi(-_omega, _omega_j[0], _Gamma[0]))
    mu2 = chi(_omega, _omega_j[1], _Gamma[1]) - np.conj(chi(-_omega, _omega_j[1], _Gamma[1]))
    mu3 = chi(_omega, _omega_j[2], _Gamma[2]) - np.conj(chi(-_omega, _omega_j[2], _Gamma[2]))
    return [mu1, mu2, mu3]

### NOISES ###
# optical
def Q_opt(_omega, _detuning, _kappa, _phi):
#_omega: array
#        freq range in which spectrum is to be computed
#_phi: float
#        phase
# 
    return np.exp(-1j*_phi)*np.einsum('i, ji', chi(_omega, -_detuning, _kappa), a_in(_omega)) + np.exp(1j*_phi)*np.einsum('i, ji', np.conj(chi(-_omega, -_detuning, _kappa)),a_in_d(_omega))

# mechanical
def Q_mech(_omega, _omega_j, _Gamma, _phi):
    Q1 = np.einsum('i, ji-> ji', chi(_omega, _omega_j[0], _Gamma[0]),b1_in(_omega)) + np.einsum('i,ji -> ji',np.conj(chi(-_omega, _omega_j[0], _Gamma[0])),b1_in_d(_omega))
    Q2 = np.einsum('i, ji-> ji', chi(_omega, _omega_j[1], _Gamma[1]),b2_in(_omega)) + np.einsum('i,ji -> ji',np.conj(chi(-_omega, _omega_j[1], _Gamma[1])),b2_in_d(_omega))
    Q3 = np.einsum('i, ji-> ji', chi(_omega, _omega_j[2], _Gamma[2]),b3_in(_omega)) + np.einsum('i,ji -> ji',np.conj(chi(-_omega, _omega_j[2], _Gamma[2])),b3_in_d(_omega))
    return [Q1, Q2, Q3]

### normalization factor
def M(_omega, _omega_j, _detuning, _phi, _Gamma, _kappa, _g):
    M1 = 1+ _g[0]**2 *mu(_omega, _omega_j, _Gamma)[0]*eta(_omega, _detuning, _phi[0], _kappa)
    M2 = 1+ _g[1]**2 *mu(_omega, _omega_j, _Gamma)[1]*eta(_omega, _detuning, _phi[1], _kappa)
    M3 = 1+ _g[2]**2 *mu(_omega, _omega_j, _Gamma)[2]*eta(_omega, _detuning, _phi[2], _kappa)
    return [M1, M2, M3]

### displacement operator
def q(_omega, _omega_j, _detuning, _g, _Gamma, _kappa, _phi):
    _M = M(_omega, _omega_j, _detuning, _phi, _Gamma, _kappa, _g)
    _Q_mech = Q_mech(_omega, _omega_j, _Gamma, _phi)
    _mu = mu(_omega, _omega_j, _Gamma)
    
    q1 = np.sqrt(_Gamma[0])*np.einsum('i,ji -> ji',1/_M[0], _Q_mech[0]) + 1j*np.sqrt(_kappa)*_g[0]*np.einsum('i, i, j->ji',1/_M[0], _mu[0], Q_opt(_omega, _detuning, _kappa, _phi[0]))
    q2 = np.sqrt(_Gamma[1])*np.einsum('i,ji -> ji',1/_M[1], _Q_mech[1]) + 1j*np.sqrt(_kappa)*_g[1]*np.einsum('i, i, j->ji',1/_M[1], _mu[1], Q_opt(_omega, _detuning, _kappa, _phi[1]))
    q3 = np.sqrt(_Gamma[2])*np.einsum('i,ji -> ji',1/_M[2], _Q_mech[2]) + 1j*np.sqrt(_kappa)*_g[2]*np.einsum('i, i, j->ji',1/_M[2], _mu[2], Q_opt(_omega, _detuning, _kappa, _phi[2]))

    
    return [q1, q2, q3]
    
### helper
def expectation_value(_operator, _n, _pair):
    return _n* np.abs(_operator[_pair])**2 + (_n+1)* np.abs(_operator[_pair+1])**2

def spectrum_2(_operator, _n_opt, _n_mech):
    s_a = expectation_value(_operator, _n_opt, 0)
    s_b1 = expectation_value(_operator, _n_mech[0], 1)
    s_b2 = expectation_value(_operator, _n_mech[1], 2)
    s_b3 = expectation_value(_operator, _n_mech[2], 3)
    return s_a + s_b1 + s_b2 + s_b3


def spectrum_output(_i, param):
    # define phases
    _phi = np.array([0,0,np.pi/2])
    
    # unpack parameters
    omega = param[0]
    omega_j = param[1]
    detuning = param[2]
    g = param[3]
    Gamma = param[4]
    kappa = param[5]
    n_opt = param[6]
    n_mech = param[7]
    
    #print(param[4])
    
    
    # calculate q operator
    operator = q(omega, omega_j, detuning, g, Gamma, kappa, _phi)[_i]
    
    _spectrum = spectrum_2(operator, n_opt, n_mech)
    return g[_i]**2 * np.abs(eta(omega, detuning, _phi[_i], kappa))**2 * _spectrum


