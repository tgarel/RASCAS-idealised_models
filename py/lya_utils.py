import numpy as np

# physical constants, all in cgs.
kb_cgs    = 1.3806200e-16  # [erg/K]
h_cgs     = 6.6260755e-27  # [erg s]
mp_cgs    = 1.6726231e-24  # [g]
mh_cgs    = 1.6749286e-24  # [g]
clight    = 2.99792458e10  # [cm/s]
lambda0   = 1215.6701e-8   # [cm]
nu0       = clight/lambda0 # [Hz]
A21       = 6.265e8        # [s^-1]
f12       = 0.416          # [-]
me_cgs    = 9.1093897e-28  # [g]
e_cgs     = 4.8032068e-10  # [esu] election charge
# conversion factors 
pc2cm     = 3.08567758e18 

# ----------------------------------------------
# all usefull functions below assume cgs units.
# ----------------------------------------------
def sigma_0(temperature):
    # Eq. 15 from Dijkstra 2014. 
    dnud = delta_nu_d(temperature)
    s    = f12 * np.sqrt(np.pi) * e_cgs**2 / dnud / me_cgs / clight
    return s

def thermal_velocity(temperature):
    v = np.sqrt( 2.*kb_cgs*temperature/mp_cgs ) # cm/s  
    return v

def delta_nu_d(temperature):
    vth = thermal_velocity(temperature)
    dnud = vth / lambda0
    return dnud

def voigt_param(temperature): 
    a = A21 / 4. / np.pi / delta_nu_d(temperature)
    return a

def x_from_nu(nu,temperature):
    x = (nu - nu0) / delta_nu_d(temperature)
    return x

def x_from_v(v,temperature):
    x = - v / thermal_velocity(temperature) # Verhamme, Schaerer, Maselli 2006, Eq. 4.
    return x

def v_from_x(x,temperature): 
    v = - x * thermal_velocity(temperature) 
    return v 

def v_from_nu(nu):
    v = lambda0 * (nu0 - nu)
    return v

def nu_from_x(x,temperature):
    nu = x * delta_nu_d(temperature) + nu0
    return nu

def nu_from_v(v):
    nu = nu0 * (1. - v / clight)
    return nu




 
