import numpy as np
import scipy.special as ss
from scipy.integrate import quad

from integrals import *
from nonthermal import *

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

const_c1 = 6.27e18		# page 85-86
const_c2 = 2.37e-3		# page 86 [for electrons & positrons]
const_c3 = 1.87e-23		# page 90
const_c4 = 4.20e7		# page 93
const_c5 = 9.68e-24     # page 232  gamma dependent
const_c6 = 8.10e-41     # page 232  gamma dependent

sin_th_av = 0.785
const_c_cgs = 2.99e10	# cm/s

# ------------------------------------------------------------------------------

def funcF(x):

    # Modified Bessel function second kind: K_{5/3}
    K = lambda z : ss.kv(1.67, z)

    if (type(x) is np.float64):
        integral = logquad(K, x, 1e3)[0]
        F = x*integral
    else:
        F=[]
        for xmin in x:
            integral = logquad(K, xmin, 1e3)[0]
            F.append(xmin*integral)

    F = np.array(F)

    return F

# ------------------------------------------------------------------------------

def E_integrand(E_erg,nu,H_G):


	#E_J  = E*1e6*1.6e-19		# MeV ---> J
	#E_erg = E_J*1e7 			# J ---> erg
	
	x = nu / (const_c1*H_G*sin_th_av*E_erg**2)	# CGS units

	#integrand = n_E(E)
	integrand = n_E_erg(E_erg)
	integrand*= funcF(x)
	
	return integrand

# ------------------------------------------------------------------------------

def E_integrand2(E_erg,nu,H_G):


	#E_J = E*1e6*1.6e-19		#MeV ---> J
	#E_erg = E_J*1e7 			# J ---> erg
	
	x = nu / (const_c1*H_G*sin_th_av*E_erg**2)	# CGS units

	integrand = -1.*dn_E_E2_erg(E_erg)
	integrand*= funcF(x)
	integrand*= E_erg**2
	
	return integrand

# ------------------------------------------------------------------------------

def E_integrand_AGN(E_erg,nu,H_G):
    
    
    #E_J  = E*1e6*1.6e-19        # MeV ---> J
    #E_erg = E_J*1e7             # J ---> erg
    
    x = nu / (const_c1*H_G*sin_th_av*E_erg**2)    # CGS units
    
    #integrand = n_E(E)
    integrand = n_E_AGN_erg(E_erg)
    integrand*= funcF(x)
    
    return integrand

# ------------------------------------------------------------------------------

def E_integrand2_AGN(E_erg,nu,H_G):
    
    
    #E_J = E*1e6*1.6e-19        #MeV ---> J
    #E_erg = E_J*1e7             # J ---> erg
    
    x = nu / (const_c1*H_G*sin_th_av*E_erg**2)    # CGS units
    
    integrand = -1.*dn_E_E2_AGN_erg(E_erg)
    integrand*= funcF(x)
    integrand*= E_erg**2
    
    return integrand

# ------------------------------------------------------------------------------

def calc_em_synch(nu,H_G,N0):

	em_nu = []
	for freq in nu:

	    E_int = logquad(E_integrand,1e-9,1e14,args=(freq,H_G))[0]
	    em_nu.append(0.5*const_c3*H_G*sin_th_av*N0*E_int)

	em_nu = np.array(em_nu)

	return em_nu

# ------------------------------------------------------------------------------

def calc_kappa_synch(nu,H_G,N0):

	kappa_nu = []
	for freq in nu:

		E_int2 = logquad(E_integrand2,1e-9,1e14,args=(freq,H_G))[0]
		kappa_nu.append(0.5*const_c_cgs**2*const_c3*freq**(-2)*H_G*sin_th_av*N0*E_int2)
		
	kappa_nu = np.array(kappa_nu)

	return kappa_nu

# ------------------------------------------------------------------------------

def calc_em_synch_AGN(nu,H_G,N0):
    
    """
    For AGN, the cosmic ray electron energy sp[ectrum starts
    at 100 MeV = 1.6e-4 erg
    """
    
    em_nu = []
    for freq in nu:
        
        E_int = logquad(E_integrand_AGN,1e-19,1e14,args=(freq,H_G))[0]
        em_nu.append(const_c3*H_G*sin_th_av*N0*E_int)  # removed a pre-factor of 0.5 here
    
    em_nu = np.array(em_nu)

    return em_nu

# ------------------------------------------------------------------------------

def calc_kappa_synch_AGN(nu,H_G,N0):
    
    """
    For AGN, the cosmic ray electron energy sp[ectrum starts
    at 100 MeV = 1.6e-4 erg
    """
    
    kappa_nu = []
    for freq in nu:
        
        E_int2 = logquad(E_integrand2_AGN,1e-19,1e14,args=(freq,H_G))[0]
        kappa_nu.append(0.5*const_c_cgs**2*const_c3*freq**(-2)*H_G*sin_th_av*N0*E_int2)
    
    kappa_nu = np.array(kappa_nu)

    return kappa_nu


# ------------------------------------------------------------------------------

def eq_3pt50(nu,H_G,N0):

    gamma = 2.5
    em = const_c5*N0*(H_G*sin_th_av)**(0.5*(gamma+1))*(0.5*nu/const_c1)**(0.5*(1-gamma))

    return em

# ------------------------------------------------------------------------------

