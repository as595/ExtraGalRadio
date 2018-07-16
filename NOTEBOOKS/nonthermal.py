import numpy as np

const_e = 1.6e-19	# Joules [SI units]
m_e_kg = 9.1e-31	# kg 	 [SI units]
eV2erg = 1.602e-12  # erg    [CGS units]

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

def n_E_erg(E):

	if (type(E) is np.float64):
		if (E>6.4e-4):
			n_E = (E/6.4e-4)**(-2.8)
		if ((E>8.0e-5) and (E<6.4e-4)):
			n_E = (E/6.4e-4)**(-1.8)
		if (E<8.0e-5):
			n_E = (8.0e-5/6.4e-4)**(-1.8)
	else:
		n_E = []
		for val in E:
			if (val>6.4e-4):
				n_E.append((val/6.4e-4)**(-2.8))
			elif ((val>8.0e-5) and (val<400.)):
				n_E.append((val/6.4e-4)**(-1.8))
			elif (val<8.0e-5):
				n_E.append((8.0e-5/6.4e-4)**(-1.8))

		n_E = np.array(n_E)

	
	return n_E

# ------------------------------------------------------------------------------

def dn_E_E2_erg(E):

	if (type(E) is np.float64):
		if (E>6.4e-4):
			n_E = (-4.8/(6.4e-4)**3)*(E/6.4e-4)**(-5.8)
		if ((E>8.0e-5) and (E<6.4e-4)):
			n_E = (-3.8/(6.4e-4)**3)*(E/6.4e-4)**(-4.8)
		if (E<8.0e-5):
			n_E = -2.*(8.0e-5/6.4e-4)**(-1.8)*(E**-3)
	else:
		n_E = []
		for val in E:
			if (val>6.4e-4):
				n_E.append((-4.8/(6.4e-4)**3)*(val/6.4e-4)**(-5.8))
			elif ((val>8.0e-5) and (val<6.4e-4)):
				n_E.append((-3.8/(6.4e-4)**3)*(val/6.4e-4)**(-4.8))
			elif (val<8.0e-5):
				n_E.append(-2.*(8.0e-5/6.4e-4)**(-1.8)*(E**-3))

		n_E = np.array(n_E)

	
	return n_E

# ------------------------------------------------------------------------------

def n_E_AGN_erg(E):
    
    E0 = 100e6*eV2erg
    #E0 = 10e6*eV2erg
    if (type(E) is np.float64 or float):
        if (E>=E0):
            n_E = (E/E0)**(-2.5)
        if (E<E0):
            n_E = 0.
    else:
        n_E = []
        for val in E:
            if (val>=E0):
                n_E.append((val/E0)**(-2.5))
            elif (val<E0):
                n_E.append(0.0)

    n_E = np.array(n_E)


    return n_E

# ------------------------------------------------------------------------------

def dn_E_E2_AGN_erg(E):
    
    E0 = 100e6*eV2erg
    #E0 = 10e6*eV2erg
    if (type(E) is np.float64 or float):
        if (E>=E0):
            n_E = (-4.5/(E0)**3)*(E/E0)**(-5.5)
        if (E<E0):
            n_E = 0.
    else:
        n_E = []
        for val in E:
            if (val>=E0):
                n_E.append((-4.5/(E0)**3)*(val/E0)**(-5.5))
            elif (val<E0):
                n_E.append(0.0)

    n_E = np.array(n_E)

    return n_E

# ------------------------------------------------------------------------------

def L_abs_ML(nu,H_T,A_m2):

	"""
	Calculate SSA luminosity following Longair Eq 18.61 and P&B96 Eq. 9 for conversion
	to luminosity.

	Variables are in SI units.

	nu in Hertz.
	"""

	nu_g = const_e*H_T/(2*np.pi*m_e_kg)
	
	L = (8.*np.pi/3.)*m_e_kg*A_m2*nu**(2.5)*nu_g**(-0.5)

	return L


