# ------------------------------------------------------------------------------

def tau_ff(nu,ne,Te,dl):

	"""
	Cane 1979: Equation in Sec 7(c)
	"""

	EM = ne**2*dl
	tau = 1.64e5*Te**(-1.35)*(nu/1e6)**(-2.1)*EM

	return tau

# ------------------------------------------------------------------------------

def kappa_ff(nu,ne,Te):

	"""
	Cane 1979: Equation in Sec 7(c)
	"""

	kappa = 1.64e5*Te**(-1.35)*(nu/1e6)**(-2.1)*ne**2

	return kappa
