# Calculate gamma-ray interaction length for photon-photon pair production on the cosmic microwave background.
# [180607 - Justin Bray] Original code
# [180608 - Anna Scaife] Updated to include scipy integration [not always successfully]
# [180618 - Justin Bray] Implemented efficient scipy integration; pruned dead code; changed structure of functions underneath get_intlen to be a bit more general.
# [180618 - Justin Bray] Implemented muon-muon pair production.



import pylab as pl
import numpy as np
from scipy.integrate import quad

# ------------------------------------------------------------------------------

# Physical constants.
m_e = 9.11e-31 # kg; electron mass
m_mu = 1.88e-28 # kg; muon mass
q_e = 1.67e-19 # C; electron charge
c = 3e8 # m/s; speed of light
sigma_t = 6.65e-29 # m^2; Thomson cross-section
T_CMB = 2.725 # K; temperature of cosmic microwave background
h = 6.626e-34 # Js; Planck constant
k = 1.38e-23 # J/K; Boltzmann constant
Mpc_in_m = 3.086e22 # m; 1 Mpc

# ------------------------------------------------------------------------------

def logquad(func, a, b, **kwargs):

  """
  Perform scipy.integrate.quad more efficiently on a function that's better-behaved in log space.
  """

  # Perform integral $\int y dx$ by transforming variables x->log(x), y->log(y).

  return quad(lambda x,*args: np.exp(x + np.log(func(np.exp(x),*args))), np.log(a), np.log(b), **kwargs)

# ------------------------------------------------------------------------------

def get_sigma_pp(s, m_l):

  """
  Calculate p-p pair production cross-section, per berestetskii1982,
  as cited in eqns 1&2 of
  https://www.astro.umd.edu/~miller/teaching/astr688m/lecture03.pdf
  """

  v = np.sqrt(1. - (m_l*c**2)**2 / (s/4.))
  sigma_pp = 3./16. * sigma_t * (m_e/m_l)**2 * (1. - v**2) * ( (3.-v**4)\
            * np.log( (1.+v)/(1.-v) ) - 2.*v*(2.-v**2)) # m^2; p-p cross-section
  # Note: function above has been modified to scale sigma_t by (m_e/m_l)**2.
  # Justification: Thomson cross-section has a 1/m_e^2 factor in it, so this
  #  scales correctly for other lepton masses.

  return sigma_pp

# ------------------------------------------------------------------------------

def s_integrand(s, m_l):

  """
  Define the final integrand, over s, from eqn 23 of protheroe1996.
  """

  integrand = s * get_sigma_pp(s, m_l) # J^2 m^2

  return integrand

# ------------------------------------------------------------------------------

def s_integral(E, eps, m_l):

  """
  Calculate the final integral, over s, from eqn 23 of protheroe1996.
  """

  # protheroe1996, eqn 26:
  smin = (2.*m_l*c**2)**2

  # J^2; maximum CoM energy^2; protheroe1996, eqn 28:
  smax = 4*E*eps

  # protheroe1996, eqn 23:
  integral = logquad(s_integrand, smin, smax, args=(m_l))[0] # J^4 m^2

  return integral

# ------------------------------------------------------------------------------

def eps_integrand(eps, E, get_Neps, m_l):

  """
  Define the integrand, over epsilon, from eqn 23 of protheroe1996.
  """

  # get the differential photon number density, protheroe1996, eqn 24:
  Neps = get_Neps(eps) # photons m^-3 J^-1

  integrand = Neps * eps**-2 * s_integral(E, eps, m_l) # J/m

  return integrand

# ------------------------------------------------------------------------------

def eps_integral(E, get_Neps, m_l):

  """
  Calculate the final integral, over epsilon, from eqn 23 of protheroe1996.
  """

  # protheroe1996, eqn 27:
  epsmin = (2*m_l*c**2)**2 / (4*E) # Joules
  epsmax = 1e10*epsmin # J; target photon energy above majority of effect
  #epsmax = np.inf

  # protheroe1996, eqn 23:
  integral = logquad(eps_integrand, epsmin, epsmax, args=(E, get_Neps, m_l))[0] # J^2/m

  return integral

# ------------------------------------------------------------------------------

def get_intlen_pp(E, get_Neps, m_l=m_e): # default lepton mass m_l is m_e

  """
  Find interaction length for pair production for photons of given energy,
  from eqn 23 of protheroe1996.
  """

  return 8*E**2 / eps_integral(E, get_Neps, m_l) # m

# ------------------------------------------------------------------------------

def get_intlen(E, get_Neps, use_pp=True, use_mu=False, use_4e=False):

  """
  Find interaction length for photons of given energy, using provided get_Neps
  function to determine photon background, and incorporating interaction types
  specified by flags (use_pp for pair production, use_muon for muon pair
  production, use_4e for double pair production).
  """

  if use_pp:
    intlen_pp = get_intlen_pp(E, get_Neps)
  else:
    intlen_pp = 1e99

  if use_mu:
    intlen_mu = get_intlen_pp(E, get_Neps, m_l=m_mu)
  else:
    intlen_mu = 1e99

  if use_4e:
    assert False, 'Double pair production not yet implemented.' # TODO: implement
  else:
    intlen_4e = 1e99

  return 1./(1./intlen_pp + 1./intlen_mu + 1./intlen_4e)

# ------------------------------------------------------------------------------

def get_Neps_CMB(eps):

  """
  Find differential photon number density of CMB at given photon energy.
  """

  nu = eps / h # Hz; photon frequency

  # Calculate spectral radiance of the CMB via Planck's Law.
  I = 2*h*nu**3 / c**2 * 1./(np.exp((h*nu)/(k*T_CMB)) - 1.) # W Hz^-1 m^-2 sr^-1
  # Calculate differential photon number density via eqn 24 of protheroe1996.
  Neps = (4*np.pi) / (h*c) * I / (h*nu) # photons m^-3 J^-1

  return Neps

# ------------------------------------------------------------------------------

def get_Neps_IR(eps):

  """
  Find differential photon number density from the optical and infrared
  background at a given photon energy.
  """

  return zeros(len(eps)) # TODO: proper optical/IR background

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

if __name__ == "__main__":

    # define array of energies:
    E = 10**np.linspace(13,22,20) # eV; photon energy
    E *= q_e # eV->J

    # calculate mean interaction length for each energy, for e+/e- production:
    intlen_pp = [get_intlen(x, get_Neps_CMB) for x in E]
    intlen_pp = np.array(intlen_pp)
    intlen_pp/=Mpc_in_m

    # calculate mean interaction length for each energy, for mu+/mu- production:
    intlen_mu = [get_intlen(x, get_Neps_CMB, use_pp=False, use_mu=True) for x in E]
    intlen_mu = np.array(intlen_mu)
    intlen_mu/=Mpc_in_m

    # plot results:
    pl.subplot(111)
    pl.plot(E/q_e, intlen_pp)
    pl.plot(E/q_e, intlen_mu)
    pl.xlabel('Photon energy (eV)')
    pl.ylabel('Interaction length (Mpc)')
    pl.loglog()
    ymin,ymax = pl.ylim()
    pl.ylim(ymin,1e7*ymin)
    pl.show()
