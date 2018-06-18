# Calculate gamma-ray interaction length for photon-photon pair production on the cosmic microwave background.
# [180607 - Justin Bray] Original code
# [180608 - Anna Scaife] Updated to include scipy integration [not always successfully]



import pylab as pl
import numpy as np
import scipy.integrate as integrate
from scipy.integrate import quad

# ------------------------------------------------------------------------------

# Physical constants.
m_e = 9.11e-31 # kg; electron mass
q_e = 1.67e-19 # C; electron charge
c = 3e8 # m/s; speed of light
sigma_t = 6.65e-29 # m^2; Thomson cross-section
T_CMB = 2.725 # K; temperature of cosmic microwave background
h = 6.626e-34 # Js; Planck constant
k = 1.38e-23 # J/K; Boltzmann constant
Mpc_in_m = 3.086e22 # m; 1 Mpc

# ------------------------------------------------------------------------------

def get_sigma_pp(s):

  """
  Calculate p-p pair production cross-section, per berestetskii1982,
  as cited in eqns 1&2 of
  https://www.astro.umd.edu/~miller/teaching/astr688m/lecture03.pdf
  """

  v = np.sqrt(1. - (m_e*c**2)**2 / (s/4.))
  sigma_pp = 3./16. * sigma_t * (1. - v**2) * ( (3.-v**4)\
            * np.log( (1.+v)/(1.-v) ) - 2.*v*(2.-v**2)) # m^2; p-p cross-section

  return sigma_pp

# ------------------------------------------------------------------------------

def s_integrand(s):

  """
  Define the final integrand, over s, from eqn 23 of protheroe1996.
  """

  integrand = s * get_sigma_pp(s) # J^2 m^2

  return integrand

# ------------------------------------------------------------------------------

def s_integral(E, eps):

  """
  Calculate the final integral, over s, from eqn 23 of protheroe1996.
  """

  # protheroe1996, eqn 26:
  smin = (2.*m_e*c**2)**2
  smin*= 1. + 1e-12 # small fudge factor to avoid infinities

  # J^2; maximum CoM energy^2; protheroe1996, eqn 28:
  smax = 4*E*eps

  # check limits:
  assert np.log10(smax/smin) < 10, 'Probably too coarse in the smin->smax integral.'
  if smax < smin: return 0.

  # numerical integration for checking:
  #s = np.exp( np.linspace(np.log(smin), np.log(smax), 500) )
  #integrand = s * get_sigma_pp(s)
  #integral = np.sum(0.5*(integrand[:-1]+integrand[1:])*np.diff(s))

  # protheroe1996, eqn 23:
  integral = quad(s_integrand, smin, smax, limit=500)[0]

  return integral

# ------------------------------------------------------------------------------

def eps_integrand(eps, E):

  """
  Define the integrand, over epsilon, from eqn 23 of protheroe1996.
  """

  # get the differential photon number density, protheroe1996, eqn 24:
  Neps = get_Neps(eps)

  integrand = Neps * eps**-2 * s_integral(E, eps)

  return integrand

# ------------------------------------------------------------------------------

def eps_integral(E):

  """
  Calculate the final integral, over epsilon, from eqn 23 of protheroe1996.
  """

  # protheroe1996, eqn 27:
  epsmin = (2*m_e*c**2)**2 / (4*E) # Joules
  epsmax = 1e10*epsmin # J; target photon energy above majority of effect
  #epsmax = np.inf

  # protheroe1996, eqn 23:
  eps = np.exp( np.linspace(np.log(epsmin), np.log(epsmax), 500) )
  Neps = get_Neps(eps)

  integrand = Neps * eps**-2 * np.array([s_integral(E, x) for x in eps])
  integral = np.sum(0.5*(integrand[:-1]+integrand[1:])*np.diff(eps))

  # protheroe1996, eqn 23:
  # scipy integration doesn't work because integrand is very highly localised.
  #integral = quad(eps_integrand, epsmin, epsmax, args=(E), limit=500)[0]

  return integral

# ------------------------------------------------------------------------------

def get_intlen(E):

  """
  Find interaction length for photons of given energy,
  from eqn 23 of protheroe1996.
  """

  return 8*E**2 / eps_integral(E) # m

# ------------------------------------------------------------------------------

def get_Neps(eps):

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

def get_Neps_muon(eps):

  """
  Find differential photon number density from muon pair production
  in the microwave background at a given photon energy.
  """

  nu = eps / h # Hz; photon frequency

  return

 # ------------------------------------------------------------------------------

 def get_Neps_4e(eps):

   """
   Find differential photon number density from double pair production
   in the microwave background at a given photon energy.
   """

   nu = eps / h # Hz; photon frequency

   return


# ------------------------------------------------------------------------------

def get_Neps_IR(eps):

  """
  Find differential photon number density from the optical and infrared
  background at a given photon energy.
  """

  nu = eps / h # Hz; photon frequency

  return

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

if __name__ == "__main__":

    # define array of energies:
    E = 10**np.linspace(13,22,20) # eV; photon energy
    E *= q_e # eV->J

    # calculate mean interaction length for each energy:
    intlen = [get_intlen(x) for x in E]
    intlen = np.array(intlen)
    intlen/=Mpc_in_m

    # plot results:
    pl.subplot(111)
    pl.plot(E/q_e, intlen)
    pl.xlabel('Photon energy (eV)')
    pl.ylabel('Interaction length (Mpc)')
    pl.loglog()
    ymin,ymax = pl.ylim()
    pl.ylim(ymin,1e5*ymin)
    pl.show()
