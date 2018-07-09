from scipy.integrate import quad
import numpy as np

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

def logquad(func, a, b, **kwargs):

  """
  Perform scipy.integrate.quad more efficiently on a function that's better-behaved in log space.
  """

  # Perform integral $\int y dx$ by transforming variables x->log(x), y->log(y).

  return quad(lambda x,*args: np.exp(x + np.log(func(np.exp(x),*args))), np.log(a), np.log(b), **kwargs)
