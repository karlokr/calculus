import numpy as np

class derive:

  def sevenpoint(f, x, h=0.000855):
      """
      Uses the 7-point central-divided difference 
      to approximate the value of df/dx where the default
      value for h is h=0.000855
      """
      return f(x+3*h)/(60*h) + f(x+2*h)*(-3)/(20*h) + f(x+h)*3/(4*h) + f(x-h)*3/(-4*h) + f(x-2*h)*3/(20*h) + f(x-3*h)/(-60*h) 

  def richardson(f, x, h=1e-3, N_max=3):
      """
      The Richardson extrapolation algorithm for evaluating the derivative of f(x) at the
      point x, where ab is specified as (for example) ab=[0,1]. N_max is the maximum number 
      of iterations for the algorithm.
      """
      R = np.zeros((N_max, N_max))
      R[0,0] = cdd(f, x, h)
      for i in range(N_max-1):
          R[i+1,0] = cdd(f, x, h/(2**(i+1)))
          for j in range(i+1):
              R[i+1,j+1] = (4**(j+1) * R[i+1,j] - R[i,j])/(4**(j+1) - 1)
      return R[N_max-1,N_max-1]


  def cdd(f, x, h=1e-8):
      """
      Uses the central-divided difference (f(x+h)-f(x-h))/(2 * h) 
      to approximate the value of df/dx where the default
      value for h is h=1e-8
      """
      return (f(x+h) - f(x-h))/(2*h)
