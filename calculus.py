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
      The Richardson extrapolation algorithm for evaluating 
      the derivative of f(x) at thepoint x. N_max is the 
      maximum number of iterations for the algorithm.
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
    
    
class integrate:

  def romberg(f, ab, n=1, N_max=5):
      """
      The Romberg integration algorithm for evaluating the integral f(x) over the
      interval ab, where ab is specified as (for example) ab=[0,1]. n is the number of 
      subintervals to start with (default n=1) and N_max is the maximum number of iterations
      for the algorithm
      """
      R = np.zeros((N_max, N_max))  # define a matrix which will contain the approximations
      R[0,0] = trapezoid(f, ab, n)  # initiate the first approximation
      for i in range(N_max-1):
          R[i+1,0] = trapezoid(f, ab, (2**(i+1))*n) # fill in the remaining approximations for increasing number of subintervals
          for j in range(i+1):
              R[i+1,j+1] = (4**(j+1) * R[i+1,j] - R[i,j])/(4**(j+1) - 1) # iterate over the previous approximations
                                                                       # to eliminate successive errors
      return R[N_max-1,N_max-1]




  def trapezoid(f, ab, n):
      '''
      The n-segment trapezoid rule function for evaluating the integral of f(x).
      ab is the interval which we want to integrate over (eg. ab=[0,1])
      n is the number of subintervals 
      '''
      sub_intvls = divide_intvl(ab[0],ab[1],n) # divide the 
      I = 0 # initiate the value of the integral calculation
      for intvl in sub_intvls:
          I = I + (intvl[1]-intvl[0])/2*(f(intvl[0])+f(intvl[1])) # add the evaluation 
                                                                # of the trapezoid rule
                                                                # for each subinterval
      return I


  def divide_intvl(xL, xR, n):
      """
      Divides the interval (xL, xR) into n subintervals.
      Returns an array of each n subinterval (also stored as an array)\n
      eg. divide(-1,1,2) returns [[-1,0],[0,1]]
      """
      width = (xR - xL)/n
      intvl = []
      for i in range(n):
          intvl.append([xL, xL+width])
          xL = xL+width
return intvl
