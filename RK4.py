# Advances the solution of a set of differential equations by a small time step
# using the fourth-order Runge-Kutta algorithm.
#
# deriv:  function returning an array that is derivatives of x[0], x[1], etc.
# x:  array of current values of x[0], x[1], etc.
# h:  time step
#
# returns an array of new values of x[0], x[1], etc.
#
def RK4(deriv,x,h):
   k1 = h*deriv(x)
   x2 = x + 0.5*k1
   k2 = h*deriv(x2)
   x3 = x + 0.5*k2
   k3 = h*deriv(x3)
   x4 = x + k3
   k4 = h*deriv(x4)
   return x +(k1 +2.0*(k2+k3) + k4)/6.0


