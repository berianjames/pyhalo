from pylab import logspace, loglog, log10
from pyhalo import Cosmology
from math import pi,sin
from scipy.integrate import quad

# Parameters
Delta_c = 180
M = 1e13
c = 1

C = Cosmology()
r_h = ( 3*M / (4*pi*2.78e11*C.om_0) / Delta_c ) ** (1./3.)

ks = logspace(-2,2,100)
y = 0 * ks
for i in range(len(ks)):
    k = ks[i]
    core = lambda r: 4*pi*r**2 / ( r*c/r_h * (1+r*c/r_h)**2 ) * sin(k*r) / (k*r)
    x = quad(core,0,10*r_h)
    y[i] = x[0] / M

y = y / y[0]
loglog(ks,y)





    

