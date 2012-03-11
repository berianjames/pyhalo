from pylab import *
from pyhalo import Cosmology, power_spectrum, transfer, ps_to_xi
from scipy import pi,sin
from time import sleep

close('all')

# Parameters
C = Cosmology()
z = 0
k = logspace(-3,3,1000)
T = transfer(k,C)
loglog(k,T)
#ps = power_spectrum(k,z,C)
#r = logspace(-2,2,100)
#xi = ps_to_xi(k,ps,r,'high')
#loglog(r,xi,k,ps)
#xlim(1e-2,1e2)
#ylim(1e-4,1e2)

#fname = '/Users/berian/Research/Tools/cmbfast/berian/lcdm_berian.dat'
#data = load(fname)
#kd = data[:,0]
#P = data[:,1]
#P = P / max(P)
#Td = interp(kd,k,T)
#semilogx(kd,Td/P)
#axis([0.001,1,0.95,1.05])

    

