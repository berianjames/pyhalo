from pylab import *
from pyhalo import Cosmology, power_spectrum, transfer, ps_to_xi
from scipy import pi,sin
from time import sleep

close('all')

# Parameters
C = Cosmology()
z = 0.3
k = logspace(-1,1,100)
#r = logspace(-2,2,100)

ps = power_spectrum(k,z,C)

#xi = ps_to_xi(k,ps,r,'high')

#loglog(r,xi,color='b')
#loglog(k,ps2,color=[1 0 0])
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

    

