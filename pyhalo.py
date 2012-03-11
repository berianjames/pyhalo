"""pyhalo.py

Python implementation of routines to compute a variety
of geometric statistics using the halo paradigm."""

def ps_to_xi(k,ps,r,precision='mid'):
    xi = 0 * r
    if precision=='low':
        from scipy.integrate import trapz
        from scipy import sin
        
        for i in range(len(r)):
            xi[i] = trapz((ps/k) * sin(k*r[i]) / (k*r[i]),k)
    
    elif precision=='high':
        from scipy import sin
        from scipy.integrate import quad,quadrature
        from scipy.interpolate import interp1d
        for i in range(len(r)):
            psi = interp1d(k,ps)
            core = lambda k: (psi(k)/k) * sin(k*r[i]) / (k*r[i])
            A = quadrature(core,min(k),2.0/r[i],tol=1e-3)
            xi[i] = A[0]
        
    else:
        from scipy.integrate import trapz
        from scipy import sin
        from pylab import find
        import numpy as np
        for i in range(len(r)):
            cutoff = np.min(find(k>2.0/r[i]))
            kl = k[1:cutoff]
            psl = ps[1:cutoff]
            i1 = trapz((psl/kl) * sin(kl*r[i]) / (kl*r[i]),kl)
            psh = ps[cutoff+1:len(k)]
            kh = k[cutoff+1:len(k)]
            i2 = trapz((psh/kh) * sin(kh*r[i]) / (kh*r[i]),kh)
            xi[i] = i1+i2
                         
    return xi

#def ps_to_xi(ks,ps,r,precision='lower'):
#    from scipy.integrate import quad
#    from math import sin
#    xi = 0*r
#    for i in range(len(r)):
#        core = lambda k: (power_spectrum(k,z,C) / k) * (sin(k*r[i])/(k*r[i]))
#        A = quad(core,a,b)
#        xi[i] = A[0]
#    return xi                                      

def power_spectrum(k,z,C,precision='lower'):
    if precision=='higher':
        # Calls CMBFast -- Not yet implemented
        print 'Not yet implemented'
    
    # Uses Eisenstein & Hu (1999) transfer function fitting
    g = growth(z,C,precision) / growth(0,C,precision)
    T = 0*k + transfer(k,C)          
    D2m = 0*k + ps_norm(C) * (T*g)**2 \
        * (k / C.k_WMAP)**(C.n_s-1) * (k * C.cH0)**4
    
    return D2m

def ps_norm(C):
    norm = 4./25. * C.D2_R
    return norm

def transfer(ks,C):
    """Transfer function for the baryon/dark matter fluid, as
    per Eisenstien & Hu (1998) ApJ 496 605. Comments reference
    equations in this work."""
    from math import log,sqrt,e,sin
    
    # General definitions
    om0h2 = C.om_0 * C.h**2   # Omega_m * h^2
   
    # Baryon-to-photon momentum ratio at redshift z (Eq. 5)
    R = lambda z: 31.5*C.om_b*C.h**2*C.thet**-4*(z/1e3)**-1

    # Compton drag epoch (Eq. 4)
    b1 = 0.313*om0h2**-0.419 * (1 + 0.607*om0h2**0.674)
    b2 = 0.238*om0h2**0.223
    zd = 1291*om0h2**0.251/(1+0.659*om0h2**0.828)*(1+b1*(C.om_b*C.h**2)**b2)

    # Sound horizon at drag epoch (Eqs 5,6)
    Rd = R(zd)
    Req = R(C.z_eq)
    s = (2 / (3*C.k_eq)) * sqrt(6 / Req) * \
        log( (sqrt(1+Rd) + sqrt(Rd+Req)) / (1 + sqrt(Req)) )

    # Silk damping scale (Eq. 7)
    k_Silk = 1.6*(C.om_b*C.h**2)**0.52 * om0h2**0.73 * \
        (1 + (10.4*om0h2)**-0.95)

    # CDM transfer component parameters (Eqs. 9 - 12)
    b1 = 0.944*(1 + (458*om0h2)**-0.708)**-1
    b2 = (0.395*om0h2)**-0.0266
    beta_c = (1 + b1*((C.om_c/C.om_0)**b2 - 1))**-1.

    a1 = (46.9*om0h2)**0.670*(1+(32.1*om0h2)**-0.532)
    a2 = (12.0*om0h2)**0.424*(1+(45.0*om0h2)**-0.582)
    alpha_c = a1**-(C.om_b/C.om_0) * a2**-(C.om_b/C.om_0)**3

    # Baryon transfer component parameters (Eqs. 14,15,21 - 24)
    G = lambda y: y*( -6*sqrt(1+y) + (2+3*y)*log( (sqrt(1+y)+1)/(sqrt(1+y)-1) ) )
    alpha_b = 2.07*C.k_eq*s*(1+Rd)**(-0.75) * G( (1+C.z_eq) / (1+zd) )
    beta_b = 0.5 + C.om_b/C.om_0 + (3-2*C.om_b/C.om_0) * sqrt( (17.2*om0h2)**2+1 )
    beta_node = 8.41*om0h2**0.435

    T = 0 * ks
    for i in range(len(ks)):
        # Effective wavenumber (Eq. 3)
        k = ks[i] * C.h
        q = k / (13.41*C.k_eq)    
    
        # Transition switch (Eq. 18)
        f = 1.0 / (1 + (k*s/5.4)**4)

        # Grand interpolation function (Eq. 19)
        def T_interp(q,alpha_c,beta_c):
            C = 14.2/alpha_c + 386.0/(1+69.9*q**1.08)
            T_0 = log(e+1.8*beta_c*q) / (log(e+1.8*beta_c*q) + C*q**2)
            return T_0
        
        T_c = f*T_interp(q,1,beta_c) + (1-f)*T_interp(q,alpha_c,beta_c)

        stilde = s / (1 + (beta_node/(k*s))**3)**(1.0/3.0)

        T_b = (T_interp(q,1,1) / (1 + (k*s/5.2)**2) + \
                   alpha_b / (1 + (beta_b/(k*s))**3) * e**(-(k/k_Silk)**1.4) ) * \
                   sin(k * stilde) / (k * stilde)

        # Transfer function for baryons + CDM (Eq. 16)
        T[i] = C.om_b/C.om_0 * T_b + C.om_c/C.om_0 * T_c

    return T
    

def growth(z,C,precision='lower'):
    """Growth function, takes redshift and cosmology as input."""
    if precision=='higher':
        # Not yet implemented
        print 'Not yet implemented'

    """Approximation from Carroll, Press & Turner (1992)."""
    g2 = C.om_0*(1+z)**3 + (1 - C.om_0 - C.om_v)*(1+z)**2 + C.om_v
    om = C.om_0*(1+z)**3 / g2
    om_L = C.om_v / g2
    D1 = (5.*om/2.) / (om**(4./7.) - om_L + (1+om/2.)*(1+om_L/70.)) / (1+z)
    
    return D1 

class Cosmology:
    """Defines the cosmological parameter set and
    associated quantities. This is only for the
    Lambda CDM cosmology, using WMAP-7 figures."""
    def __init__(self):
        from math import sqrt

        # Fundamental parameters
        self.h = 0.702     # Hubble parameter
        self.om_c = 0.227  # CDM density
        self.om_b = 0.0455 # Baryonic matter density
        self.om_m = 0.272  # Total matter density
        self.om_v = 0.728  # Dark energy density
        self.n_s = 0.961   # Spectral index of density perturbations
        self.D2_R = 2.45e-9 # Amplitude of curvature perturbation at k_WMAP
        
        # Auxiliary parameters
        self.k_WMAP = 0.002*self.h # Wavenumber of normalisation
        self.sig8 = 0.807     # Amplitude of linear density fluctuation at 8 Mpc/h
        #self.z_eq = 3196      # Redshift of matter-radiation equality
        
        # Derived parameters
        self.om_0 = self.om_c + self.om_b
        self.f_c = self.om_c / self.om_0
        self.f_b = self.om_b / self.om_0
        self.f_cb = (self.om_c + self.om_b) / self.om_0
        self.H0 = 100*self.h
        self.cH0 = 299792.458 / self.H0   # c/H_0 in Mpc/h
        self.thet = 2.728 / 2.7 
        self.k_eq = 7.46e-2 * self.om_0 * self.h**2 * self.thet**-2
        self.z_eq = 2.5e4 * self.om_0 *self.h**2 * self.thet**-4.

    # Function to update auxiliary and derived parameters
    def update(C):
        # Auxiliary parameters
        C.k_WMAP = 0.002*C.h # Wavenumber of normalisation
        
        # Derived parameters
        C.om_0 = C.om_c + C.om_b
        C.f_c = C.om_c / C.om_0
        C.f_b = C.om_b / C.om_0
        C.f_cb = (C.om_c + C.om_b) / C.om_0
        C.H0 = 100*C.h
        C.cH0 = 299792.458 / C.H0   # c/H_0 in Mpc/h
        C.thet = 2.728 / 2.7 
        C.k_eq = 7.46e-2 * C.om_0 * C.h**2 * C.thet**-2
        C.z_eq = 2.5e4 * C.om_0 *C.h**2 * C.thet**-4.

        return C
        



        
        
