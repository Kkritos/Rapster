'''
 Copyright (C) 2026  Konstantinos Kritos <kkritos1@jhu.edu>

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.

'''

from .constants import *

def E_cosmo(z):
    """
    @in z: redshift
    
    @out: auxiliary cosmological function
    """
    
    return np.sqrt(Omega_R * (1 + z)**4 + Omega_M * (1 + z)**3 + Omega_K * (1 + z)**2 + Omega_V)

def lookback(z):
    """
    @in z: redshift
    
    @out: lookback time [Myr]
    """

    tL = 0.0
    zz = 0.0
    dzz = 1e-3
    while zz < z:
        tL += dzz / (1 + zz) / E_cosmo(zz) * t_Hubble
        zz += dzz
    
    return tL

def redshift(t):
    """
    @in t: lookback time [Myr]
    
    @out: redshift
    """
    
    z = 0.0
    
    dz = 1e-3

    while lookback(z) < t:
        z += dz
    
    return z

def T_GW(m1, m2, a0, e0):
    """
    I. Mandel (2021) fit to Peters timescale.
    
    @in m1: primary mass of BBH [Msun]
    @in m2: secondary mass of BBH [Msun]
    @in a0: initial semimajor axis of BBH [pc]
    @in e0: initial eccentricity of BBH
    
    @out: GW coalescence timescale [Myr]
    """
    
    # Coalescence time for circular orbits:
    Tc = 5 * c_light**5 * a0**4 / 256 / G_Newton**3 / m1 / m2 / (m1 + m2)
    
    # Eccentricity factor:
    factor_e = (1 + 0.27 * e0**10 + 0.33 * e0**20 + 0.2 * e0**1000) * (1 - e0**2)**(7/2)
    
    return Tc * factor_e

def t_relax(Mcl, rh, m_avg, psi, logL):
    """
    @in Mcl: cluster's mass [Msun]
    @in rh: half-mass radius [pc]
    @in m_avg: average mass [Msun]
    @in psi: multimass relaxation factor
    @in logL: Coulomb logarithm
    
    @out: half-mass relaxation timescale [Myr]
    """
    
    return 0.138 * np.sqrt(Mcl * rh**3 / G_Newton) / m_avg / psi / logL

def v_esc(Mcl, rh):
    """
    @in Mcl: cluster mass [Msun]
    @in rh: half-mass radius [pc]
    
    @out: escape velocity [km/s]
    """
    
    return 2 * np.sqrt(0.4 * G_Newton * Mcl / rh)

def sample_angles():
    """
    @out theta1: [rad]
    @out theta2: [rad]
    @out dPhi: [rad]
    """
    
    cos_theta1 = np.random.uniform(-1, 1)
    cos_theta2 = np.random.uniform(-1, 1)
    
    theta1 = np.arccos(cos_theta1)
    theta2 = np.arccos(cos_theta2)
    
    dPhi = np.random.uniform(0, 2*np.pi)
    
    return theta1, theta2, dPhi

def Rate_3bb(m, n, v):
    """
    @in m: mass scale [Msun]
    @in n: number density [pc^-3]
    @in v: velocity dispersion [km/s]
    
    @out: three-body binary formation rate [1/Myr]
    """
    
    return 8 * np.pi / np.sqrt(3) * n**3 * (G_Newton * m)**5 / v**9 * eta_min**(-11/2) * (1 + 3 * eta_min) * (1 + 6 * eta_min) * P_3bb

def Rate_cap(m, n, v):
    """
    @in m: mass [Msun]
    @in n: number density [pc^-3]
    @in v: velocity dispersion [km/s]
    
    @out: 2-body capture rate [1/Myr]
    """

    # relative velocity:
    v_rel = np.sqrt(2) * v
    
    # capture cross section:
    Sigma_cap = 2 * np.pi * (85 * np.pi / 6 / np.sqrt(2))**(2/7) * G_Newton**2 * (2 * m)**(10/7) * m**(4/7) / c_light**(10/7) / v_rel**(18/7)
    
    return n**2 * Sigma_cap * v_rel

def Rate_int(m, n, v, rp):
    """
    @in m: total mass [Msun]
    @in n: number density [pc^-3]
    @in v: velocity dispersion [km/s]
    @in rp: pericenter of interaction [pc]
    
    @out: interaction rate [1/Myr]
    """
    
    # relative velocity:
    v_rel = np.sqrt(2) * v
    
    # interaction cross section:
    Sigma_int = 2 * np.pi * G_Newton * m * rp / v_rel**2
    
    return n * Sigma_int * v_rel

def Rate_exc(m1, m2, m3, n3, v_rel, a):
    """
    @in m1: to-be-exchanged mass of binary
    @in m2: primary mass of binary
    @in m3: incoming mass
    @in n3: number density of incoming mass
    @in v_rel: binary-single relative velocity before interaction
    @in a: initial binary semimajor axis
    
    @out: 1-2 -> 3-2 exchange rate [1/Myr]
    """

    m12 = m1 + m2
    m123 = m12 + m3
    m13 = m1 + m3
    m23 = m2 + m3
    mu1 = m1 / m12
    mu2 = m3 / m123
    
    # exchange cross section [Heggie et al. (1996)]:
    Sigma_exc = 6.76e-5 * a * (10 / v_rel)**2 * m123 * (m23/m123)**(1/6) * (m3/m13)**(7/2) * (m123/m12)**(1/3) * (m13/m123) \
    * np.exp(3.70 + 7.49*mu1 - 1.89*mu2 - 15.49*mu1**2 - 2.93*mu1*mu2 - 2.92*mu2**2 + 3.07*mu1**3 + 13.15*mu1**2*mu2 - 5.23*mu1*mu2**2 + 3.12*mu2**3)
    
    return n3 * Sigma_exc * v_rel
    
def sample_hardness():
    """
    @out: 3bb hardness parameter
    """
    
    return eta_min * (1 - np.random.rand())**(-2/7)

def CDF_maxwell(w):
    """
    Returns cumulative density function of Maxwellian distribution.

    @in w: velocity value, normalized to the 1D velocity dispersion
    """

    return erf(w/np.sqrt(2)) - np.sqrt(2/np.pi)*w*np.exp(-w**2/2)

w_s = np.linspace(0, 10000, 10**6) # normalized velocities (assuming 1D velo. disp. =1)
CDF_maxwell_s = np.vectorize(CDF_maxwell)(w_s)

def get_maxwell_sample(sigma):
    """
    Returns a sample from the Maxwellian with the inverse sampling method.

    @in sigma: 1D velocity dispersion parameter
    """
    
    return sigma * np.interp(np.random.rand(), CDF_maxwell_s, w_s)

def DW(x, m0, m1, m2, a1, a2, e1, e2, cosi1, cosi2, w):
    """
    Change in the double averaged Hamiltonian of a hierarchical triple as inner eccentricity is varied.
    Notation from M.C.Miller & D.P.Hamilton ApJ 576, 894 (2002).
    
    @in x: 1-e1^2 (see definition of e1 below)
    @in m0: primary mass of inner binary
    @in m1: secondary mass of inner binary
    @in m2: tertiary mass
    @in a1: inner binary semimajor axis
    @in a2: outer binary semimajor axis
    @in e1: inner binary eccentricity
    @in e2: outer binary eccentricity
    @in cosi1: cosine of inner binary's inclination relative to invariant plane
    @in cosi2:    "   of outer    "          "          "     "     "       "
    @in w: argument of periapsis for inner binary
    """
 
    E = 1 - e1**2

    # inclinations:
    i1, i2 = np.arccos(cosi1), np.arccos(cosi2)
 
    sinw = np.sin(w)
    
    M1 = m0 + m1 # total inner binary mass
    M2 = M1 + m2 # total outer binary mass
    
    mu1 = m0*m1/M1 # reduced mass of inner binary
    mu2 = M1*m2/M2 # reduced mass of outer binary

    # integrals of motion:
    B = mu2/mu1*np.sqrt(M1*a1/M2/a2*(1-e2**2))
    A = np.sqrt(E)*np.cos(i1) + B*np.cos(i2)

    # general relativistic correction term at zero inner eccentricity:
    tPN = 8*M1/m2*(a2/a1*np.sqrt(1-e2**2))**3*M1/232/a1/(3e5)**2
    
    cosi = (A**2-B**2-E)/2/B/np.sqrt(E)
    W = - 2*E + E*cosi**2 + 5*(1-E)*sinw**2*(cosi**2-1) + tPN/np.sqrt(E)
    
    cosI = (A**2-B**2-x)/2/B/np.sqrt(x)
    
    return -2*E + E*cosI**2 + 5*(1-E)*(cosI**2-1) + tPN/np.sqrt(x) - W

def find_eMAX(m0, m1, m2, a1, a2, e1, e2, cosi1, cosi2, w, x0=1e-10, tol=1e-10):
    """
    Maximum eccentricity during a von Zeipel-Lidov-Kozai cycle.
    Same inputs as function DW in this file.
    """
 
    x_min = root(fun=DW, x0=x0, args=(m0, m1, m2, a1, a2, e1, e2, cosi1, cosi2, w), tol=tol)['x'][0]
    
    eMAX = np.sqrt(1 - x_min)
    
    return eMAX

def X_isco(x, s):
    """
    Returns the isco in units of gravitational radius (G*m/c^2).
    Inputs:
    x: dimensionless spin parameter in [0, 1]
    s: 1 for prograde and -1 for retrograde orbits.
    """
    
    Z1 = 1 + (1 - x**2)**(1/3)*((1+x)**(1/3) + (1-x)**(1/3))
    Z2 = np.sqrt(3*x**2 + Z1**2)
    XX = 3 + Z2 - s*np.sqrt((3-Z1)*(3+Z1+2*Z2))
    return XX

def dxdM(M, x, s):
    """
    Bardeen formula; Bardeen 1970 ("Kerr Metric Black Holes")
    M: mass
    x: dimensionless spin parameter in [0, 1]
    """
    
    return 2/3/np.sqrt(3)*s/M*(1 + 2*np.sqrt(3*X_isco(x, s) - 2))/np.sqrt(1 - 2/3/X_isco(x, s)) - 2*x/M

def evolve_spin_RungeKutta(Mi, Mf, xi, s, dM):
    """
    Returns the final spin of an accretion episode. Uses Runge-Kutta's 4th-order method.
    Mi: initial compact object mass
    Mf: final compact object mass so that (Mf-Mi)=accreted mass
    xi: initial compact object dimensionless spin parameter in [0, 1]
    s: +1 for prograde and -1 for retrograde accretion
    dM: numerical integration mass-step
    """
    
    M = Mi
    x = xi
    while M<Mf:
        k1 = dM*dxdM(M, x, s)
        k2 = dM*dxdM(M+dM/2, x+k1/2, s)
        k3 = dM*dxdM(M+dM/2, x+k2/2, s)
        k4 = dM*dxdM(M+dM, x+k3, s)
        dx = (k1 + 2*k2 + 2*k3 + k4)/6
        if x + dx > 1:
            x = 1
            break
        elif x + dx < 0:
            s = -s
        else:
            x = x + dx
            M = M + dM
    xf = x
    return xf

# end of file
