'''
 Copyright (C) 2023  Konstantinos Kritos <kkritos1@jhu.edu>

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

from constants2 import *

def lookback_astropy(z):
    """
    Assuming Planck (2018) cosmological parameters. Using astropy.
    
    @in z: redshift
    
    @out: lookback time [Myr]
    """
    
    return Planck18.lookback_time(z).value * 1e3

def redshift_astropy(t):
    """
    Assuming Planck (2018) cosmological parameters. Using astropy.
    
    @in t: lookback time [Myr]
    
    @out: redshift
    """
    
    return z_at_value(Planck18.lookback_time, t * u.Myr)

@njit
def E_cosmo(z):
    """
    @in z: redshift
    
    @out: auxiliary cosmological function
    """
    
    return np.sqrt(Omega_R * (1 + z)**4 + Omega_M * (1 + z)**3 + Omega_K * (1 + z)**2 + Omega_V)

@njit
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

@njit
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

def IMF_kroupa(m):
    """
    Kroupa (2002) initial mass function.

    @in m : stellar mass array [Msun]

    @out: number dN of stars in mass bin (m,m+dm) [1/Msun]
    """

    # mass boundaries (in solar masses):
    m1 = 0.08
    m2 = 0.50
    m3 = 1.00

    # spectral indices (broken power law; central values):
    a0 = -0.3
    a1 = -1.3
    a2 = -2.3
    a3 = alphaIMF

    # normalization constants:
    c1 = m1**a0 / m1**a1
    c2 = c1 * m2**a1 / m2**a2
    c3 = c2 * m3**a2 / m3**a3

    out = np.zeros(m.size)

    for i in range(0,m.size):
        
        if  (m[i] <= m1):
            
            out[i] = m[i]**a0
            
        elif(m[i] <= m2 and m[i] > m1):
            
            out[i] = c1 * m[i]**a1
            
        elif(m[i] <= m3 and m[i] > m2):
            
            out[i] = c2 * m[i]**a2
            
        elif(m[i] >= m3):
            
            out[i] = c3 * m[i]**a3
                    
    return out

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

def merger_remnant(m1, m2, chi1, chi2, theta1, theta2, dPhi):
    """
    Final mass, final spin parameter and GW kick velocity of a merger remnant
    calculated with the PRECESSION package.
    Ref: D.Gerosa & M.Kesden, PRD 93 (2016), 124066.
    
    @in  m1: first  BH mass [<units>]
    @in  m2: second BH mass [<units>]
    @in  chi1: first  BH spin parameter in [0,1]
    @in  chi2: second BH spin parameter in [0,1]
    @in  theta1: angle of spin 1 with angular momentum [rad]
    @in  theta2: angle of spin 2 with angular momentum [rad]
    @in  dPhi: angle between the orbital plane projections of the spins [rad]
    
    @out mRem    : merger remnant BH mass [<units>]
    @out chiRem  : merger remnant BH spin parameter in [0,1]
    @out vGWkick : merger remnant GW kick velocity [km/s]
    """
    
    # BBH mass ratio:
    q = np.min([m1, m2]) / np.max([m1, m2])
    
    chi_p   = chi1
    chi_s   = chi2
    
    theta_p = theta1
    theta_s = theta2
    
    # order `1`-> primary, `2`-> secondary:
    if m2 > m1:
        
        chi_p   = chi2
        chi_s   = chi1
        
        theta_p = theta2
        theta_s = theta1
        
    M, m_1, m_2, S_1, S_2 = pre.get_fixed(q, chi_p, chi_s) # units M=1
    
    # Final mass of merger remnant in solar masses:
    m_rem = pre.finalmass(theta_p, theta_s, dPhi, q, S_1, S_2) * (m1 + m2)
    
    # Final spin of merger remnant:
    chi_rem = pre.finalspin(theta_p, theta_s, dPhi, q, S_1, S_2)
    
    # Final GW kick:
    vGW_kick = pre.finalkick(theta_p, theta_s, dPhi, q, S_1, S_2, maxkick=False, kms=True, more=False)
    
    return m_rem, chi_rem, vGW_kick

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
    
N_grid = 500
M_grid = np.linspace(20, 340, N_grid)
Z_grid = np.logspace(np.log10(1e-4), np.log10(2e-2), N_grid)

Mremnants_F12d = np.loadtxt('./MzamsMrem/MzamsMrem_F12d.txt', unpack=True)

Mremnants_F12d = np.transpose(Mremnants_F12d)

MremInterpol_F12d = interpolate.interp2d(M_grid, Z_grid, Mremnants_F12d, kind='linear', bounds_error=True)

def Mrem_F12d(M, Z):
    """
    Fryer et al. (2002) delayed remnant mass prescription model.
    
    @in M: ZAMS mass [Msun]
    @in Z: metallicity

    @out : remnant mass [Msun]
    """

    return MremInterpol_F12d(M, Z)

# Reading files exported from SEVN code and stored according to metallicity for various ZAMS masses:
MzamsMrem1  = np.load('./MzamsMrem/MzamsMrem1.npz' ); Mrem_delayed_1  = MzamsMrem1 ['Mrem1' ]
MzamsMrem2  = np.load('./MzamsMrem/MzamsMrem2.npz' ); Mrem_delayed_2  = MzamsMrem2 ['Mrem2' ]
MzamsMrem3  = np.load('./MzamsMrem/MzamsMrem3.npz' ); Mrem_delayed_3  = MzamsMrem3 ['Mrem3' ]
MzamsMrem4  = np.load('./MzamsMrem/MzamsMrem4.npz' ); Mrem_delayed_4  = MzamsMrem4 ['Mrem4' ]
MzamsMrem5  = np.load('./MzamsMrem/MzamsMrem5.npz' ); Mrem_delayed_5  = MzamsMrem5 ['Mrem5' ]
MzamsMrem6  = np.load('./MzamsMrem/MzamsMrem6.npz' ); Mrem_delayed_6  = MzamsMrem6 ['Mrem6' ]
MzamsMrem7  = np.load('./MzamsMrem/MzamsMrem7.npz' ); Mrem_delayed_7  = MzamsMrem7 ['Mrem7' ]
MzamsMrem8  = np.load('./MzamsMrem/MzamsMrem8.npz' ); Mrem_delayed_8  = MzamsMrem8 ['Mrem8' ]
MzamsMrem9  = np.load('./MzamsMrem/MzamsMrem9.npz' ); Mrem_delayed_9  = MzamsMrem9 ['Mrem9' ]
MzamsMrem10 = np.load('./MzamsMrem/MzamsMrem10.npz'); Mrem_delayed_10 = MzamsMrem10['Mrem10']
MzamsMrem11 = np.load('./MzamsMrem/MzamsMrem11.npz'); Mrem_delayed_11 = MzamsMrem11['Mrem11']
MzamsMrem12 = np.load('./MzamsMrem/MzamsMrem12.npz'); Mrem_delayed_12 = MzamsMrem12['Mrem12']

# collect remnant masses with various metallicity values in a single array:
Mrem_delayed = np.array([Mrem_delayed_1, Mrem_delayed_2, Mrem_delayed_3, Mrem_delayed_4, Mrem_delayed_5, Mrem_delayed_6, Mrem_delayed_7, \
                         Mrem_delayed_8, Mrem_delayed_9, Mrem_delayed_10, Mrem_delayed_11, Mrem_delayed_12])

# Metallicity should not be out of this range: [1e-4,1.7e-2]:
Zvalues = np.array([1.0e-4, 2.0e-4, 5.0e-4, 1.0e-3, 2.0e-3, 4.0e-3, 6.0e-3, 8.0e-3, 1.0e-2, 1.4e-2, 1.7e-2, 2.0e-2])

# Mass should not be out of this range: [20, 340] solar masses

Npoints = 100
Mzams = np.linspace(20, 340, Npoints)

# interpolate:
MremInterpol = interpolate.interp2d(Mzams, Zvalues, Mrem_delayed, kind='linear', bounds_error=True)

def Mrem_SEVN(M, Z):
    '''
    Remnant mass as a function of progenitor metallicity and ZAMS mass.

    @in M: ZAMS mass in solar masses ; in range [20, 340]
    @in Z: absolute metallicity      ; in range [1e-4, 1.7e-2]

    @out : remnant mass in solar masses (scalar or array depending on M)
    '''
    
    M_lowerEdge = 60  # absolute lower edge of the upper mass gap (in solar masses)
    M_upperEdge = 120 # absolute upper edge of the upper mass gap (in solar masses)
    
    # check if mass input is an array or not:
    if isinstance(M, np.ndarray): # M is array
        
        out = MremInterpol(M, Z) * (np.heaviside(M_lowerEdge * np.ones(M.size) - MremInterpol(M, Z), 0) \
            + np.heaviside(MremInterpol(M, Z) - M_upperEdge * np.ones(M.size), 0))

    else: # M is not array

        out = MremInterpol(M, Z) * (np.heaviside(M_lowerEdge - MremInterpol(M, Z), 0) \
            + np.heaviside(MremInterpol(M, Z) - M_upperEdge, 0))

    return out

def sample_hardness():
    """
    @out: 3bb hardness parameter
    """
    
    return eta_min * (1 - np.random.rand())**(-2/7)

def f_fb(Mzams):
    """
    Fraction of ejected supernova mass that falls back onto the newly-borned proto-comapct object.
    
    @in Mzams: ZAMS star mass [Msun]
    
    @out: fall-back fraction
    """

    # Chandrasekhar limit in solar masses:
    Mch = 1.4
    
    b_36 = 4.36e-4
    b_37 = 5.22
    b_38 = 6.84e-2
    
    # core mass at the Base of the Asymptotic Giant Branch:
    McBAGB = (b_36*Mzams**b_37 + b_38)**(1/4)
    
    # Carbon/Oxygen core mass:
    M_CO = np.max([Mch,0.773*McBAGB-0.35])
    
    # Determine proto-compact object mass:
    if   M_CO<=3.5:
        
        M_proto = 1.2
        
    elif M_CO>=3.5 and M_CO<6.0:
        
        M_proto = 1.3
        
    elif M_CO>=6.0 and M_CO<11.:
        
        M_proto = 1.4
        
    else:
        
        M_proto = 1.6

    # Determine fall-back fraction:
    a2 = 0.133 - 0.093/(Mzams - M_proto)
    b2 = -11*a2 + 1

    if   M_CO<2.5:

        Mfb = 0.2
        ffb = Mfb/(Mzams-M_proto)

    elif M_CO>=2.5 and M_CO<3.5:

        Mfb = 0.5*M_CO - 1.05
        ffb = Mfb/(Mzams-M_proto)

    elif M_CO>=3.5 and M_CO<11:

        ffb = a2*M_CO + b2

    else:
        
        ffb = 1
        
    return ffb
    
# end of file