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

    '''
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
    '''
    
    m_rem = remnant_mass(m1, m2, chi1, chi2, theta1, theta2, dPhi)
    chi_rem = remnant_spin(m1, m2, chi1, chi2, theta1, theta2, dPhi)
    vGW_kick = remnant_kick(m1, m2, chi1, chi2, theta1, theta2, dPhi)
    
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
    
N_grid = 700
M_grid = np.linspace(10, 340, N_grid)
Z_grid = np.logspace(np.log10(1e-4), np.log10(2e-2), N_grid)

Mremnants_F12d = np.loadtxt('./MzamsMrem/MzamsMrem_F12d.txt', unpack=True)
Mremnants_F12r = np.loadtxt('./MzamsMrem/MzamsMrem_F12r.txt', unpack=True)

Mremnants_F12d = np.transpose(Mremnants_F12d)
Mremnants_F12r = np.transpose(Mremnants_F12r)

MremInterpol_F12d = interpolate.interp2d(M_grid, Z_grid, Mremnants_F12d, kind='linear', bounds_error=True)
MremInterpol_F12r = interpolate.interp2d(M_grid, Z_grid, Mremnants_F12r, kind='linear', bounds_error=True)

def Mrem_F12d(M, Z):
    """
    Fryer et al. (2002) delayed remnant mass prescription model.
    
    @in M: ZAMS mass [Msun]
    @in Z: metallicity

    @out : remnant mass [Msun]
    """

    M_lowerEdge = 45  # absolute lower edge of the upper mass gap (in solar masses)
    M_upperEdge = 120 # absolute upper edge of the upper mass gap (in solar masses)
    
    # check if mass input is an array or not:
    if isinstance(M, np.ndarray): # M is array
        
        out = MremInterpol_F12d(M, Z) * (np.heaviside(M_lowerEdge * np.ones(M.size) - MremInterpol_F12d(M, Z), 0) \
            + np.heaviside(MremInterpol_F12d(M, Z) - M_upperEdge * np.ones(M.size), 0))
    else: # M is not array
        
        out = MremInterpol_F12d(M, Z) * (np.heaviside(M_lowerEdge - MremInterpol_F12d(M, Z), 0) \
            + np.heaviside(MremInterpol_F12d(M, Z) - M_upperEdge, 0))
    
    return out

def Mrem_F12r(M, Z):
    """
    Fryer et al. (2002) rapid remnant mass prescription model.
    
    @in M: ZAMS mass [Msun]
    @in Z: metallicity

    @out : remnant mass [Msun]
    """

    M_lowerEdge = 45  # absolute lower edge of the upper mass gap (in solar masses)
    M_upperEdge = 120 # absolute upper edge of the upper mass gap (in solar masses)
    
    # check if mass input is an array or not:
    if isinstance(M, np.ndarray): # M is array
        
        out = MremInterpol_F12r(M, Z) * (np.heaviside(M_lowerEdge * np.ones(M.size) - MremInterpol_F12r(M, Z), 0) \
            + np.heaviside(MremInterpol_F12r(M, Z) - M_upperEdge * np.ones(M.size), 0))
    else: # M is not array
        
        out = MremInterpol_F12r(M, Z) * (np.heaviside(M_lowerEdge - MremInterpol_F12r(M, Z), 0) \
            + np.heaviside(MremInterpol_F12r(M, Z) - M_upperEdge, 0))
    
    return out

# Reading ``delayed'' files exported from SEVN code and stored according to metallicity for various ZAMS masses:
MzamsMrem1  = np.load('./MzamsMrem/MzamsMrem1_delayed.npz' ); Mrem_delayed_1  = MzamsMrem1 ['Mrem1' ]
MzamsMrem2  = np.load('./MzamsMrem/MzamsMrem2_delayed.npz' ); Mrem_delayed_2  = MzamsMrem2 ['Mrem2' ]
MzamsMrem3  = np.load('./MzamsMrem/MzamsMrem3_delayed.npz' ); Mrem_delayed_3  = MzamsMrem3 ['Mrem3' ]
MzamsMrem4  = np.load('./MzamsMrem/MzamsMrem4_delayed.npz' ); Mrem_delayed_4  = MzamsMrem4 ['Mrem4' ]
MzamsMrem5  = np.load('./MzamsMrem/MzamsMrem5_delayed.npz' ); Mrem_delayed_5  = MzamsMrem5 ['Mrem5' ]
MzamsMrem6  = np.load('./MzamsMrem/MzamsMrem6_delayed.npz' ); Mrem_delayed_6  = MzamsMrem6 ['Mrem6' ]
MzamsMrem7  = np.load('./MzamsMrem/MzamsMrem7_delayed.npz' ); Mrem_delayed_7  = MzamsMrem7 ['Mrem7' ]
MzamsMrem8  = np.load('./MzamsMrem/MzamsMrem8_delayed.npz' ); Mrem_delayed_8  = MzamsMrem8 ['Mrem8' ]
MzamsMrem9  = np.load('./MzamsMrem/MzamsMrem9_delayed.npz' ); Mrem_delayed_9  = MzamsMrem9 ['Mrem9' ]
MzamsMrem10 = np.load('./MzamsMrem/MzamsMrem10_delayed.npz'); Mrem_delayed_10 = MzamsMrem10['Mrem10']
MzamsMrem11 = np.load('./MzamsMrem/MzamsMrem11_delayed.npz'); Mrem_delayed_11 = MzamsMrem11['Mrem11']
MzamsMrem12 = np.load('./MzamsMrem/MzamsMrem12_delayed.npz'); Mrem_delayed_12 = MzamsMrem12['Mrem12']
# collect remnant masses with various metallicity values in a single array:
Mrem_delayed = np.array([Mrem_delayed_1, Mrem_delayed_2, Mrem_delayed_3, Mrem_delayed_4, Mrem_delayed_5, Mrem_delayed_6, Mrem_delayed_7, \
                         Mrem_delayed_8, Mrem_delayed_9, Mrem_delayed_10, Mrem_delayed_11, Mrem_delayed_12])

# Reading ``rapid'' files exported from SEVN code and stored according to metallicity for various ZAMS masses:
MzamsMrem1  = np.load('./MzamsMrem/MzamsMrem1_rapid.npz' ); Mrem_rapid_1  = MzamsMrem1 ['Mrem1' ]
MzamsMrem2  = np.load('./MzamsMrem/MzamsMrem2_rapid.npz' ); Mrem_rapid_2  = MzamsMrem2 ['Mrem2' ]
MzamsMrem3  = np.load('./MzamsMrem/MzamsMrem3_rapid.npz' ); Mrem_rapid_3  = MzamsMrem3 ['Mrem3' ]
MzamsMrem4  = np.load('./MzamsMrem/MzamsMrem4_rapid.npz' ); Mrem_rapid_4  = MzamsMrem4 ['Mrem4' ]
MzamsMrem5  = np.load('./MzamsMrem/MzamsMrem5_rapid.npz' ); Mrem_rapid_5  = MzamsMrem5 ['Mrem5' ]
MzamsMrem6  = np.load('./MzamsMrem/MzamsMrem6_rapid.npz' ); Mrem_rapid_6  = MzamsMrem6 ['Mrem6' ]
MzamsMrem7  = np.load('./MzamsMrem/MzamsMrem7_rapid.npz' ); Mrem_rapid_7  = MzamsMrem7 ['Mrem7' ]
MzamsMrem8  = np.load('./MzamsMrem/MzamsMrem8_rapid.npz' ); Mrem_rapid_8  = MzamsMrem8 ['Mrem8' ]
MzamsMrem9  = np.load('./MzamsMrem/MzamsMrem9_rapid.npz' ); Mrem_rapid_9  = MzamsMrem9 ['Mrem9' ]
MzamsMrem10 = np.load('./MzamsMrem/MzamsMrem10_rapid.npz'); Mrem_rapid_10 = MzamsMrem10['Mrem10']
MzamsMrem11 = np.load('./MzamsMrem/MzamsMrem11_rapid.npz'); Mrem_rapid_11 = MzamsMrem11['Mrem11']
MzamsMrem12 = np.load('./MzamsMrem/MzamsMrem12_rapid.npz'); Mrem_rapid_12 = MzamsMrem12['Mrem12']
# collect remnant masses with various metallicity values in a single array:
Mrem_rapid = np.array([Mrem_rapid_1, Mrem_rapid_2, Mrem_rapid_3, Mrem_rapid_4, Mrem_rapid_5, Mrem_rapid_6, Mrem_rapid_7, \
                         Mrem_rapid_8, Mrem_rapid_9, Mrem_rapid_10, Mrem_rapid_11, Mrem_rapid_12])

# Metallicity should not be out of this range: [1e-4, 1.7e-2]:
Zvalues = np.array([1.0e-4, 2.0e-4, 5.0e-4, 1.0e-3, 2.0e-3, 4.0e-3, 6.0e-3, 8.0e-3, 1.0e-2, 1.4e-2, 1.7e-2, 2.0e-2])

# Mass should not be out of this range: [15, 340] solar masses
Npoints = 500
Mzams = np.linspace(15, 340, Npoints)

# interpolate:
MremInterpol_delayed = interpolate.interp2d(Mzams, Zvalues, Mrem_delayed, kind='linear', bounds_error=True)
MremInterpol_rapid   = interpolate.interp2d(Mzams, Zvalues, Mrem_rapid  , kind='linear', bounds_error=True)

def Mrem_SEVNdelayed(M, Z):
    '''
    Remnant mass as a function of progenitor metallicity and ZAMS mass.
    Delayed SN engine assumed.
    
    @in M: ZAMS mass in solar masses ; in range [15, 340]
    @in Z: absolute metallicity      ; in range [1e-4, 1.7e-2]
    
    @out : remnant mass in solar masses (scalar or array depending on M)
    '''
    
    M_lowerEdge = 55  # absolute lower edge of the upper mass gap (in solar masses)
    M_upperEdge = 120 # absolute upper edge of the upper mass gap (in solar masses)
    
    # check if mass input is an array or not:
    if isinstance(M, np.ndarray): # M is array
        
        out = MremInterpol_delayed(M, Z) * (np.heaviside(M_lowerEdge * np.ones(M.size) - MremInterpol_delayed(M, Z), 0) \
            + np.heaviside(MremInterpol_delayed(M, Z) - M_upperEdge * np.ones(M.size), 0))

    else: # M is not array

        out = MremInterpol_delayed(M, Z) * (np.heaviside(M_lowerEdge - MremInterpol_delayed(M, Z), 0) \
            + np.heaviside(MremInterpol_delayed(M, Z) - M_upperEdge, 0))

    return out

def Mrem_SEVNrapid(M, Z):
    '''
    Remnant mass as a function of progenitor metallicity and ZAMS mass.
    Rapid SN engine assumed.
    
    @in M: ZAMS mass in solar masses ; in range [15, 340]
    @in Z: absolute metallicity      ; in range [1e-4, 1.7e-2]
    
    @out : remnant mass in solar masses (scalar or array depending on M)
    '''
    
    M_lowerEdge = 55  # absolute lower edge of the upper mass gap (in solar masses)
    M_upperEdge = 120 # absolute upper edge of the upper mass gap (in solar masses)
    
    # check if mass input is an array or not:
    if isinstance(M, np.ndarray): # M is array
        
        out = MremInterpol_rapid(M, Z) * (np.heaviside(M_lowerEdge * np.ones(M.size) - MremInterpol_rapid(M, Z), 0) \
            + np.heaviside(MremInterpol_rapid(M, Z) - M_upperEdge * np.ones(M.size), 0))

    else: # M is not array

        out = MremInterpol_rapid(M, Z) * (np.heaviside(M_lowerEdge - MremInterpol_rapid(M, Z), 0) \
            + np.heaviside(MremInterpol_rapid(M, Z) - M_upperEdge, 0))

    return out

def sample_hardness():
    """
    @out: 3bb hardness parameter
    """
    
    return eta_min * (1 - np.random.rand())**(-2/7)

def M_CO_SSE(M, Z):
    """
    Carbon-oxygen mass, from Hurley et al. (2000).
    
    @in M: ZAMS star mass [Msun]
    @in Z: absolute metallicity
    
    @out: CO core mass [Msun]
    """
    
    # Chandrasekhar limit in solar masses:
    Mch = 1.4
    
    zeta = np.log(Z / Z_sun)
    
    bp_36 = (1.445216e-1) + (-6.180219e-2) * zeta + (3.093878e-2) * zeta**2 + (+1.567090e-2) * zeta**3
    bp_37 = (1.304129e+0) + (+1.395919e-1) * zeta + (4.142455e-3) * zeta**2 + (-9.732503e-3) * zeta**3
    bp_38 = (5.114149e-1) + (-1.160850e-2) * zeta + (0.000000e+0) * zeta**2 + (+0.000000e+0) * zeta**3
    
    b_36 = bp_36**4
    b_37 = 4.0 * bp_37
    b_38 = bp_38**4
    
    # core mass at the Base of the Asymptotic Giant Branch:
    McBAGB = (b_36 * M**b_37 + b_38)**(1/4)
    
    # Carbon/Oxygen core mass:
    M_CO = np.max([Mch, 0.773 * McBAGB - 0.35])
    
    return M_CO

path = './MzamsMrem/'
MCO1  = np.load(path + 'MCO1.npz' )['MCO']
MCO2  = np.load(path + 'MCO2.npz' )['MCO']
MCO3  = np.load(path + 'MCO3.npz' )['MCO']
MCO4  = np.load(path + 'MCO4.npz' )['MCO']
MCO5  = np.load(path + 'MCO5.npz' )['MCO']
MCO6  = np.load(path + 'MCO6.npz' )['MCO']
MCO7  = np.load(path + 'MCO7.npz' )['MCO']
MCO8  = np.load(path + 'MCO8.npz' )['MCO']
MCO9  = np.load(path + 'MCO9.npz' )['MCO']
MCO10 = np.load(path + 'MCO10.npz')['MCO']
MCO11 = np.load(path + 'MCO11.npz')['MCO']
MCO12 = np.load(path + 'MCO12.npz')['MCO']
# collect CO core masses with various metallicity values in a single array:
MCO = np.array([MCO1, MCO2, MCO3, MCO4, MCO5, MCO6, MCO7, MCO8, MCO9, MCO10, MCO11, MCO12])

# Metallicity should not be out of this range: [1e-4, 1.7e-2]:
Zvalues = np.array([1.0e-4, 2.0e-4, 5.0e-4, 1.0e-3, 2.0e-3, 4.0e-3, 6.0e-3, 8.0e-3, 1.0e-2, 1.4e-2, 1.7e-2, 2.0e-2])

# Mass should not be out of this range: [15, 340] solar masses
Npoints = 500
Mzams = np.linspace(15, 340, Npoints)

# interpolate:
McoInterpol = interpolate.interp2d(Mzams, Zvalues, MCO, kind='linear', bounds_error=True)

def M_CO_SEVN(M, Z):
    """
    Carbon-oxygen mass, from Spera & Mapelli (2017).
    
    @in M: ZAMS star mass [Msun]
    @in Z: absolute metallicity
    
    @out: CO core mass [Msun]
    """
    
    M_CO = McoInterpol(M, Z)+0
    
    return M_CO

def f_fb_delayed(M, M_CO):
    """
    Fraction of ejected supernova mass that falls back onto the newly-born proto-compact object
    Delayed SN engine assumed.
    
    @in M: ZAMS star mass [Msun]
    @in M_CO: CO core mass [Msun]
    
    @out: fall-back fraction
    """
    
    # Determine proto-compact object mass:
    if   M_CO <= 3.5:
        
        M_proto = 1.2
        
    elif M_CO >= 3.5 and M_CO < 6.0:
        
        M_proto = 1.3
        
    elif M_CO >= 6.0 and M_CO < 11.0:
        
        M_proto = 1.4
        
    else:
        
        M_proto = 1.6

    # Determine fall-back fraction:
    a2 = 0.133 - 0.093 / (M - M_proto)
    b2 = -11 * a2 + 1

    if   M_CO < 2.5:

        Mfb = 0.2
        ffb = Mfb / (M - M_proto)

    elif M_CO >= 2.5 and M_CO < 3.5:

        Mfb = 0.5 * M_CO - 1.05
        ffb = Mfb / (M - M_proto)

    elif M_CO >= 3.5 and M_CO < 11:

        ffb = a2 * M_CO + b2

    else:
        
        ffb = 1.0
        
    return ffb

def f_fb_rapid(M, M_CO):
    """
    Fraction of ejected supernova mass that falls back onto the newly-born proto-compact object
    Rapid SN engine assumed.
    
    @in M: ZAMS star mass [Msun]
    @in M_CO: CO core mass [Msun]
    
    @out: fall-back fraction
    """
    
    # Determine proto-compact object mass:
    M_proto = 1.0

    # Determine fall-back fraction:
    a1 = 0.25 - 1.275 / (M - M_proto)
    b1 = -11 * a1 + 1

    if   M_CO < 2.5:

        Mfb = 0.2
        ffb = Mfb / (M - M_proto)
     
    elif M_CO < 6.0 and M_CO >= 2.5:

        Mfb = 0.286 * M_CO - 0.514
        ffb = Mfb / (M - M_proto)
     
    elif M_CO < 7.0 and M_CO >= 6.0:
     
        ffb = 1.0
     
    elif M_CO < 11.0 and M_CO >= 7.0:
     
        ffb = a1 * M_CO + b1
     
    else:
     
        ffb = 1.0
 
    return ffb

def get_SN_kick(mBH, wSN_kick):
    """
    Returns supernova kick drawn from Maxwellian with parameter wSN_kick * 1.4 / mBH.

    @in mBH: BH mass (Msun), scalar
    @in wSN_kick: neutron star 1D kick parameter (km/s), scalar
    
    @out: SN kick (km/s), scalar
    """
    
    return get_maxwell_sample(np.sqrt(3) * wSN_kick * 1.4/mBH)

def remnant_kick(m1, m2, chi1, chi2, theta1, theta2, dPhi):
    """
    Merger remnant kick velocity in km/s.
    Assumes isotropic spin directions.
    From Gerosa & Kesdsen (2016) and references therein.
    
    Inputs:
    @in m1: primary mass
    @in m2: secondary mass
    @in chi1: dimensionless spin parameter of primary
    @in chi2: dimensionless spin parameter of secondary
    @in theta1: angle between orbital ang. mom. and primary spin vector
    @in theta2: angle between orbital ang. mom. and secondary spin vector
    @in dPhi: angle between spin projections in orbital plane
    """
    
    q = m2 / m1 if m2 < m1 else m1 / m2
    
    eta = q / (1 + q)**2
    
    A, B, H, V11, VA, VB, VC, C2, C3, zeeta = \
    1.2e4, -0.93, 6.9e3, 3677.76, 2481.21, 1792.45, 1506.52, 1140, 2481, 145*np.pi/180
    
    cost1, cost2 = np.cos(theta1), np.cos(theta2)
    sint1, sint2 = np.sqrt(1 - cost1**2), np.sqrt(1 - cost2**2)
    
    chi1_para, chi2_para = chi1 * cost1, chi2 * cost2
    
    chi1_perp, chi2_perp = chi1 * sint1, chi2 * sint2
    
    vm = A * eta**2 * (1 - q) / (1 + q) * (1 + B * eta)
    
    Delta_para = (chi1_para - q * chi2_para               ) / (1 + q)
    Delta_perp = (chi1_perp - q * chi2_perp * np.cos(dPhi)) / (1 + q)
    
    chi_para = (chi1_para + q**2 * chi2_para) / (1 + q)**2
    chi_perp = (chi1_perp + q**2 * chi2_perp * np.cos(dPhi)) / (1 + q)**2
    
    vs_perp = H * eta**2 * Delta_para
    
    vs_para = 16 * eta**2 * (Delta_perp * (V11 + 2 * VA * chi_para + 4 * VB * chi_para**2 + 8 * VC * chi_para**3) \
                             + 2 * chi_perp * Delta_para * (C2 + 2 * C3 * chi_para)) \
    * np.cos(np.random.rand() * 2 * np.pi)
    
    vk = np.sqrt(vm**2 + 2 * vm * vs_perp * np.cos(zeeta) + vs_perp**2 + vs_para**2)
    
    return vk

def remnant_spin(m1, m2, chi1, chi2, theta1, theta2, dPhi):
    """
    Dimensionles merger remnant spin parameter in [0, 1].
    From Gerosa & Kesdsen (2016) and references therein.
    
    Inputs:
    @in m1: primary mass
    @in m2: secondary mass
    @in chi1: dimensionless spin parameter of primary
    @in chi2: dimensionless spin parameter of secondary
    @in theta1: angle between orbital ang. mom. and primary spin vector
    @in theta2: angle between orbital ang. mom. and secondary spin vector
    @in dPhi: angle between spin projections in orbital plane
    """
    
    q = m2 / m1 if m2 < m1 else m1 / m2
    
    eta = q / (1 + q)**2
    
    cost1, cost2 = np.cos(theta1), np.cos(theta2)
    sint1, sint2 = np.sqrt(1 - cost1**2), np.sqrt(1 - cost2**2)
    
    chi1_para, chi2_para = chi1 * cost1, chi2 * cost2
    
    chi1_perp, chi2_perp = chi1 * sint1, chi2 * sint2
    
    Delta_para = (chi1_para - q * chi2_para               ) / (1 + q)
    Delta_perp = (chi1_perp - q * chi2_perp * np.cos(dPhi)) / (1 + q)
    
    chi_para = (chi1_para + q**2 * chi2_para) / (1 + q)**2
    chi_perp = (chi1_perp + q**2 * chi2_perp * np.cos(dPhi)) / (1 + q)**2
    
    chi1_chi2 = chi1_para * chi2_para + chi1_perp * chi2_perp * np.cos(dPhi)
    
    chi_squared = (chi1**2 + q**4 * chi2**2 + 2 * q**2 * chi1_chi2) / (1 + q)**4
    
    t0, t2, t3, s4, s5 = -2.8904, -3.51712, 2.5763, -0.1229, 0.4537
    
    el = 2 * np.sqrt(3) + t2 * eta + t3 * eta**2 + s4 * (1 + q)**4 / (1 + q**2) * chi_squared \
    + (s5 * eta + t0 + 2) * (1 + q)**2 / (1 + q**2) * chi_para
    
    chi_f = np.sqrt(q**2 * el**2 / (1 + q)**4 + chi_squared + 2 * q * el * chi_para / (1 + q)**2)
    
    chi_f = chi_f if chi_f < 1 else 1
    
    return chi_f

def remnant_mass(m1, m2, chi1, chi2, theta1, theta2, dPhi):
    """
    Merger remnant mass.
    From Gerosa & Kesdsen (2016) and references therein.
    
    Inputs:
    @in m1: primary mass
    @in m2: secondary mass
    @in chi1: dimensionless spin parameter of primary
    @in chi2: dimensionless spin parameter of secondary
    @in theta1: angle between orbital ang. mom. and primary spin vector
    @in theta2: angle between orbital ang. mom. and secondary spin vector
    @in dPhi: angle between spin projections in orbital plane
    """
    
    q = m2 / m1 if m2 < m1 else m1 / m2
    
    eta = q / (1 + q)**2
    
    cost1, cost2 = np.cos(theta1), np.cos(theta2)
    sint1, sint2 = np.sqrt(1 - cost1**2), np.sqrt(1 - cost2**2)
    
    chi1_para, chi2_para = chi1 * cost1, chi2 * cost2
    
    chi_para = (chi1_para + q**2 * chi2_para) / (1 + q)**2
    
    Z1 = 1 + (1 - chi_para**2)**(1 / 3) * ((1 + chi_para)**(1 / 3) + (1 - chi_para)**(1 / 3))
    Z2 = np.sqrt(3 * chi_para**2 + Z1**2)
    
    r_isco = 3 + Z2 - np.sign(chi_para) * np.sqrt((3 - Z1) * (3 + Z1 + 2 * Z2))
    
    E_isco = np.sqrt(1 - 2 / 3 / r_isco)
    
    p0, p1 = 0.04827, 0.01707
    
    m_f = (m1 + m2) * (1 - eta * (1 - 4 * eta) * (1 - E_isco) - 16 * eta**2 * (p0 + 4 * p1 * chi_para * (chi_para + 1)))
    
    return m_f

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


# end of file
