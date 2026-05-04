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
from .functions import *

def IMF_kroupa(m):
    """
    Kroupa (2002) initial mass function.

    @in m : stellar mass or stellar mass array [Msun]

    @out: number dN of stars in mass bin (m,m+dm) [1/Msun]
    """

    m = np.asarray(m, dtype=float)
    scalar_input = (m.ndim == 0)
    m = np.atleast_1d(m)

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
                    
    return float(out[0]) if scalar_input else out

N_grid = 700
M_grid = np.linspace(10, 340, N_grid)
Z_grid = np.logspace(np.log10(1e-4), np.log10(2e-2), N_grid)

# Get the directory where constants.py lives (the 'rapster' folder):
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Fryer+2012 remnant mass lookup tables:
DATA_PATH_F12d = os.path.join(BASE_DIR, '..', 'Data', 'MzamsMrem', 'MzamsMrem_F12d.txt')
DATA_PATH_F12r = os.path.join(BASE_DIR, '..', 'Data', 'MzamsMrem', 'MzamsMrem_F12r.txt')

Mremnants_F12d = np.loadtxt(DATA_PATH_F12d, unpack=True)
Mremnants_F12r = np.loadtxt(DATA_PATH_F12r, unpack=True)

MremInterpol_F12d = interpolate.RegularGridInterpolator((M_grid, Z_grid), Mremnants_F12d, method='linear', bounds_error=True)
MremInterpol_F12r = interpolate.RegularGridInterpolator((M_grid, Z_grid), Mremnants_F12r, method='linear', bounds_error=True)

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
        
        out = MremInterpol_F12d((M, Z)) * (np.heaviside(M_lowerEdge * np.ones(M.size) - MremInterpol_F12d((M, Z)), 0) \
            + np.heaviside(MremInterpol_F12d((M, Z)) - M_upperEdge * np.ones(M.size), 0))
    else: # M is not array
        
        out = MremInterpol_F12d((M, Z)) * (np.heaviside(M_lowerEdge - MremInterpol_F12d((M, Z)), 0) \
            + np.heaviside(MremInterpol_F12d((M, Z)) - M_upperEdge, 0))
        out = float(out)
    
    return out
Mrem_F12d = np.vectorize(Mrem_F12d)

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
        
        out = MremInterpol_F12r((M, Z)) * (np.heaviside(M_lowerEdge * np.ones(M.size) - MremInterpol_F12r((M, Z)), 0) \
            + np.heaviside(MremInterpol_F12r((M, Z)) - M_upperEdge * np.ones(M.size), 0))
    else: # M is not array
        
        out = MremInterpol_F12r((M, Z)) * (np.heaviside(M_lowerEdge - MremInterpol_F12r((M, Z)), 0) \
            + np.heaviside(MremInterpol_F12r((M, Z)) - M_upperEdge, 0))
        out = float(out)
    
    return out
Mrem_F12r = np.vectorize(Mrem_F12r)

DATA_DIR = os.path.join(BASE_DIR, '..', 'Data', 'MzamsMrem')

def _load_grid(pattern, key_prefix, numbered_keys=True, n=12):
    """Load n npz files matching pattern and stack into a (Npoints, n) array.

    Args:
        pattern: Filename pattern with {} for the index (e.g., 'MCO{}.npz').
        key_prefix: Key prefix inside each npz file.
        numbered_keys: If True, keys are 'prefix1', 'prefix2', etc.
            If False, all files use the same key 'prefix'.
        n: Number of files to load.
    """
    arrays = []
    for i in range(1, n + 1):
        path = os.path.join(DATA_DIR, pattern.format(i))
        key = f'{key_prefix}{i}' if numbered_keys else key_prefix
        arrays.append(np.load(path)[key])
    return np.array(arrays).T

# SEVN delayed and rapid remnant mass grids (12 metallicity files each):
Mrem_delayed = _load_grid('MzamsMrem{}_delayed.npz', 'Mrem')
Mrem_rapid   = _load_grid('MzamsMrem{}_rapid.npz',   'Mrem')

# Metallicity should not be out of this range: [1e-4, 1.7e-2]:
Zvalues = np.array([1.0e-4, 2.0e-4, 5.0e-4, 1.0e-3, 2.0e-3, 4.0e-3, 6.0e-3, 8.0e-3, 1.0e-2, 1.4e-2, 1.7e-2, 2.0e-2])

# Mass should not be out of this range: [15, 340] solar masses
Npoints = 500
Mzams = np.linspace(15, 340, Npoints)

# interpolate:
MremInterpol_delayed = interpolate.RegularGridInterpolator((Mzams, Zvalues), Mrem_delayed, method='linear', bounds_error=True)
MremInterpol_rapid   = interpolate.RegularGridInterpolator((Mzams, Zvalues), Mrem_rapid  , method='linear', bounds_error=True)

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
        
        out = MremInterpol_delayed((M, Z)) * (np.heaviside(M_lowerEdge * np.ones(M.size) - MremInterpol_delayed((M, Z)), 0) \
            + np.heaviside(MremInterpol_delayed((M, Z)) - M_upperEdge * np.ones(M.size), 0))

    else: # M is not array

        out = MremInterpol_delayed((M, Z)) * (np.heaviside(M_lowerEdge - MremInterpol_delayed((M, Z)), 0) \
            + np.heaviside(MremInterpol_delayed((M, Z)) - M_upperEdge, 0))
        out = float(out)

    return out
Mrem_SEVNdelayed = np.vectorize(Mrem_SEVNdelayed)

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
        
        out = MremInterpol_rapid((M, Z)) * (np.heaviside(M_lowerEdge * np.ones(M.size) - MremInterpol_rapid((M, Z)), 0) \
            + np.heaviside(MremInterpol_rapid((M, Z)) - M_upperEdge * np.ones(M.size), 0))

    else: # M is not array

        out = MremInterpol_rapid((M, Z)) * (np.heaviside(M_lowerEdge - MremInterpol_rapid((M, Z)), 0) \
            + np.heaviside(MremInterpol_rapid((M, Z)) - M_upperEdge, 0))
        out = float(out)

    return out
Mrem_SEVNrapid = np.vectorize(Mrem_SEVNrapid)

def M_CO_SSE(M, Z):
    """
    Carbon-oxygen mass, from Hurley et al. (2000).
    
    @in M: ZAMS star mass [Msun]
    @in Z: absolute metallicity
    
    @out: CO core mass [Msun]
    """
    
    # Chandrasekhar limit in solar masses:
    Mch = M_Chandrasekhar
    
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

# SEVN CO core mass grid (12 metallicity files):
MCO = _load_grid('MCO{}.npz', 'MCO', numbered_keys=False)

# Metallicity should not be out of this range: [1e-4, 1.7e-2]:
Zvalues = np.array([1.0e-4, 2.0e-4, 5.0e-4, 1.0e-3, 2.0e-3, 4.0e-3, 6.0e-3, 8.0e-3, 1.0e-2, 1.4e-2, 1.7e-2, 2.0e-2])

# Mass should not be out of this range: [15, 340] solar masses
Npoints = 500
Mzams = np.linspace(15, 340, Npoints)

# interpolate:
McoInterpol = interpolate.RegularGridInterpolator((Mzams, Zvalues), MCO, method='linear', bounds_error=True)

def M_CO_SEVN(M, Z):
    """
    Carbon-oxygen mass, from Spera & Mapelli (2017).
    
    @in M: ZAMS star mass [Msun]
    @in Z: absolute metallicity
    
    @out: CO core mass [Msun]
    """
    
    M_CO = float(McoInterpol((M, Z)))
    
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
    Returns supernova kick drawn from Maxwellian with parameter wSN_kick * neutron_star_mass / mBH.

    @in mBH: BH mass (Msun), scalar
    @in wSN_kick: neutron star 1D kick parameter (km/s), scalar
    
    @out: SN kick (km/s), scalar
    """
    
    return get_maxwell_sample(np.sqrt(3) * wSN_kick * neutron_star_mass/mBH)

def R_WhiteDwarf(M_wd=white_dwarf_mass, M_Ch=M_Chandrasekhar):
    """
    White-dwarf mass-radius relation; from Nauenberg 1972.
    Returns the radius in parsec.

    M_wd: mass of the white dwarf (solar masses)
    M_Ch: Chandrasekhar mass (solar masses)
    """

    return 7.80e6*(M_wd/M_Ch)**(-1/3)*np.sqrt(1 - (M_wd/M_Ch)**(4/3))/3.086e16

# Sampling a stellar mass from evolving mass function at simulation time t:

m_fine_grid = np.logspace(np.log10(0.08), np.log10(340.0), 10**5)
pdf_values = IMF_kroupa(m_fine_grid)*m_fine_grid**(star_mass_bias_index) # includes IMF and rate-dependent mass factor

# calculate the CDF:
cdf_values = np.cumsum(pdf_values)
cdf_values /= cdf_values[-1] # normalize

# create an "Inverse CDF" function: Maps [0, 1] -> Mass
inv_cdf = interpolate.interp1d(cdf_values, m_fine_grid, bounds_error=False, fill_value=(0.08, 340.0))

def get_star(t, tBH_form, m_min, m_max):
    """
    Returns a main-sequence star (mass/Msun and radius/pc) from an evolving mass function at time t.

    @in t: current simulation time (Myr)
    @in tBH_form: time of BH-subsystem formation (Myr)
    @in m_min: ZAMS minimum mass (Msun)
    @in m_max: ZAMS maximum mass (Msun)
    """
    
    # current mass limit:
    current_m_max = (1e4/t)**(1/2.5) if t > tBH_form else m_max
    
    # find the CDF value corresponding to the current mass limit:
    cutoff_percentile = np.interp(current_m_max, m_fine_grid, cdf_values)
    
    # draw a random number scaled to the available percentile range:
    u = np.random.random() * cutoff_percentile
    
    # map back to mass:
    mstar = inv_cdf(u)
    
    # radius calculation: (Ryu et al. 2020)
    rstar = 0.93 * mstar**0.88 * R_sun # in pc
    return mstar, rstar

# End of file.
