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

from constants import *

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

N_grid = 700
M_grid = np.linspace(10, 340, N_grid)
Z_grid = np.logspace(np.log10(1e-4), np.log10(2e-2), N_grid)

Mremnants_F12d = np.loadtxt('./MzamsMrem/MzamsMrem_F12d.txt', unpack=True)
Mremnants_F12r = np.loadtxt('./MzamsMrem/MzamsMrem_F12r.txt', unpack=True)

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

# Reading ``delayed'' files exported from SEVN code and stored according to metallicity for various ZAMS masses:
MzamsMrem1  = np.load('../Data/MzamsMrem/MzamsMrem1_delayed.npz' ); Mrem_delayed_1  = MzamsMrem1 ['Mrem1' ]
MzamsMrem2  = np.load('../Data/MzamsMrem/MzamsMrem2_delayed.npz' ); Mrem_delayed_2  = MzamsMrem2 ['Mrem2' ]
MzamsMrem3  = np.load('../Data/MzamsMrem/MzamsMrem3_delayed.npz' ); Mrem_delayed_3  = MzamsMrem3 ['Mrem3' ]
MzamsMrem4  = np.load('../Data/MzamsMrem/MzamsMrem4_delayed.npz' ); Mrem_delayed_4  = MzamsMrem4 ['Mrem4' ]
MzamsMrem5  = np.load('../Data/MzamsMrem/MzamsMrem5_delayed.npz' ); Mrem_delayed_5  = MzamsMrem5 ['Mrem5' ]
MzamsMrem6  = np.load('../Data/MzamsMrem/MzamsMrem6_delayed.npz' ); Mrem_delayed_6  = MzamsMrem6 ['Mrem6' ]
MzamsMrem7  = np.load('../Data/MzamsMrem/MzamsMrem7_delayed.npz' ); Mrem_delayed_7  = MzamsMrem7 ['Mrem7' ]
MzamsMrem8  = np.load('../Data/MzamsMrem/MzamsMrem8_delayed.npz' ); Mrem_delayed_8  = MzamsMrem8 ['Mrem8' ]
MzamsMrem9  = np.load('../Data/MzamsMrem/MzamsMrem9_delayed.npz' ); Mrem_delayed_9  = MzamsMrem9 ['Mrem9' ]
MzamsMrem10 = np.load('../Data/MzamsMrem/MzamsMrem10_delayed.npz'); Mrem_delayed_10 = MzamsMrem10['Mrem10']
MzamsMrem11 = np.load('../Data/MzamsMrem/MzamsMrem11_delayed.npz'); Mrem_delayed_11 = MzamsMrem11['Mrem11']
MzamsMrem12 = np.load('../Data/MzamsMrem/MzamsMrem12_delayed.npz'); Mrem_delayed_12 = MzamsMrem12['Mrem12']
# collect remnant masses with various metallicity values in a single array:
Mrem_delayed = np.array([Mrem_delayed_1, Mrem_delayed_2, Mrem_delayed_3, Mrem_delayed_4, Mrem_delayed_5, Mrem_delayed_6, Mrem_delayed_7, \
                         Mrem_delayed_8, Mrem_delayed_9, Mrem_delayed_10, Mrem_delayed_11, Mrem_delayed_12]).T

# Reading ``rapid'' files exported from SEVN code and stored according to metallicity for various ZAMS masses:
MzamsMrem1  = np.load('../Data/MzamsMrem/MzamsMrem1_rapid.npz' ); Mrem_rapid_1  = MzamsMrem1 ['Mrem1' ]
MzamsMrem2  = np.load('../Data/MzamsMrem/MzamsMrem2_rapid.npz' ); Mrem_rapid_2  = MzamsMrem2 ['Mrem2' ]
MzamsMrem3  = np.load('../Data/MzamsMrem/MzamsMrem3_rapid.npz' ); Mrem_rapid_3  = MzamsMrem3 ['Mrem3' ]
MzamsMrem4  = np.load('../Data/MzamsMrem/MzamsMrem4_rapid.npz' ); Mrem_rapid_4  = MzamsMrem4 ['Mrem4' ]
MzamsMrem5  = np.load('../Data/MzamsMrem/MzamsMrem5_rapid.npz' ); Mrem_rapid_5  = MzamsMrem5 ['Mrem5' ]
MzamsMrem6  = np.load('../Data/MzamsMrem/MzamsMrem6_rapid.npz' ); Mrem_rapid_6  = MzamsMrem6 ['Mrem6' ]
MzamsMrem7  = np.load('../Data/MzamsMrem/MzamsMrem7_rapid.npz' ); Mrem_rapid_7  = MzamsMrem7 ['Mrem7' ]
MzamsMrem8  = np.load('../Data/MzamsMrem/MzamsMrem8_rapid.npz' ); Mrem_rapid_8  = MzamsMrem8 ['Mrem8' ]
MzamsMrem9  = np.load('../Data/MzamsMrem/MzamsMrem9_rapid.npz' ); Mrem_rapid_9  = MzamsMrem9 ['Mrem9' ]
MzamsMrem10 = np.load('../Data/MzamsMrem/MzamsMrem10_rapid.npz'); Mrem_rapid_10 = MzamsMrem10['Mrem10']
MzamsMrem11 = np.load('../Data/MzamsMrem/MzamsMrem11_rapid.npz'); Mrem_rapid_11 = MzamsMrem11['Mrem11']
MzamsMrem12 = np.load('../Data/MzamsMrem/MzamsMrem12_rapid.npz'); Mrem_rapid_12 = MzamsMrem12['Mrem12']
# collect remnant masses with various metallicity values in a single array:
Mrem_rapid = np.array([Mrem_rapid_1, Mrem_rapid_2, Mrem_rapid_3, Mrem_rapid_4, Mrem_rapid_5, Mrem_rapid_6, Mrem_rapid_7, \
                         Mrem_rapid_8, Mrem_rapid_9, Mrem_rapid_10, Mrem_rapid_11, Mrem_rapid_12]).T

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

path = '../Data/MzamsMrem/'
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
MCO = np.array([MCO1, MCO2, MCO3, MCO4, MCO5, MCO6, MCO7, MCO8, MCO9, MCO10, MCO11, MCO12]).T

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

end of file.
