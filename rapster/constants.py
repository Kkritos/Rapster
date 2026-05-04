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

import numpy as np
import argparse
from scipy.stats import poisson
from scipy.stats import maxwell
import scipy.integrate as integrate
from scipy import interpolate
import time
import pickle
import pandas as pd
from math import erf
from scipy.optimize import root
from scipy.interpolate import CloughTocher2DInterpolator
import os

# Use astrophysical units throughout:
# mass -> solar masses (M_sun)
# time -> million years (Myr)
# distance -> parsec (pc)
# velocity -> km/s=pc/Myr

# Gravitational constant:
G_Newton = 1 / 232

# speed of light:
c_light = 3.0e5

# Hubble parameter:
h = 0.70

# Hubble time:
t_Hubble = 9.78e3 / h

# Matter density parameter:
Omega_M = 0.30

# Radiation density parameter:
Omega_R = 0.0

# Curvature density parameter:
Omega_K = 0.0

# Vacuum density parameter:
Omega_V = 0.70

# Hardening constant:
Hardening_constant = 4.0 /7

# Average number of intermediate BBH states during binary-single interaction:
N_IMS = 20

# Maximum pericenter distance for interaction normalized to semimajor axis:
kp_max = 2.0

# Minimum binary hardness ratio:
eta_min = 5.0

# Probability for 3bb formation:
P_3bb = 0.8

# core collapse factor:
k_cc = 3.21

# Coulomb logarithm prefactor:
lc = 0.02

# Solar metallicity:
Z_sun = 0.014

# Solar radius:
R_sun = 2.25e-8

# Probability for star ejection from isolated cluster (assuming Maxwellian):
xi_e0 = 0.0074

# Stellar evolution start time:
t_sev = 2.0

# Black hole formation time:
tBH_form = 3.5

# Maximum fraction of cluster's mass locked in BHs:
fBH_max = 0.20

# Stellar evolution parameter:
nu_sev = 0.07

# Burning coefficient:
zeta = 0.08

# BH burning coefficient:
zeta_BH = zeta

# High-mass IMF power index:
alphaIMF = -2.3

# minimum BH mass:
mBH_min = 3.0

# minimum ZAMS mass of a massive star that probably results in BH:
mM_min = 15.0

# Chandrasekhar mass:
M_Chandrasekhar = 1.40

# neutron star mass:
neutron_star_mass = M_Chandrasekhar

# white dwarf mass:
white_dwarf_mass = 0.60

# Bias power-law index for sampling stellar mass for TDEs:
star_mass_bias_index = 2.55
# effective sampling pdf: p(m_star) ~ p_IMF(m_star)*m_star**(star_mass_bias_index)

evolution_keys = [
    'seed', 't', 'z', 'dt', 'm_avg', 'M_cl', 'r_h', 'R_gal', 'v_gal', 't_rh', 't_rhBH', 'n_star', 'N_BH', 'mBH_avg', 'mBH_max', 'r_hBH', 'r_cBH', 'S', 'xi', 'psi', 'psi_BH', 't_3bb', 't_2cap', 'k_3bb', 'k_2cap', 'N_me', 'N_BBH', 'N_meRe', 
    'N_meEj', 'v_star', 'v_BH', 'n_hBH', 'n_cBH', 'n_aBH', 'N_3bb', 'N_2cap', 'N_3cap', 'N_BHej', 'N_BBHej', 'N_dis', 'N_ex', 't_bb', 'N_bb', 'N_meFi', 'N_me2b', 't_ex1', 't_ex2', 'k_ex1', 'k_ex2', 'N_ex1', 'N_ex2', 'N_BHstar', 
    't_pp', 'k_pp', 'N_pp', 'v_esc', 'v_escBH', 'N_triples', 'N_ZLK', 'N_WD', 'v_WD', 'k_tdeBHWD', 'N_tdeBHWD', 'dN_WDformdt', 'dN_WDevdt', 'dN_tdeBHWDdt', 'dN_WDdt', 'dN_tdeBHstardt', 'N_tdeBHstar'
]
hardening_keys = [
    't', 'dt', 't_local', 'dt_local', 'ind', 'a', 'e', 'm1', 'm2', 'q', 'condition', 'N_ex'
]
merger_keys = [
    'seed', 'ind', 'channel', 'a', 'e', 'm1', 'm2', 'chi1', 'chi2', 'g1', 'g2', 'theta1', 'theta2', 
    'dPhi', 'tForm', 'zForm', 't', 'z', 'mRem', 'chiRem', 'gRem', 'vGW', 'chiEff', 'q',  'v_esc',  'h1', 'h2', 
    'Mcl0', 'rh0', 'Z', 'zClForm', 'Rgal0', 'Mcl', 'rh', 'Rgal'
]
tdes_keys = ['seed', 't', 'z', 'type', 'm_star', 'R_star', 'm_BH', 's_BH', 'g_BH', 'r_t', 'r_p', 'beta', 'iota', 'r_mb', 'dm', 's_new', 'v_rel', 'h_BH']

# End of file.
