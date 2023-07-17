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

import numpy as np
import argparse
import precession as pre
from scipy.stats import poisson
from scipy.stats import maxwell
import scipy.integrate as integrate
from scipy import interpolate
import astropy.units as u
from astropy.cosmology import Planck18, z_at_value
import pandas as pd
import time
from numba import njit

# Use astrophysical units throughout:

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

# load Planck18 lookup table:
Planck18_lookup_table = np.load('./Planck18_lookup_table.npz')

# interpolate redshift-lookback and lookback-redshift relations:
lookback_interp = interpolate.interp1d(Planck18_lookup_table['z'], Planck18_lookup_table['lookback'])
redshift_interp = interpolate.interp1d(Planck18_lookup_table['t'], Planck18_lookup_table['redshift'])

# Hardening constant:
Hardening_constant = 4/7

# Average number of intermediate BBH states during binary-single interaction:
N_IMS = 20

# Maximum pericenter distance for interaction normalized to semimajor axis:
kp_max = 2.0

# Minimum binary hardness ratio:
eta_min = 5.0

# Probability for 3bb formation:
P_3bb = 0.80

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

# minimum ZAMS mass of massive star that probably results into BH:
mM_min = 20.0

# end of file
